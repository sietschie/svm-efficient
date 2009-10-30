#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

#define INFINITY	__builtin_inf() // todo: get rid of that
//#defince INFINITY	HUGE_VAL //aus libsvm

struct svm_node
{
	int index;
	double value;
};

struct svm_problem
{
	int l;
	double *y;
	struct svm_node **x;
};

struct svm_problem prob;		// set by read_problem

struct svm_problem prob_a;
struct svm_problem prob_b;

static char *line = NULL;
static int max_line_len;
struct svm_node *x_space;
static double C = 5.5;



void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	
	exit(1);
}

static char* readline(FILE *input)
{
	int len;
	
	if(fgets(line,max_line_len,input) == NULL)
		return NULL;
	
	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			brea	;
	}
	return line;
}

void read_problem(const char *filename)
{
	int elements, max_index, inst_max_index, i, j;
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;
	
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}
	
	prob.l = 0;
	elements = 0;
	
	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	while(readline(fp)!=NULL)
	{
		char *p = strtok(line," \t"); // label
		
		// features
		while(1)
		{
			p = strtok(NULL," \t");
			if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			++elements;
		}
		++elements;
		++prob.l;
	}
	rewind(fp);
	
	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node,elements);
	
	max_index = 0;
	j=0;
	for(i=0;i<prob.l;i++)
	{
		inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
		readline(fp);
		prob.x[i] = &x_space[j];
		label = strtok(line," \t");
		prob.y[i] = strtod(label,&endptr);
		if(endptr == label)
			exit_input_error(i+1);
		
		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");
			
			if(val == NULL)
				break;
			
			errno = 0;
			x_space[j].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
				exit_input_error(i+1);
			else
				inst_max_index = x_space[j].index;
			
			errno = 0;
			x_space[j].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(i+1);
			
			++j;
		}
		
		if(inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[j++].index = -1;
	}

	//if(param.gamma == 0 && max_index > 0)
	//	param.gamma = 1.0/max_index;
	
	/*if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if (prob.x[i][0].index != 0)
			{
				fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
				exit(1);
			}
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
				exit(1);
			}
		}*/
	
	fclose(fp);
}

double k(struct svm_node* px, struct svm_node* py)
{
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			sum += px->value * py->value;
			++px;
			++py;
		}
		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}			
	}
	return sum;
}

void generate_K(double *K)
{
	
	int i, j;
	for(i = 0; i < prob.l; i++)
	for(j = 0; j < prob.l; j++)
	{
		K[i + prob.l * j] = prob.y[i] * prob.y[j] * k(prob.x[i], prob.x[j]) + prob.y[i] * prob.y[j];
		if(i == j)
			K[i + prob.l * j] += 1/C;
		//printf(" K[%d + %d * %d] = %f \n", i, prob.l, j, K[i + prob.l * j]);
	}
}

double f(double* K, double *alpha)
{
	double sum = 0;
		
	int i, j;
	for(i = 0; i < prob.l; i++)
	for(j = 0; j < prob.l; j++)
	{
		sum += - alpha[i] * alpha[j] * K[i + prob.l * j];
		//printf("sum = %f, alpha[i] = %f, alpha[j] = %f, K[i + prob.l * j] = %f, i=%d, j=%d \n", sum, alpha[i], alpha[j], K[i + prob.l * j], i, j);
	}
	return sum;
}

double zeile_aus_matrix_mal_vektor(double *K, double *alpha, int index)
{
	double sum = 0;
	int i;
	for(i=0;i<prob.l;i++)
	{
		sum += K[i + prob.l * index] * alpha[i];
	}
}

int argmax_nabla_f(double* K, double* alpha)
{
	int i = 0;
	double max = -2 * zeile_aus_matrix_mal_vektor(K, alpha, 0);
	int max_index = 0;
	
	for(i=1; i < prob.l; i++)
	{
		double current = -2 * zeile_aus_matrix_mal_vektor(K, alpha, i);
		if(current > max)
		{
			max = current;
			max_index = i;
		}
	}
	//printf(" max = %f  max_index = %d \n", max, max_index);
	return max_index;
}

int argmax_f_ei(double *K)
{
	int i;
	int min_index = 0;
	double min_value = K[0];
	for (i=1; i<prob.l; i++) {
		if(K[i + prob.l * i] < min_value)
		{
			min_index = i;
			min_value = K[i + prob.l * i];
		}
	}
	
	return min_index;
}

void generate_einheitsvektor(int index, double *alpha)
{
	int i;
	for(i=0;i<prob.l;i++)
		if(i == index)
			alpha[i]=1.0;
		else
			alpha[i]=0.0;
}

void compute_new_alpha(double* alpha, double t, int is, double *new_alpha)
{
	//printf("new_alpha = ");
	int i;
	for(i=0;i<prob.l;i++)
	{
		if(i==is)
			new_alpha[i] = alpha[i] + t * (1 - alpha[i]);
		else
			new_alpha[i] = alpha[i] - t * ( alpha[i] );
		//printf(" %f ", new_alpha[i]);
	}
	//printf("\n");
}



double find_t(double *K, double *alpha, int is, double * new_alpha)
{

	
	
	double t = 0.0;

	compute_new_alpha(alpha, t, is, new_alpha);	
	double max_value = f(K, new_alpha);
	double max_t = 0.0;

	for (; t<=1.0; t+=0.0001) {
		compute_new_alpha(alpha, t, is, new_alpha);	
		double new_value = f(K, new_alpha);
		//printf("t = %f, value = %f \n", t, new_value);
		if (new_value > max_value) {
			max_value = new_value;
			max_t = t;
		}
	}
		
	printf(" empirical max_t = %f \n", max_t);
	
	return max_t;
}

double e(int i, int j)
{	
	if(i==j)
		return 1.0;
	else 
		return 0.0;
}

double find_t2(double *K, double *alpha, int is)
{
	int i,j;
	
	double zaehler = 0.0;
	double nenner = 0.0;
	
	for (i=0; i<prob.l; i++)
	for (j=0; j<prob.l; j++) {

		zaehler += ( -alpha[i]*e(is, j) + alpha[i]*alpha[j] - alpha[j]*e(is,i) + alpha[i]*alpha[j]) * K[i + prob.l * j];
		nenner += 2 * ( e(is,i)* e(is,j) - alpha[j] * e(is,i) - alpha[i] * e(is,j) + alpha[i] * alpha[j]) * K[i + prob.l * j];
	}
	
	double t = zaehler / nenner;

	//printf("berechnetes t = %f \n" , t );
	
	return t;
}

double compute_absolute_duality_gap(double *K, double *alpha)
{
	int i = argmax_nabla_f(K, alpha);
	double res =  -2*zeile_aus_matrix_mal_vektor(K, alpha, i) - 2 * f(K, alpha);
	return res;
}

double compute_relative_duality_gap(double *K, double *alpha)
{
	double absolute_gap = compute_absolute_duality_gap( K, alpha);

	double nenner = absolute_gap + f(K,alpha); 
	double res;

	if(nenner >= 0.0)
		res = INFINITY; 
	else
		res = - absolute_gap / (absolute_gap + f(K,alpha));
	
	return res;
}

int main (int argc, const char ** argv)
{  
    const char* filename;
    
    if(argc < 2)
		filename = "data/heart_scale.plus";
    else
    	filename = argv[1];

    read_problem(filename);
    
    prob_a = prob;

    if(argc < 3)
		filename = "data/heart_scale.minus";
    else
    	filename = argv[2];

    read_problem(filename);

    prob_b = prob;
    
    printf("prob_a.l = %d   prob_b.l = %d \n", prob_a.l, prob_b.l);

/* 	double *K = (double*) malloc(prob.l * prob.l * sizeof(double));
	//double *alpha = malloc(prob.l * sizeof(double));
   
    generate_K(K);
	
	int start_index = argmax_f_ei(K);
	
	printf("start_index = %d \n", start_index);

	double* alpha = (double*) malloc(prob.l * sizeof(double));
	double *new_alpha = (double*) malloc(prob.l * sizeof(double));

	generate_einheitsvektor(start_index, alpha);
//    double test = f()
	
	double min_error= 5;
	double rdg = INFINITY; //todo:besseren initialen wert finden
	int counter = 0;
	
	//while (rdg < min_error) {
	while (rdg > min_error) { // todo: ist kleiner als ein schreibfehler?
	
		int is = argmax_nabla_f(K, alpha);
	
		//printf(" is = %d\n", is);

		double t = find_t2(K, alpha, is);
		
		compute_new_alpha(alpha, t, is, new_alpha);	
		
		//t = find_t(K, alpha, is, new_alpha);
	
		double *temp = alpha;
		alpha = new_alpha;
		new_alpha = temp;
	
		double adg = compute_absolute_duality_gap(K, alpha);
		rdg = compute_relative_duality_gap(K,alpha);
	
		counter++;
		
		//printf("adg = %f, rdg = %f f(alpha) = %f is = %d count = %d\n", adg, rdg, f(K, alpha), is, counter);
	
	}

	int i;
	for(i = 0; i < prob.l; i++)
	{
		//printf(" alpha[%3d} = %f \n",i, alpha[i]);
	}*/

    return 0;
}
