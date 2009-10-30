#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

#define INFINITY	__builtin_inf() // todo: get rid of that

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

struct svm_problem prob_p;
struct svm_problem prob_q;

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
			break;
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

int compute_max_index(const struct svm_problem prob)
{
	int i;
	int max_index = 1;
	for(i=0;i<prob.l;i++)
	{
		struct svm_node* px = prob.x[i]; 
		while(px->index != -1)
		{
			if(px->index > max_index)
				max_index = px->index;
			px++;
		}
	}
	return max_index;
}

void init_nodes(struct svm_node* vector, int elements)
{
	int i;
	for(i=0;i< elements-1 ; i++)
		vector[i].index = i+1;
	vector[elements-1].index = -1;
}

void compute_vector_from_weights(double* weights, struct svm_node* vector, const struct svm_problem prob)
{
	int i;
	for(i=0;i<prob.l;i++)
	{
		struct svm_node* out= vector;
		struct svm_node* in = prob.x[0];
		while(in->index != -1 && out->index != -1)
		{
			if(in->index == out->index)
			{
				if(i==0)
					out->value = 0; //todo: schoener machen
				//printf("gleiche Indizes: %d  %d  , value = %f, weight = %f, result = %f outvalue = %f\n", in->index, out->index, in->value, weights[i], in->value*weights[0], out->value);
				out->value += in->value * weights[i];
				//printf("geschriebener Wert = %f \n", out->value);
				
				++in;
				++out;
			}
			else
			{
				//printf("missing element: %d  %d  \n", in->index, out->index);
				if(in->index > out->index)
					++out;
				else
					++in;
			}			
		}
	}	
}

double compute_dot_product_with_differences(const struct svm_node* vector, const struct svm_node* x, const struct svm_node* y) // <p - y_i, x_i - y_i>
{
	double sum = 0;
	while(x->index != -1 && vector->index != -1)
	{
		if(x->index == vector->index)
		{
			//printf("sum = %f  vector = %f  x = %f  y = %f \n", sum, vector->value, x->value, y->value);
			sum += (vector->value - y->value)*(x->value - y->value);
			++x;
			++y;
			++vector;
		}
		else
		{
			//printf("missing element 2.. \n");
			if(x->index > vector->index)
				++vector;
			else {
				++x; ++y; }
		}			
	}
	return sum;
}



int find_max_dotproduct(const struct svm_node* x, const struct svm_node* y, struct svm_problem prob, double &max)
{
	int i;
	max = -HUGE_VAL; //compute_dot_product_with_differences(prob.x[0], x, y);
	int max_index = -1;
	for(i=0;i<prob.l;i++)
//	i=1;
	{
		double dotproduct = compute_dot_product_with_differences(prob.x[i], x, y);
		printf("index = %d dotproduct = %f  max = %f \n", i, dotproduct, max);
		if(dotproduct > max) {
			max = dotproduct;
			max_index = i;
		}
	}
	return max_index;
}

void print_vector( const struct svm_node* vector)
{
	while(vector->index != -1)
	{
		printf("%d: %f \n", vector->index, vector->value);
		++vector;
	}
}

double compute_dot_product_with_differences(const struct svm_node* a, const struct svm_node* b, const struct svm_node* c, const struct svm_node* d) // <a - b, c - d>
{
	double sum = 0;
	while(a->index != -1 && b->index != -1 && c->index != -1 && d->index != -1)
	{
		if(a->index == b->index == c->index == d->index == )
		{
			//printf("sum = %f  vector = %f  x = %f  y = %f \n", sum, vector->value, x->value, y->value);
			sum += (a->value - b->value)*(c->value - d->value);
			++a;
			++b;
			++c;
			++d;
		}
		else
		{
			//printf("missing element 2.. \n"); // todo: was tun bei fehlenden werten?

			struct svm_node* min = a;
			if(min->index > b->index)	min=b;
			if(min->index > c->index)	min=c;
			if(min->index > d->index)	min=d;
		
			++min;

			//if(x->index > vector->index)
			//	++vector;
			//else {
			//	++x; ++y; }
		}			
	}
	return sum;
}

void add_proportional(const struct svm_node* a, const struct svm_node* b, double lambda)
{
	while(a->index != -1 && b->index != -1) {
		if(a->index == b->index)
			a->value = lambda*a->value + (1.0 - lambda)* b->value;
			++a;
			++b;
		} else {
			if(a->index > b->index)
				++b;
			else {
				++a; ++b; }
		}
	}
}

void add_to_weights(double* weights, double lambda, int index, struct svm_problem prob)
{
	int i;
	for(i=0;i<prob.l;i++)
	{
		if(i!=index)
			weights[i] *= lambda;
		else
			weights[i] = weights[i]*lambda + (1.0 - lambda)*1;
	}
}

int main (int argc, const char ** argv)
{  
    const char* filename;
    
    if(argc < 2)
		filename = "data/heart_scale.plus";
    else
    	filename = argv[1];

    read_problem(filename);
    
    prob_p = prob;

    if(argc < 3)
		filename = "data/heart_scale.minus";
    else
    	filename = argv[2];

    read_problem(filename);

    prob_q = prob;
    
    printf("prob_p.l = %d   prob_q.l = %d \n", prob_p.l, prob_q.l);
    
    int max_index_q = compute_max_index(prob_q);
    int max_index_p = compute_max_index(prob_p);
    
    printf("max_index P = %d \n", max_index_p);
    printf("max_index Q = %d \n", max_index_q);

	int elements = max_index_p;
	if(elements < max_index_q)	elements = max_index_q;
    elements++; //add one for the stopping index -1
    
    double* x_weights = (double *) malloc(prob_p.l * sizeof(double));
    struct svm_node* x = (struct svm_node *) malloc(elements * sizeof(struct svm_node));
	init_nodes(x, elements);
    
    double* y_weights = (double *) malloc(prob_q.l * sizeof(double));
    struct svm_node* y = (struct svm_node *) malloc(elements * sizeof(struct svm_node));
	init_nodes(y, elements);

	int i;
	x_weights[0] = 1.0;
	for(i=1;i<prob_p.l;i++)
		x_weights[i]=0.0;

	y_weights[0] = 1.0;
	for(i=1;i<prob_q.l;i++)
		y_weights[i]=0.0;
		
	

	compute_vector_from_weights(x_weights, x, prob_p);

//	print_vector(x);
//	print_vector(prob_p.x[0]);
	

	compute_vector_from_weights(y_weights, y, prob_q);


	// Step 1
	printf("erster Aufruf... \n");
	double max_p;
	int max_p_index = find_max_dotproduct( x, y, prob_p, max_p);

	printf("zweiter Aufruf... \n");
	double max_q;
	int max_q_index = find_max_dotproduct( y, x, prob_q, max_q);

	printf("max_p = %d   max_q = %d \n", max_p, max_q);

	double lambda;

	if(max_p > max_q) {
		double zaehler = compute_dot_product_with_differences(y, prob_p.x[max_p_index], x, prob_p.x[max_p_index]);
		double nenner = compute_dot_product_with_differences(x, prob_p.x[max_p_index], x, prob_p.x[max_p_index]);

		lambda = zaehler / nenner;		

		add_proportional(x,prob_p.x[max_p_index], lambda);
		add_to_weights(x_weights, lambda, max_p_index, prob_p);

		max_p_index = find_max_dotproduct( x, y, prob_p, max_p); // max_p updaten
	} else
	{
		// Gilbertschritt im Polytop Q
		// \lambda = <p-q_i, q-q_i> / <q-q_i, q-q_i>

		// <x - max_q_index, y - max_q_index> /  <y - max_q_index, y - max_q_index>

		double zaehler = compute_dot_product_with_differences(x, prob_q.x[max_q_index], y, prob_q.x[max_q_index]);
		double nenner = compute_dot_product_with_differences(y, prob_q.x[max_q_index], y, prob_q.x[max_q_index]);

		lambda = zaehler / nenner;
		add_proportional(y,prob_q.y[max_q_index], lambda);
		add_to_weights(y_weights, lambda, max_q_index, prob_q);
		
		max_q_index = find_max_dotproduct( y, x, prob_q, max_q); // max_q updaten
	}

	//duality gap
	// absolute duality gap
	// berechne max_p und max_q neu, dann adg = max_p+max_q

	double adg = max_p + max_q;

	// relative duality gap
	// adg / ||p-q||^2 - adg
	// adg / <p-q, p-q> - adg

	double rdg_nenner = compute_dot_product_with_differences(x,y,x,y) - adg;
	double rdg;

	if(rdg_nenner <= 0)
	{
		rdg = HUGE_VAL;
	} else
	{
		rdg = adg / rdg_nenner;
	}

	printf("relative duality gap = %f \n", rdg);

    return 0;
}
