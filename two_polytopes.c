#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "readsvm.h"
//#define INFINITY	__builtin_inf() // todo: get rid of that



struct svm_parameter param;		// set by parse_command_line

void exit_with_help()
{
	printf(
	"Usage: svm-train [options] training_set_file [model_file]\n"
	"options:\n"
	"-s svm_type : set type of SVM (default 0)\n"
	"	0 -- C-SVC\n"
	"	1 -- nu-SVC\n"
	"	2 -- one-class SVM\n"
	"	3 -- epsilon-SVR\n"
	"	4 -- nu-SVR\n"
	"-t kernel_type : set type of kernel function (default 2)\n"
	"	0 -- linear: u'*v\n"
	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	"	4 -- precomputed kernel (kernel values in training_set_file)\n"
	"-d degree : set degree in kernel function (default 3)\n"
	"-g gamma : set gamma in kernel function (default 1/k)\n"
	"-r coef0 : set coef0 in kernel function (default 0)\n"
	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-m cachesize : set cache memory size in MB (default 100)\n"
	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
	"-v n: n-fold cross validation mode\n"
	"-q : quiet mode (no outputs)\n"
	);
	exit(1);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name)
{
	int i;

	// default values
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0;	// 1/k
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
//	cross_validation = 0;

	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 's':
				param.svm_type = atoi(argv[i]);
				break;
			case 't':
				param.kernel_type = atoi(argv[i]);
				break;
			case 'd':
				param.degree = atoi(argv[i]);
				break;
			case 'g':
				param.gamma = atof(argv[i]);
				break;
			case 'r':
				param.coef0 = atof(argv[i]);
				break;
			case 'n':
				param.nu = atof(argv[i]);
				break;
			case 'm':
				param.cache_size = atof(argv[i]);
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'e':
				param.eps = atof(argv[i]);
				break;
			case 'p':
				param.p = atof(argv[i]);
				break;
			case 'h':
				param.shrinking = atoi(argv[i]);
				break;
			case 'b':
				param.probability = atoi(argv[i]);
				break;
//			case 'q':
//				svm_print_string = &print_null;
//				i--;
//				break;
//			case 'v':
//				cross_validation = 1;
//				nr_fold = atoi(argv[i]);
//				if(nr_fold < 2)
//				{
//					fprintf(stderr,"n-fold cross validation: n must >= 2\n");
//					exit_with_help();
//				}
//				break;
			case 'w':
				++param.nr_weight;
				param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
				param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
				param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
				param.weight[param.nr_weight-1] = atof(argv[i]);
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}

	// determine filenames

	if(i>=argc)
		exit_with_help();

	strcpy(input_file_name, argv[i]);

	if(i<argc-1)
		strcpy(model_file_name,argv[i+1]);
	else
	{
		char *p = strrchr(argv[i],'/');
		if(p==NULL)
			p = argv[i];
		else
			++p;
		sprintf(model_file_name,"%s.model",p);
	}
}


inline double powi(double base, int times)
{
	double tmp = base, ret = 1.0;

    int t;
	for(t=times; t>0; t/=2)
	{
		if(t%2==1) ret*=tmp;
		tmp = tmp * tmp;
	}
	return ret;
}

int compute_max_index(const struct svm_problem prob)
{
    int i;
    int max_index = 1;
    for (i=0;i<prob.l;i++)
    {
        struct svm_node* px = prob.x[i];
        while (px->index != -1)
        {
            if (px->index > max_index)
                max_index = px->index;
            px++;
        }
    }
    return max_index;
}

void print_vector( const struct svm_node* vector)
{
    while (vector->index != -1)
    {
        printf("%d: %f  ", vector->index, vector->value);
        ++vector;
    }
    printf("\n");
}

double dot(const struct svm_node *px, const struct svm_node *py)
{
//    print_vector(px);
//    print_vector(py);
	double sum = 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			sum += px->value * py->value;
            //printf("  dot: %f = %f * %f, sum = %f, index = %d \n", px->value * py->value, px->value, py->value, sum, px->index);
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

//double kernel(int set1, int element1, int set2, int element2)
//{
//    return dot(prob[set1].x[element1], prob[set2].x[element2]);
//}

double kernel_linear(int set1, int element1, int set2, int element2)
//double kernel(int set1, int element1, int set2, int element2)
{
    return dot(prob[set1].x[element1], prob[set2].x[element2]);
}

double kernel_poly(int set1, int element1, int set2, int element2)
{
    return powi(param.gamma*dot(prob[set1].x[element1], prob[set2].x[element2])+param.coef0,param.degree);
}

double kernel_rbf(int set1, int element1, int set2, int element2)
{
//    return exp(-param.gamma*(x_square[i]+x_square[j]-2*dot(prob[set1].x[element1], prob[set2].x[element2])));

//    return exp(-param.gamma*( dot(prob[set1].x[element1], prob[set1].x[element1])+
//                        dot(prob[set1].x[element1], prob[set2].x[element2])-2*
//                        dot(prob[set1].x[element1], prob[set2].x[element2])));

    double dots = ( dot(prob[set1].x[element1], prob[set1].x[element1])+
                        dot(prob[set1].x[element1], prob[set2].x[element2])-2*
                        dot(prob[set1].x[element1], prob[set2].x[element2]));
    double wgamma = -param.gamma*dots;
    double wexp = exp(wgamma);

    return wexp;

}

double kernel_sigmoid(int set1, int element1, int set2, int element2)
{
    return tanh(param.gamma*dot(prob[set1].x[element1], prob[set2].x[element2])+param.coef0);
}

/*double kernel_precomputed(int set1, int element1, int set2, int element2)
{
    return x[i][(int)(x[j][0].value)].value;
}*/

double (*kernel)(int, int, int, int);


void add_to_weights(double* weights, double lambda, int index, int set)
{
    int i;
    for (i=0;i<prob[set].l;i++)
    {
        if (i!=index)
            weights[i] *= lambda;
        else
            weights[i] = weights[i]*lambda + (1.0 - lambda)*1;
    }
}

void print_weights(double* weights, int set)
{
	int i;
	for(i=0;i<prob[set].l;i++)
	{
			printf("%d:\t%f  ",i,weights[i]);
	}
	printf("\n");
}

int find_max(int p, double *dot_yi_x, double* dot_xi_x, double dot_xi_yi, double dot_xi_xi, double *max_p) {
    // find max
    int max_p_index = -1;
    *max_p = -HUGE_VAL;
    int i;
    for (i=0;i<prob[p].l;i++) {
        double sum = dot_yi_x[i] - dot_xi_x[i] - dot_xi_yi + dot_xi_xi;
        if(sum > *max_p)
        {
            *max_p = sum;
            max_p_index = i;
        }
    }
    return max_p_index;
}

double compute_zaehler(double dot_xi_yi, double* dot_yi_x, double* dot_xi_x, int p, int max_p_index ) {
    double zaehler = dot_xi_yi - dot_yi_x[max_p_index] - dot_xi_x[max_p_index] + kernel(p,max_p_index, p, max_p_index); //todo: samevector
    return zaehler;
}

double compute_nenner(double dot_xi_xi, double* dot_xi_x, int p, int max_p_index) {
    double nenner = dot_xi_xi - 2* dot_xi_x[max_p_index] +  kernel(p, max_p_index, p, max_p_index);
    return nenner;
}

double update_xi_xi(double dot_xi_xi, double* dot_xi_x, int p, int max_p_index, double lambda) {
    dot_xi_xi = lambda * lambda * dot_xi_xi
            + 2 * lambda * (1.0 - lambda) * dot_xi_x[max_p_index]
            + (1.0 - lambda)*(1.0 - lambda)*kernel(p, max_p_index, p ,max_p_index );
    return dot_xi_xi;
}

double update_xi_yi(double dot_xi_yi, double* dot_yi_x, int max_p_index, double lambda) {
    dot_xi_yi = lambda * dot_xi_yi + (1.0 - lambda) * dot_yi_x[max_p_index];
    return dot_xi_yi;
}

void update_xi_x(double* dot_xi_x, int p, int p2, int max_p_index, double lambda) {
    int i;
    for (i=0;i<prob[p2].l;i++) {
        dot_xi_x[i]= dot_xi_x[i] * lambda + (1.0 - lambda) * kernel(p, max_p_index, p2, i);
    }
}

void compute_weights(double *x_weights, double* y_weights)
{
    printf("Gewichtsvektoren initialisieren.. \n");

    // initialize weights
    int i;
    for (i=0;i<prob[0].l;i++)
        x_weights[i] = 0.0;

    for (i=0;i<prob[1].l;i++)
        y_weights[i] = 0.0;

    x_weights[0] = 1.0;
    y_weights[0] = 1.0;


    // deklaration der variablen die werte zwischenspeichern
    double *dot_xi_x; // < x_i, x> \forall x \in P
    double *dot_yi_x;  // < y_i, x> \forall x \in P
    double dot_xi_yi; // <x_i, y_i >
    double dot_xi_xi; // <x_i, x_i >

    double *dot_yi_y; // < y_i, y> \forall y \in Q
    double *dot_xi_y;  // < x_i, y> \forall y \in Q
    double dot_yi_yi; // <y_i, y_i >

    // speicher anfordern
    dot_xi_x = (double *) malloc(prob[0].l * sizeof(double));
    dot_yi_x = (double *) malloc(prob[0].l * sizeof(double));

    dot_xi_y = (double *) malloc(prob[1].l * sizeof(double));
    dot_yi_y = (double *) malloc(prob[1].l * sizeof(double));

    // initialisieren
    for (i=0;i<prob[0].l;i++) {
        dot_xi_x[i]=kernel(0, 0, 0, i);
        dot_yi_x[i]=kernel(1, 0, 0, i);
    }

    for (i=0;i<prob[1].l;i++) {
        dot_xi_y[i]=kernel(0, 0, 1, i);
        dot_yi_y[i]=kernel(1, 0, 1, i);
    }

    dot_xi_xi = kernel(0, 0, 0, 0);
    dot_xi_yi = kernel(0, 0, 1, 0);
    dot_yi_yi = kernel(1, 0, 1, 0);



    // find max
    int max_p_index;
    double max_p;
    max_p_index = find_max(0, dot_yi_x, dot_xi_x, dot_xi_yi, dot_xi_xi, &max_p);

    int max_q_index;
    double max_q;
    max_q_index = find_max(1, dot_xi_y, dot_yi_y, dot_xi_yi, dot_yi_yi, &max_q);

    int j;

    for (j=0;j<10  ;j++)
    {
        double lambda;
        if (max_p >= max_q)
        {
            double zaehler = compute_zaehler(dot_xi_yi, dot_yi_x, dot_xi_x, 0, max_p_index);
            double nenner = compute_nenner(dot_xi_xi, dot_xi_x, 0, max_p_index);

            lambda = zaehler / nenner;

            if(zaehler == 0.0 && nenner == 0.0) lambda = 0.0;
            if(lambda < 0.0)	lambda = 0.0;
            if(lambda > 1.0)	lambda = 0.0;

            add_to_weights(x_weights, lambda, max_p_index, 0);

            // update dotproducts

            dot_xi_xi = update_xi_xi(dot_xi_xi, dot_xi_x, 0, max_p_index, lambda);

            dot_xi_yi = update_xi_yi(dot_xi_yi, dot_yi_x, max_p_index, lambda);

            update_xi_x(dot_xi_x, 0, 0, max_p_index, lambda);

            update_xi_x(dot_xi_y, 0, 1, max_p_index, lambda);
        }
        else
        {
            double zaehler = compute_zaehler(dot_xi_yi, dot_xi_y, dot_yi_y, 1, max_q_index);
            double nenner = compute_nenner(dot_yi_yi, dot_yi_y, 1, max_q_index);

            lambda = zaehler / nenner;

            if(zaehler == 0.0 && nenner == 0.0) lambda = 0.0;
            if(lambda < 0.0)	lambda = 0.0;
            if(lambda > 1.0)	lambda = 0.0;

            add_to_weights(y_weights, lambda, max_q_index, 1);

            // update dotproducts

            dot_yi_yi = update_xi_xi(dot_yi_yi, dot_yi_y, 1, max_q_index, lambda);

            dot_xi_yi = update_xi_yi(dot_xi_yi, dot_xi_y, max_q_index, lambda);

            update_xi_x(dot_yi_y, 1, 1, max_q_index, lambda);

            update_xi_x(dot_yi_x, 1, 0, max_q_index, lambda);
        }
        // find max
        max_q_index = find_max(1, dot_xi_y, dot_yi_y, dot_xi_yi, dot_yi_yi, &max_q);
        max_p_index = find_max(0, dot_yi_x, dot_xi_x, dot_xi_yi, dot_xi_xi, &max_p);


        //duality gap
        // absolute duality gap

        double adg = max_p + max_q;

        //printf("max_p = %f  max_q = %f ", max_p, max_q);
        //printf("adg = %f ", adg);

        // relative duality gap
        // adg / ||p-q||^2 - adg
        // adg / <p-q, p-q> - adg

        double distance = dot_xi_xi + dot_yi_yi - 2 * dot_xi_yi;

		printf("<x-y,x-y> = %e " , distance);
		printf("adg = %e " , adg);

        double rdg_nenner = distance - adg;
        double rdg;

        if (rdg_nenner <= 0)
        {
            //printf("set huge value... ");
            rdg = HUGE_VAL;
        }
        else
        {
            rdg = adg / rdg_nenner;
        }

        printf("rdg = %e \n", rdg);
		//print_weights(x_weights, prob[0]);
		//print_weights(y_weights, prob[1]);

    }
}

int main (int argc, char ** argv)
{

    kernel = &kernel_linear;

    char input_filename[1024];
    char model_filename[1024];

    parse_command_line(argc, argv, input_filename, model_filename);

    printf("gamma = %e  max_index = %d \n", param.gamma, max_index);
//    if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

    param.gamma = 1.0 / 3.0;

    printf("gamma = %e \n", param.gamma);

	switch(param.kernel_type)
	{
		case LINEAR:
			kernel = &kernel_linear;
			break;
		case POLY:
			kernel = &kernel_poly;
			break;
		case RBF:
			kernel = &kernel_rbf;
			break;
		case SIGMOID:
			kernel = &kernel_sigmoid;
			break;
//		case PRECOMPUTED:
//			kernel = &kernel_precomputed;
//			break;
	}


    printf("Dateien einlesen... \n");
//    if (argc < 2)
        //filename = "data/heart_scale.plus";
//        filename = "data/triangle.txt";
//    else
//        filename = argv[1];

    read_problem(input_filename);

    printf("Anzahl der eingelesenen Datensaetze: P = %d   Q = %d \n", prob[0].l, prob[1].l);

    int max_index_q = compute_max_index(prob[1]); //Todo: besseren namen finden
    int max_index_p = compute_max_index(prob[0]);

    double *x_weights, *y_weights;
    x_weights = (double *) malloc(prob[0].l * sizeof(double));
    y_weights = (double *) malloc(prob[1].l * sizeof(double));

    compute_weights(x_weights, y_weights);

    struct svm_model model;
    model.param = param;

    model.weights[0] = x_weights;
    model.weights[1] = y_weights;

    model.SV[0] = prob[0].x;
    model.SV[1] = prob[1].x;


    model.l=0;
    model.nSV[0]=0;
    model.nSV[1]=0;

    int i;
    int j;
    for(j=0;j<2;j++)
    for(i=0;i<prob[j].l;i++)
    {
        if(model.weights[j][i] != 0.0)
        {
            model.l++;
            model.nSV[j]++;
        }
    }

    svm_save_model(model_filename, &model);

    return 0;
}
