#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kernel.h"
#include "svm.h"

#include "readsvm.h"
//#define INFINITY	__builtin_inf() // todo: get rid of that



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

double compute_wxi(int p, int index, double *x_weights, double *y_weights) // w \cdot x_i
{
    double res=0.0;
    int i;
    for(i=0;i<prob[0].l;i++)
    {
//        printf("xw = %f   kernel = %f sum = %f \n ", x_weights[i], kernel(0,i,p,index), kernel(0,i,p,index));
        res += x_weights[i] * kernel(0,i,p,index);

    }

    //printf("res1 = %f \n", res);

    for(i=0;i<prob[1].l;i++)
    {
//        printf("yw = %f   kernel = %f sum = %f \n ", y_weights[i], kernel(1,i,p,index), kernel(1,i,p,index));
        res -= y_weights[i] * kernel(1,i,p,index);
    }

    //printf("res2 = %f \n", res);
    return res;
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


//    int max_index_q = compute_max_index(prob[1]); //Todo: besseren namen finden
//    int max_index_p = compute_max_index(prob[0]);

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

    model.label[0] = 1;
    model.label[1] = -1; //todo: richtig machen. (namen aus daen file lesen?)

    // <q, p-q> - <p-q, p-q>/2
    model.rho = rho;

    //printf("rho = %f \n", rho);

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

//    rho = 1.0;
/*
    // compute beta
    int k;
    double beta_average = 0.0;

    for(i=0;i<prob[0].l;i++)
    {
        if(x_weights[i] != 0.0)
        {
            double beta = compute_wxi(0,i, x_weights, y_weights) - 1;
            beta_average += beta / (model.nSV[0] + model.nSV[1]);
            //printf(" beta = %f \n", beta);
        }
    }

    for(i=0;i<prob[1].l;i++)
    {
        if(y_weights[i] != 0.0)
        {
            double beta = compute_wxi(1,i, x_weights, y_weights) + 1;
            beta_average += beta / (model.nSV[0] + model.nSV[1] );
            //printf(" beta = %f \n", beta);
        }
    }

    printf("beta average = %f %d %d\n", beta_average, model.nSV[0], model.nSV[1]);

    double beta = beta_average;

    //classify

    int counter = 0;
    int count_correct = 0;

    for(k=0;k<2;k++)
    for(i=0;i<prob[k].l;i++)
    {
        counter++;
        double res = compute_wxi(k,i,x_weights, y_weights) - beta;
        //printf("res = %f \n", res);
        if(k==0 && res > 0.0) count_correct++;
        if(k==1 && res < 0.0) count_correct++;

    }

    printf("counter = %d, count_correct = %d \n", counter, count_correct);
	*/
    return 0;
}

