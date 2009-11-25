#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "readsvm.h"
//#define INFINITY	__builtin_inf() // todo: get rid of that


struct svm_problem prob_p;
struct svm_problem prob_q;

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

int find_max_dotproduct_weights(double* weights_x, double* weights_y, struct svm_problem prob_x, struct svm_problem prob_y, double *max_ret)
{
    int i;
    double max = -HUGE_VAL; //compute_dot_product_with_differences(prob.x[0], x, y);
    int max_index = -1;
    for (i=0;i<prob.l;i++)
    {
        double sum = 0.0;
        int j;
        for(j=0;j<prob_y.l;j++)
        {
            sum += weights_y[j] * dot(prob_x.x[i], prob_y.x[j]);
        }
        for(j=0;j<prob_x.l;j++)
        {
            sum -= weights_x[j] * dot(prob_x.x[i], prob_x.x[j]);
        }
        for(j=0;j<prob_x.l;j++)
        {
            int k;
            for(k=0;k<prob_y.l;k++)
            {
                sum -= weights_x[j] * weights_y[k] * dot(prob_x.x[j], prob_y.x[k]);
            }
            for(k=0;k<prob_x.l;k++)
            {
                sum += weights_x[j] * weights_x[k] * dot(prob_x.x[j], prob_x.x[k]);
            }
        }

        if( max < sum) {
            max = sum;
            max_index = i;
        }
    }
    (*max_ret) = max; //todo: watch closely
    return max_index;
}

void add_to_weights(double* weights, double lambda, int index, struct svm_problem prob)
{
    int i;
    for (i=0;i<prob.l;i++)
    {
        if (i!=index)
            weights[i] *= lambda;
        else
            weights[i] = weights[i]*lambda + (1.0 - lambda)*1;
    }
}

void print_weights(double* weights, struct svm_problem prob)
{
	int i;
	for(i=0;i<prob.l;i++)
	{
			printf("%d:\t%f  ",i,weights[i]);
	}
	printf("\n");
}

double compute_zaehler(double* x_weights, double* y_weights, const struct svm_node* x,  struct svm_problem prob_y, struct svm_problem prob_x)
{
    double sum = 0.0;
    int j, k;
    for(j=0;j<prob_x.l;j++){
        for(k=0;k<prob_y.l;k++)
        {
            sum += x_weights[j] * y_weights[k] * dot(prob_x.x[j], prob_y.x[k]);
        }
    }

    for(j=0;j<prob_x.l;j++)
    {
        sum -= x_weights[j] * dot(prob_x.x[j], x);
    }

    for(k=0;k<prob_y.l;k++) {
        sum -= y_weights[k] * dot(prob_y.x[k], x);
    }

    sum += dot(x,x);
    return sum;
}


double compute_nenner(double* y_weights, const struct svm_node* x,  struct svm_problem prob_y)
{
    double sum = 0.0;
    int j, k;
    for(j=0;j<prob_y.l;j++){
        for(k=0;k<prob_y.l;k++)
        {
            sum += y_weights[j] * y_weights[k] * dot(prob_y.x[j], prob_y.x[k]);
        }
    }

    for(j=0;j<prob_y.l;j++)
    {
        sum -= 2 * y_weights[j] * dot(prob_y.x[j], x);
    }

    sum += dot(x,x);
    return sum;
}

double compute_distance(double* y_weights, double* x_weights, struct svm_problem prob_y, struct svm_problem prob_x)
{
    double sum = 0.0;
    int j, k;
    for(j=0;j<prob_x.l;j++){
        for(k=0;k<prob_y.l;k++)
        {
            sum -= 2 * x_weights[j] * y_weights[k] * dot(prob_x.x[j], prob_y.x[k]);
        }
    }

    for(j=0;j<prob_x.l;j++){
        for(k=0;k<prob_x.l;k++)
        {
            sum += x_weights[j] * x_weights[k] * dot(prob_x.x[j], prob_x.x[k]);
        }
    }

    for(j=0;j<prob_y.l;j++){
        for(k=0;k<prob_y.l;k++)
        {
            sum += y_weights[j] * y_weights[k] * dot(prob_y.x[j], prob_y.x[k]);
        }
    }

    return sum;
}

int main (int argc, const char ** argv)
{
    const char* filename;

    printf("Dateien einlesen... \n");
    if (argc < 2)
        //filename = "data/heart_scale.plus";
        filename = "data/triangle.left.txt";
    else
        filename = argv[1];

    read_problem(filename);

    prob_p = prob;

    if (argc < 3)
        //filename = "data/heart_scale.minus";
        filename = "data/triangle.right.txt";
    else
        filename = argv[2];

    read_problem(filename);

    prob_q = prob;

    printf("Anzahl der eingelesenen Datensaetze: P = %d   Q = %d \n", prob_p.l, prob_q.l);

    int max_index_q = compute_max_index(prob_q); //Todo: besseren namen finden
    int max_index_p = compute_max_index(prob_p);

    printf("Groesse der Vektoren: P = %d, Q = %d \n", max_index_p, max_index_q);

    int elements = max_index_p;
    if (elements < max_index_q)	elements = max_index_q;
    elements++; //add one for the stopping index -1

    printf("Gewichtsvektoren initialisieren.. \n");
    double* x_weights = (double *) malloc(prob_p.l * sizeof(double));

    double* y_weights = (double *) malloc(prob_q.l * sizeof(double));

    int i;
    x_weights[0] = 1.0;
    for (i=1;i<prob_p.l;i++)
        x_weights[i]=0.0;

    y_weights[0] = 1.0;
    for (i=1;i<prob_q.l;i++)
        y_weights[i]=0.0;

    // Step 1
    printf("Maximales Kreuzprodukt suchen... \n");
    double max_p;
    double max_q;

    int max_p_index = find_max_dotproduct_weights( x_weights, y_weights, prob_p, prob_q, &max_p);
    int max_q_index = find_max_dotproduct_weights( y_weights, x_weights, prob_q, prob_p, &max_q);


    int j;

    for (j=0;j<3;j++)
    {
        double lambda;
        if (max_p >= max_q)
        {
            double zaehler = compute_zaehler(y_weights, x_weights, prob_p.x[max_p_index], prob_p, prob_q);
            double nenner = compute_nenner(x_weights, prob_p.x[max_p_index], prob_p);

            lambda = zaehler / nenner;

            if(zaehler == 0.0 && nenner == 0.0) lambda = 0.0;
            if(lambda < 0.0)	lambda = 0.0;
            if(lambda > 1.0)	lambda = 0.0;

            printf("lambda = %f  zaehler = %f  nenner = %f \n", lambda, zaehler, nenner);


            add_to_weights(x_weights, lambda, max_p_index, prob_p);

            max_p_index = find_max_dotproduct_weights( x_weights, y_weights, prob_p, prob_q, &max_p); // max_p updaten

        }
        else
        {

            double zaehler = compute_zaehler(x_weights, y_weights, prob_q.x[max_q_index], prob_q, prob_p);
            double nenner = compute_nenner(y_weights, prob_q.x[max_q_index], prob_q);

            lambda = zaehler / nenner;

            if(zaehler == 0.0 && nenner == 0.0) lambda = 0.0;
            if(lambda < 0.0)	lambda = 0.0;
            if(lambda > 1.0)	lambda = 0.0;

            add_to_weights(y_weights, lambda, max_q_index, prob_q);
        }

        max_q_index = find_max_dotproduct_weights( y_weights, x_weights, prob_q, prob_p, &max_q); // max_q updaten
        max_p_index = find_max_dotproduct_weights( x_weights, y_weights, prob_p, prob_q, &max_p); // max_p updaten

        //duality gap
        // absolute duality gap

        double adg = max_p + max_q;

        printf("max_p = %f  max_q = %f ", max_p, max_q);
        printf("adg = %f ", adg);

        // relative duality gap
        // adg / ||p-q||^2 - adg
        // adg / <p-q, p-q> - adg

		printf("<x-y,x-y> = %f " , compute_distance(x_weights,y_weights,prob_p,prob_q));


        double rdg_nenner = compute_distance(x_weights,y_weights,prob_p,prob_q) - adg;
        double rdg;

        if (rdg_nenner <= 0)
        {
            printf("set huge value... ");
            rdg = HUGE_VAL;
        }
        else
        {
            rdg = adg / rdg_nenner;
        }

        printf("rdg = %f \n", rdg);
		print_weights(x_weights, prob_p);
		print_weights(y_weights, prob_q);

    }

    return 0;
}
