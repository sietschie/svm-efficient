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

void init_nodes(struct svm_node* vector, int elements)
{
    int i;
    for (i=0;i< elements-1 ; i++)
        vector[i].index = i+1;
    vector[elements-1].index = -1;
}

void compute_vector_from_weights(double* weights, struct svm_node* vector, const struct svm_problem prob)
{
    int i;
    for (i=0;i<prob.l;i++)
    {
        struct svm_node* out= vector;
        struct svm_node* in = prob.x[0];
        while (in->index != -1 && out->index != -1)
        {
            if (in->index == out->index)
            {
                if (i==0)
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
                if (in->index > out->index) {
                    out->value = 0.0;
                    ++out;
                }
                else
                    ++in;
            }
        }
    }
}

double compute_dot_product_with_differences4(const struct svm_node* a, const struct svm_node* b, const struct svm_node* c, const struct svm_node* d) // <a - b, c - d>
{
//    printf(" start of function.. \n");
//    printf(" a=%d  b=%d  c=%d  d=%d \n", a, b, c, d);
/*    printf(" a=%d\n", a->index);
    printf(" b=%d\n", b->index);
    printf(" c=%d\n", c->index);
    printf(" d=%d\n", d->index);*/
//    printf(" a=%d  b=%d  c=%d  d=%d \n", a->index, b->index, c->index, d->index);
    double sum = 0;
    while (a->index != -1 && b->index != -1 && c->index != -1 && d->index != -1)
    {
        if ((a->index == b->index) && (c->index == d->index) &&  (a->index == d->index))
        {
            //printf("sum = %f  a = %f  b = %f  c = %f  d = %f\n", sum, a->value, b->value, c->value, d->value);
            sum += (a->value - b->value)*(c->value - d->value);
            ++a;
            ++b;
            ++c;
            ++d;
        }
        else
        {
            //printf("missing element 2.. a=%d b=%d c=%d e=%d\n", a->index, b->index, c->index, d->index); // todo: was tun bei fehlenden werten?

            if (a->index < b->index || a->index < c->index || a->index < d->index)
                ++a;

            if (b->index < a->index || b->index < c->index || b->index < d->index)
                ++b;

            if (c->index < a->index || c->index < b->index || c->index < d->index)
                ++c;

            if (d->index < a->index || d->index < b->index || d->index < c->index)
                ++d;
            //printf("missing element 3.. a=%d b=%d c=%d e=%d\n", a->index, b->index, c->index, d->index); // todo: was tun bei fehlenden werten?
        }
    }
    return sum;
}


void print_vector( const struct svm_node* vector)
{
    while (vector->index != -1)
    {
        printf("%d: %f \n", vector->index, vector->value);
        ++vector;
    }
}


int find_max_dotproduct(const struct svm_node* x, const struct svm_node* y, struct svm_problem prob, double *max_ret)
{
//    printf("FIND MAX DOT PRODUCT\n");
    int i;
    double max = -HUGE_VAL; //compute_dot_product_with_differences(prob.x[0], x, y);
    int max_index = -1;
    for (i=0;i<prob.l;i++)
//	for(i=0;i<10;i++)
//	i=2;
    {
        //print_vector(x);
        double dotproduct = compute_dot_product_with_differences4(prob.x[i], x, y, x);
        if (dotproduct > max)
        {
            max = dotproduct;
            max_index = i;
        }
//        printf("index = %d dotproduct = %f  max = %f \n", i, dotproduct, max);
        //print_vector(prob.x[0]);
        //print_vector(x);
        //print_vector(y);
        //print_vector(x);
    }
    (*max_ret) = max; //todo: watch closely
    return max_index;
}


void add_proportional(struct svm_node* a, const struct svm_node* b, double lambda)
{
    while (a->index != -1 && b->index != -1)
    {
        if (a->index == b->index)
        {
            //printf("a = %f  b = %f  lambda = %f \n", a->value, b->value, lambda);
            a->value = lambda*a->value + (1.0 - lambda)* b->value;
            //printf("a = %f  \n", a->value);
            ++a;
            ++b;
        }
        else
        {
        	//printf("skip: %d  %d \n", a->index, b->index);
            if (a->index > b->index)
                ++b;
            else
            {
                ++a;
            }
        }
    }
}

void add_to_weights(double* weights, double lambda, int index, struct svm_problem prob)
{
    int i;
    for (i=0;i<prob.l;i++)
    {
    	//printf("     old=%f   ",weights[i]);
        if (i!=index)
            weights[i] *= lambda;
        else
            weights[i] = weights[i]*lambda + (1.0 - lambda)*1;
        //printf("    new=%f  \n", weights[i]);
    }
}

void print_weights(double* weights, struct svm_problem prob)
{
	int i;
	for(i=0;i<prob.l;i++)
	{
		//if(weights[i] != 0.0)
			printf("%d:\t%f\n",i,weights[i]);
	}
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
    struct svm_node* x = (struct svm_node *) malloc(elements * sizeof(struct svm_node));
    init_nodes(x, elements);

    double* y_weights = (double *) malloc(prob_q.l * sizeof(double));
    struct svm_node* y = (struct svm_node *) malloc(elements * sizeof(struct svm_node));
    init_nodes(y, elements);

    int i;
    x_weights[0] = 1.0;
    for (i=1;i<prob_p.l;i++)
        x_weights[i]=0.0;

    y_weights[0] = 1.0;
    for (i=1;i<prob_q.l;i++)
        y_weights[i]=0.0;



    printf("Vektoren aus den Gewichtsvektoren berechnen.. \n");
    compute_vector_from_weights(x_weights, x, prob_p);
    compute_vector_from_weights(y_weights, y, prob_q);


    // Step 1
    printf("Maximales Kreuzprodukt suchen... \n");
    double max_p;
    int max_p_index = find_max_dotproduct( x, y, prob_p, &max_p);

    double max_q;
    int max_q_index = find_max_dotproduct( y, x, prob_q, &max_q);

    int j;

    for (j=0;j<3;j++)
    {
//        printf(" ... P = %d (%f), Q = %d (%f)\n", max_p_index, max_p, max_q_index, max_q);
        double lambda;
            printf(" max_q_index = %d , max_p_index = %d \n", max_q_index, max_p_index);

        if (max_p >= max_q)
        {
//            printf("  max_p > max_q \n");
            double zaehler = compute_dot_product_with_differences4(y, prob_p.x[max_p_index], x, prob_p.x[max_p_index]);
            double nenner = compute_dot_product_with_differences4(x, prob_p.x[max_p_index], x, prob_p.x[max_p_index]);

//            printf("zaehler=%f  nenner=%f \n", zaehler, nenner);

            lambda = zaehler / nenner;

            if(lambda < 0.0)	lambda = 0.0;
            if(lambda > 1.0)	lambda = 1.0;


			//print_weights(x_weights,prob_p);

            printf("lambda = %f  zaehler = %f  nenner = %f \n", lambda, zaehler, nenner);

            add_proportional(x,prob_p.x[max_p_index], lambda);
            add_to_weights(x_weights, lambda, max_p_index, prob_p);

			//print_weights(x_weights,prob_p);


            max_p_index = find_max_dotproduct( x, y, prob_p, &max_p); // max_p updaten
            max_q_index = find_max_dotproduct( y, x, prob_q, &max_q); // max_q updaten
        }
        else
        {
//            printf("  max_p <= max_q \n");
            // Gilbertschritt im Polytop Q
            // \lambda = <p-q_i, q-q_i> / <q-q_i, q-q_i>

            // <x - max_q_index, y - max_q_index> /  <y - max_q_index, y - max_q_index>


            double zaehler = compute_dot_product_with_differences4(x, prob_q.x[max_q_index], y, prob_q.x[max_q_index]);
            double nenner = compute_dot_product_with_differences4(y, prob_q.x[max_q_index], y, prob_q.x[max_q_index]);

//            printf("zaehler=%f  nenner=%f \n", zaehler, nenner);

            lambda = zaehler / nenner;

//            if(zaehler == 0.0 && nenner == 0.0) lambda = 0.0;
            if(lambda < 0.0)	lambda = 0.0;
            if(lambda > 1.0)	lambda = 1.0;

            printf("lambda = %e  zaehler = %e  nenner = %e \n", lambda, zaehler, nenner);
            add_proportional(y,prob_q.x[max_q_index], lambda);
            add_to_weights(y_weights, lambda, max_q_index, prob_q);

            max_q_index = find_max_dotproduct( y, x, prob_q, &max_q); // max_q updaten
            max_p_index = find_max_dotproduct( y, x, prob_p, &max_p); // max_p updaten
        }

        //duality gap
        // absolute duality gap
        // berechne max_p und max_q neu, dann adg = max_p+max_q

        double adg = max_p + max_q;

        printf("max_p = %f  max_q = %f ", max_p, max_q);
        printf("adg = %f ", adg);

        // relative duality gap
        // adg / ||p-q||^2 - adg
        // adg / <p-q, p-q> - adg

		printf("<x-y,x-y> = %f" , compute_dot_product_with_differences4(x,y,x,y));

		//print_vector(x);
		//print_vector(y);


        double rdg_nenner = compute_dot_product_with_differences4(x,y,x,y) - adg;
        double rdg;

        if (rdg_nenner <= 0)
        {
            printf("set huge value... \n");
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
