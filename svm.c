#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "svm.h"
#include "globals.h"
#include "kernel.h"
#include "cache.h"

double* dot_same[2];

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

int find_max(int p, double *dot_yi_x, double* dot_xi_x, double dot_xi_yi, double dot_xi_xi, double *max_p) {
    // find max
    int max_p_index = -1;
    *max_p = -HUGE_VAL;
    int i;
    for (i=0;i<prob[p].l;i++) {
        double sum = dot_yi_x[i] - dot_xi_x[i] - dot_xi_yi + dot_xi_xi;
        //printf("sum:%f = dot_yi_x[%d]:%f - dot_xi_x[%d]:%f - dot_xi_yi:%f + dot_xi_xi:%f \n", sum, i, dot_yi_x[i], i,  dot_xi_x[i], dot_xi_yi, dot_xi_xi);
        if(sum > *max_p)
        {
            *max_p = sum;
            max_p_index = i;
        }
    }
    return max_p_index;
}

double compute_zaehler(double dot_xi_yi, double* dot_yi_x, double* dot_xi_x, int p, int max_p_index ) {
    double zaehler = dot_xi_yi - dot_yi_x[max_p_index] - dot_xi_x[max_p_index] + dot_same[p][max_p_index];
    return zaehler;
}

double compute_nenner(double dot_xi_xi, double* dot_xi_x, int p, int max_p_index) {
    double nenner = dot_xi_xi - 2* dot_xi_x[max_p_index] + dot_same[p][max_p_index];
    return nenner;
}

double update_xi_xi(double dot_xi_xi, double* dot_xi_x, int p, int max_p_index, double lambda) {
    dot_xi_xi = lambda * lambda * dot_xi_xi
            + 2 * lambda * (1.0 - lambda) * dot_xi_x[max_p_index]
            + (1.0 - lambda)*(1.0 - lambda) * dot_same[p][max_p_index];
    return dot_xi_xi;
}

double update_xi_yi(double dot_xi_yi, double* dot_yi_x, int max_p_index, double lambda) {
    dot_xi_yi = lambda * dot_xi_yi + (1.0 - lambda) * dot_yi_x[max_p_index];
    return dot_xi_yi;
}

void update_xi_x(double* dot_xi_x, int p, int p2, int max_p_index, double lambda) {
    //printf("update_xi_x(): %d %d %d \n", p, p2, max_p_index);
    double* computed_kernels = get_element(max_p_index, p);

    int i;
    for (i=0;i<prob[p2].l;i++) {
        //dot_xi_x[i]= dot_xi_x[i] * lambda + (1.0 - lambda) * kernel(p, max_p_index, p2, i);
        int offset = p2 * prob[0].l;
        dot_xi_x[i]= dot_xi_x[i] * lambda + (1.0 - lambda) * computed_kernels[ offset + i  ]; //(p, max_p_index, p2, i);
        //printf(" %d - %f, max_p_index = %d, offset = %d\n", i, computed_kernels[ offset + i  ], max_p_index, offset);
    }
    //printf("\n");
}

void compute_weights(double *x_weights, double* y_weights)
{
        // init cache
    init(param.cache_size, prob[0].l + prob[1].l);

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
	
    dot_same[0] = (double *) malloc(prob[0].l * sizeof(double));
    dot_same[1] = (double *) malloc(prob[1].l * sizeof(double));

    // initialisieren
    for (i=0;i<prob[0].l;i++) {
        dot_xi_x[i]=kernel(0, 0, 0, i);
        dot_yi_x[i]=kernel(1, 0, 0, i);
        dot_same[0][i]=kernel(0,i,0,i);
    }

    for (i=0;i<prob[1].l;i++) {
        dot_xi_y[i]=kernel(0, 0, 1, i);
        dot_yi_y[i]=kernel(1, 0, 1, i);
        dot_same[1][i]=kernel(1,i,1,i);
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

    for (j=0;j<param.maximum_iterations;j++)
    {
        //printf("j = %d \n", j);
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

            //printf("max_p: \n");
            update_xi_x(dot_xi_x, 0, 0, max_p_index, lambda);

            update_xi_x(dot_xi_y, 0, 1, max_p_index, lambda);
        //printf("max_p = %f  max_q = %f zaehler = %f nenner = %f lambda = %f\n", max_p, max_q, zaehler, nenner, lambda);
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

            //printf("max_q: \n");
            update_xi_x(dot_yi_y, 1, 1, max_q_index, lambda);

            update_xi_x(dot_yi_x, 1, 0, max_q_index, lambda);
        //printf("max_p = %f  max_q = %f zaehler = %f nenner = %f lambda = %f\n", max_p, max_q, zaehler, nenner, lambda);
        }
        // find max
        max_p_index = find_max(0, dot_yi_x, dot_xi_x, dot_xi_yi, dot_xi_xi, &max_p);
        max_q_index = find_max(1, dot_xi_y, dot_yi_y, dot_xi_yi, dot_yi_yi, &max_q);

        //duality gap
        // absolute duality gap

        double adg = max_p + max_q;

        //printf("max_p = %f  max_q = %f ", max_p, max_q);
        //printf("adg = %f ", adg);

        // relative duality gap
        // adg / ||p-q||^2 - adg
        // adg / <p-q, p-q> - adg

        double distance = dot_xi_xi + dot_yi_yi - 2 * dot_xi_yi;


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

		if(param.verbosity == 2) 
		{
			printf("max[0] = %d (%f)   max[1] = %d (%f)  \n", max_p_index, max_p, max_q_index, max_q);
		}

		if( param.verbosity >= 1 )
		{
			printf("iter = %d ", j);
			printf("dist = %e " , distance);
			printf("adg = %e " , adg);
			printf("rdg = %e \n", rdg);
		}

	    //rho = - dot_xi_yi + dot_xi_xi - (dot_xi_xi + dot_yi_yi - 2 * dot_xi_yi)/2;
        rho = dot_xi_yi - dot_xi_xi - (dot_xi_xi + dot_yi_yi - 2 * dot_xi_yi)/2;
        //printf("xi_xi = %f   yi_yi = %f   xi_yi = %f \n", dot_xi_xi, dot_yi_yi, dot_xi_yi);
		if( rdg < param.eps )
			break;
    }
}
