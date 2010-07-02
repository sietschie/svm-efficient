#ifndef KERNEL_H_INCLUDED
#define KERNEL_H_INCLUDED

double kernel_linear(int set1, int element1, int set2, int element2);
double kernel_poly(int set1, int element1, int set2, int element2);
double kernel_rbf(int set1, int element1, int set2, int element2);
double kernel_sigmoid(int set1, int element1, int set2, int element2);
double (*kernel)(int, int, int, int);

#endif // KERNEL_H_INCLUDED
