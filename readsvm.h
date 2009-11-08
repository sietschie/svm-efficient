#ifndef READSVM_H_INCLUDED
#define READSVM_H_INCLUDED

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

void read_problem(const char *filename);

#endif // READSVM_H_INCLUDED
