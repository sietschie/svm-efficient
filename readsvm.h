#ifndef READSVM_H_INCLUDED
#define READSVM_H_INCLUDED

#include "globals.h"


void read_problem(const char *filename);
int svm_save_model(const char *model_file_name, const struct svm_model* model);


#endif // READSVM_H_INCLUDED
