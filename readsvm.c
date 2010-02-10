#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "readsvm.h"
#include "globals.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))



static char *line = NULL;
static int max_line_len;
struct svm_node *x_space[2];

int max_index = 0;

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
	int elements[2], inst_max_index, i[2], j[2];
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;

	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	prob[0].l = 0;
	prob[1].l = 0;
	elements[0] = 0;
	elements[1] = 0;

	int current_set;

	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	while(readline(fp)!=NULL)
	{
		if(line[0] == '-')
			current_set = 1;
		else
			current_set = 0;

		char *p = strtok(line," \t"); // label

		// features
		while(1)
		{
			p = strtok(NULL," \t");
			if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			++elements[current_set];
		}
		++elements[current_set];
		++prob[current_set].l;
	}
	rewind(fp);

	prob[0].y = Malloc(double,prob[0].l);
	prob[0].x = Malloc(struct svm_node *,prob[0].l);
	x_space[0] = 	Malloc(struct svm_node,elements[0]);

	prob[1].y = Malloc(double,prob[1].l);
	prob[1].x = Malloc(struct svm_node *,prob[1].l);
	x_space[1] = 	Malloc(struct svm_node,elements[1]);

//	max_index = 0;
	j[0]=0;
	j[1]=0;

	i[0] = 0;
	i[1] = 0;

	while( i[0] < prob[0].l || i[1] < prob[1].l )
//	for(i=0;i<prob.l;i++)
	{
		inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
		readline(fp);
		if(line[0] == '-')
			current_set = 1;
		else
			current_set = 0;


		prob[current_set].x[i[current_set]] = &x_space[current_set][j[current_set]];
		label = strtok(line," \t");
		prob[current_set].y[i[current_set]] = strtod(label,&endptr);
		if(endptr == label)
			exit_input_error(i[0]+i[1]+1);

		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");

			if(val == NULL)
				break;

			errno = 0;
			x_space[current_set][j[current_set]].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[current_set][j[current_set]].index <= inst_max_index)
				exit_input_error(i[0]+i[1]+1);
			else
				inst_max_index = x_space[current_set][j[current_set]].index;

			errno = 0;
			x_space[current_set][j[current_set]].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(i[0]+i[1]+1);

			++j[current_set];
		}

		if(inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[current_set][j[current_set]++].index = -1;
		i[current_set]++;
	}

	fclose(fp);
}

const char *svm_type_table[] =
{
	"c_svc","nu_svc","one_class","epsilon_svr","nu_svr",NULL
};

const char *kernel_type_table[]=
{
	"linear","polynomial","rbf","sigmoid","precomputed",NULL
};


int svm_save_model(const char *model_file_name, const struct svm_model* model)
{
	FILE *fp = fopen(model_file_name,"w");
	if(fp==NULL) return -1;

	const struct svm_parameter param = model->param;

	fprintf(fp,"svm_type %s\n", svm_type_table[param.svm_type]);
	fprintf(fp,"kernel_type %s\n", kernel_type_table[param.kernel_type]);

	if(param.kernel_type == POLY)
		fprintf(fp,"degree %d\n", param.degree);

	if(param.kernel_type == POLY || param.kernel_type == RBF || param.kernel_type == SIGMOID)
		fprintf(fp,"gamma %g\n", param.gamma);

	if(param.kernel_type == POLY || param.kernel_type == SIGMOID)
		fprintf(fp,"coef0 %g\n", param.coef0);

//	int nr_class = model->nr_class;
	int nr_class = 2;


	int i;
	int counter[2];
	counter[0] = 0;
	for(i=0;i<prob[0].l;i++)
	{
        if(model->weights[0][i] != 0.0)
            counter[0]++;
    }

	counter[1] = 0;
	for(i=0;i<prob[1].l;i++)
	{
        if(model->weights[1][i] != 0.0)
            counter[1]++;
    }


	int l = counter[0] + counter[1];
	fprintf(fp, "nr_class %d\n", nr_class);
	fprintf(fp, "total_sv %d\n",l);

//	{
//		fprintf(fp, "rho");
//        fprintf(fp," %g",model->rho);
//		fprintf(fp, "\n");
//	}
	{
		fprintf(fp, "rho");
        fprintf(fp," %g",model->rho);
		fprintf(fp, "\n");
	}

	if(model->label)
	{
		fprintf(fp, "label");
		int i;
		for(i=0;i<nr_class;i++)
			fprintf(fp," %d",model->label[i]);
		fprintf(fp, "\n");
	}

	if(model->nSV)
	{
		fprintf(fp, "nr_sv");
		int i;
		for(i=0;i<nr_class;i++)
			fprintf(fp," %d",model->nSV[i]);
		fprintf(fp, "\n");
	}

	fprintf(fp, "SV\n");
//	const double * const *sv_coef = model->sv_coef;
//	const struct svm_node * const *SV = model->SV;

    int j;
    for(j=0;j<2;j++)
	for(i=0;i<prob[j].l;i++)
	{
        //printf("neuer vektor: \n");        
	    if(model->weights[j][i] != 0.0)
	    {
            //printf(" mit gutem coeffizienten: i = %d, j = %d\n", i, j);
	        double sign = 1.0;
	        if(j==1)
                sign = -1.0;
            fprintf(fp, "%.16g ",sign * model->weights[j][i]);

            const struct svm_node *p = model->SV[j][i];

			while(p->index != -1)
			{
				//printf("p->index = %d p->value = %f\n", p->index, p->value);
				//printf("%d:%.8g ",p->index,p->value);
				fprintf(fp,"%d:%.8g ",p->index,p->value);
				p++;
				//printf("p->index = %d p->value = %f\n", p->index, p->value);
			}
		//printf("\n");
            fprintf(fp, "\n");
	    }
	}

	if (ferror(fp) != 0 || fclose(fp) != 0) return -1;
	else return 0;
}


struct svm_model *svm_load_model(const char *model_file_name)
{
	FILE *fp = fopen(model_file_name,"rb");
	if(fp==NULL) return NULL;

	// read parameters

	struct svm_model *model = Malloc(struct svm_model,1);
	struct svm_parameter param = model->param;
//	model->rho = NULL;
//	model->probA = NULL;
//	model->probB = NULL;
//	model->label = NULL;
//	model->nSV = NULL;

	char cmd[81];
	while(1)
	{
		fscanf(fp,"%80s",cmd);

		if(strcmp(cmd,"svm_type")==0)
		{
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;svm_type_table[i];i++)
			{
				if(strcmp(svm_type_table[i],cmd)==0)
				{
					param.svm_type=i;
					break;
				}
			}
			if(svm_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown svm type.\n");
//				free(model->rho);
				free(model->label);
				free(model->nSV);
				free(model);
				return NULL;
			}
		}
		else if(strcmp(cmd,"kernel_type")==0)
		{
			fscanf(fp,"%80s",cmd);
			int i;
			for(i=0;kernel_type_table[i];i++)
			{
				if(strcmp(kernel_type_table[i],cmd)==0)
				{
					param.kernel_type=i;
					break;
				}
			}
			if(kernel_type_table[i] == NULL)
			{
				fprintf(stderr,"unknown kernel function.\n");
//				free(model->rho);
				free(model->label);
				free(model->nSV);
				free(model);
				return NULL;
			}
		}
		else if(strcmp(cmd,"degree")==0)
			fscanf(fp,"%d",&param.degree);
		else if(strcmp(cmd,"gamma")==0)
			fscanf(fp,"%lf",&param.gamma);
		else if(strcmp(cmd,"coef0")==0)
			fscanf(fp,"%lf",&param.coef0);
//		else if(strcmp(cmd,"nr_class")==0)
//			fscanf(fp,"%d",&model->nr_class);
//		else if(strcmp(cmd,"total_sv")==0)
//			fscanf(fp,"%d",&model->l);
		else if(strcmp(cmd,"rho")==0)
		{
//			int n = model->nr_class * (model->nr_class-1)/2;
//			model->rho = Malloc(double,n);
//			for(int i=0;i<n;i++)
				fscanf(fp,"%lf",&model->rho);
		}
		else if(strcmp(cmd,"label")==0)
		{
//			int n = model->nr_class;
//			model->label = Malloc(int,n);
            int i;
			for(i=0;i<2;i++)
				fscanf(fp,"%d",&model->label[i]);
		}
		else if(strcmp(cmd,"nr_sv")==0)
		{
//			int n = model->nr_class;
//			model->nSV = Malloc(int,n);
			int i;
			for(i=0;i<2;i++)
				fscanf(fp,"%d",&model->nSV[i]);
		}
		else if(strcmp(cmd,"SV")==0)
		{
			while(1)
			{
				int c = getc(fp);
				if(c==EOF || c=='\n') break;
			}
			break;
		}
		else
		{
			fprintf(stderr,"unknown text in model file: [%s]\n",cmd);
//			free(model->rho);
//			free(model->label);
//			free(model->nSV);
			free(model);
			return NULL;
		}
	}

	// read sv_coef and SV

	int elements = 0;
	long pos = ftell(fp);

	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	char *p,*endptr,*idx,*val;

	while(readline(fp)!=NULL)
	{
		p = strtok(line,":");
		while(1)
		{
			p = strtok(NULL,":");
			if(p == NULL)
				break;
			++elements;
		}
	}
	elements += model->l;

	fseek(fp,pos,SEEK_SET);

//	int m = model->nr_class - 1;
	int l = model->l;
//	model->sv_coef = Malloc(double *,m);


	int i;
	for(i=0;i<2;i++) {
		model->weights[i] = Malloc(double,model->nSV[i]);
        model->SV[i] = Malloc(struct svm_node*,model->nSV[i]);
	}

	struct svm_node *x_space = NULL;
	if(l>0) x_space = Malloc(struct svm_node,elements);

	int j=0;
	int cs = 0; // current set
	int ia[2];
	ia[0] = 0;
	ia[1] = 0;
	for(;ia[0] + ia[1] < l;)
	{
		readline(fp);

		p = strtok(line, " \t");

		double weight = strtod(p,&endptr);

		if(weight < 0.0)
		{
            weight *= -1.0;
            cs = 1;
		} else
            cs = 0;

		*(model->SV[ia[cs]]) = &x_space[j]; //TODO: warum * davor schreiben, damit es keine warnung gibt?

//		struct svm_node **SV[2]
//		struct svm_node *x_space[2];


		model->weights[cs][ia[cs]] = weight;
//		int k;
//		for(k=1;k<m;k++)
//		{
//			p = strtok(NULL, " \t");
//			model->sv_coef[k][i] = strtod(p,&endptr);
//		}

		while(1)
		{
			idx = strtok(NULL, ":");
			val = strtok(NULL, " \t");

			if(val == NULL)
				break;
			x_space[j].index = (int) strtol(idx,&endptr,10);
			x_space[j].value = strtod(val,&endptr);

			++j;
		}
		x_space[j++].index = -1;
		ia[cs]++;
	}
	free(line);

	if (ferror(fp) != 0 || fclose(fp) != 0)
		return NULL;

	model->free_sv = 1;	// XXX
	return model;
}
