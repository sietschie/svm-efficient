#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "readsvm.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))



static char *line = NULL;
static int max_line_len;
struct svm_node *x_space[2];
static double C = 5.5;

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
		if(line[0] == '+')
			current_set = 0;
		else
			current_set = 1;

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
		if(line[0] == '+')
			current_set = 0;
		else
			current_set = 1;


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
