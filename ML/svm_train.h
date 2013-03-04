#ifndef _SVM_H
#define _SVM_H

#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libsvm-3.16/svm.h"
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

namespace SVM {

template <class T>
signed label(const T & r)
{
	return r > 0 ? 1 : -1;
}

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
int cross_validation;
int nr_fold;

static char *line = NULL;
static int max_line_len;


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


void init_param()
{
	void (*print_func)(const char*) = NULL;	// default printing to stdout
	param.svm_type = C_SVC;
	param.kernel_type = LINEAR;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 2000;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	cross_validation = 0;

	svm_set_print_string_function(print_func);


}

template <class T>
void init_prob(const Matrix<T> & x, const Matrix<T> & y)
{
	prob.l = x.row;
	prob.y = Malloc(double, prob.l);
	prob.x = Malloc(struct svm_node *, prob.l);
	x_space = Malloc(struct svm_node, x.row * x.col + y.row);

	unsigned j = 0;
	int inst_max_index = 0, max_index = 0;
	for (signed i = 0; i < prob.l; ++i)
	{
		inst_max_index = -1;
		prob.x[i] = &x_space[j];
		prob.y[i] = y[i][0];
		
		for (unsigned idx = 1; idx < x.col; ++idx)
		{
			x_space[j].index = idx;
			inst_max_index = x_space[j].index;
			x_space[j].value = x[i][idx];
			++j;
		}

		if(inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[j++].index = -1;
	}

	if (param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if (param.kernel_type == PRECOMPUTED)
		for (signed i = 0; i < prob.l; ++i)
		{
			if (prob.x[i][0].index != 0)
			{
				fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
				exit(1);
			}
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
				exit(1);
			}
		}

}

// read in a problem (in svmlight format)

void read_prob(const char *filename)
{
	int elements, max_index, inst_max_index, i, j;
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;

	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	prob.l = 0;
	elements = 0;

	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	while(readline(fp)!=NULL)
	{
		char *p = strtok(line," \t"); // label

		// features
		while(1)
		{
			p = strtok(NULL," \t");
			if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			++elements;
		}
		++elements;
		++prob.l;
	}
	rewind(fp);

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node,elements);

	max_index = 0;
	j=0;
	for(i=0;i<prob.l;i++)
	{
		inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
		readline(fp);
		prob.x[i] = &x_space[j];
		label = strtok(line," \t\n");
		if(label == NULL) // empty line
			exit_input_error(i+1);
	
		prob.y[i] = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
			exit_input_error(i+1);
		
		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");
			if(val == NULL)
				break;

			errno = 0;
			x_space[j].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
				exit_input_error(i+1);
			else
				inst_max_index = x_space[j].index;

			errno = 0;
			x_space[j].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(i+1);

			++j;

		}

		if(inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if (prob.x[i][0].index != 0)
			{
				fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
				exit(1);
			}
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
				exit(1);
			}
		}

	fclose(fp);
}

int learn_sv_from_file(const char * filename, const char * model_file_name)
{
	init_param();
	read_prob(filename);
	model = svm_train(&prob, &param);

	double err_in = 0.0;
	if (model)
	{
		for (signed i = 0; i < prob.l; ++i)
			err_in += (prob.y[i] != svm_predict(model, prob.x[i]));

		if(svm_save_model(model_file_name,model))
		{
			fprintf(stderr, "can't save model to file %s\n", model_file_name);
			exit(1);
		}
	}

	err_in /= prob.l;

	cerr << "SVM ein=" << err_in << endl;

	free(prob.y);
	free(prob.x);
	free(x_space);

	return * model->nSV;
}

template <class T>
int learn_sv(const Matrix<T> & x, const Matrix<T> & y, const char * model_file_name)
{
	init_param();
	init_prob(x, y);
	model = svm_train(&prob, &param);

	double err_in = 0.0;
	if (model)
	{
		for (signed i = 0; i < prob.l; ++i)
			err_in += (prob.y[i] != svm_predict(model, prob.x[i]));

		if(svm_save_model(model_file_name,model))
		{
			fprintf(stderr, "can't save model to file %s\n", model_file_name);
			exit(1);
		}
	}

	err_in /= prob.l;

	cerr << "SVM ein=" << err_in << endl;

	free(prob.y);
	free(prob.x);
	free(x_space);

	return *model->nSV;
}

template <class T>
T predict(const Matrix<T> & x, const Matrix<T> & y, struct svm_model * model)
{
	struct svm_parameter param;
	struct svm_problem prob;		// set by read_problem
	struct svm_node * x_space = init(x, y, prob, param);

	if (svm_check_probability_model(model))
		cerr << "positive check for probability" << endl;
	double eout = 0.0;
	if (model)
	{ 
		for (signed i = 0; i < prob.l; ++i)
		{
			T g = svm_predict(model, prob.x[i]);
			//cerr << "g=" << g << " y=" << prob.y[i] << endl; 
			eout += (prob.y[i] != g);
		}
	}

	free(prob.y);
	free(prob.x);
	free(x_space);

	eout /= prob.l;

	cerr << "SVM eout=" << eout << endl;
	return eout;
}

} // namespace
#endif