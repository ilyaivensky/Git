#include "svm_train.h"
#include <cstdio>
#include <cctype>
#include <exception>
#include <iostream>
#include <sstream>

#include <conio.h>

using namespace std;

namespace SVM {

namespace {

void exit_input_error(int line_num)
{
	stringstream message;
	message << "Wrong input format at line " << line_num;
	throw exception(message.str().c_str());
}

}

SVM_Train::SVM_Train() : 
	model(NULL),
	x_space(NULL),
	nr_fold(0),
	line(NULL)
{
	// default values
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
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
}

SVM_Train::~SVM_Train()
{
	//svm_free_and_destroy_model(&model);
}

char* SVM_Train::readline(FILE *input)
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

void SVM_Train::read_problem(const char *filename)
{
	int elements, max_index, inst_max_index, i, j;
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;

	if(fp == NULL)
	{
		stringstream message;
		message << "can't open input file " << filename;
		throw exception(message.str().c_str());
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
				stringstream message;
				message << "Wrong input format: first column must be 0:sample_serial_number";
				throw exception(message.str().c_str());		
			}
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				stringstream message;
				message << "Wrong input format: sample_serial_number out of range";
				throw exception(message.str().c_str());		
			}
		}

	fclose(fp);
}

void SVM_Train::do_train(const char * input_file_name, const char * model_file_name)
{
	const char *error_msg;
	//parse_command_line(argc, argv, input_file_name, model_file_name);
	read_problem(input_file_name);
	error_msg = svm_check_parameter(&prob,&param);

	if(error_msg)
		throw exception(error_msg);

	model = svm_train(&prob,&param);
	if(svm_save_model(model_file_name,model))
	{
		stringstream message;
		message << "can't save model to file " << model_file_name;
		throw exception(message.str().c_str());		
	}
	svm_free_and_destroy_model(&model);
	svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);
	free(x_space);
	free(line);
}

double SVM_Train::do_cross_validation(const char * input_file_name, int n_fold)
{
	nr_fold = n_fold;
	read_problem(input_file_name);

	const char * error_msg = svm_check_parameter(&prob,&param);
	if(error_msg)
		throw exception(error_msg);

	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double,prob.l);

	double retval = 0.0;

	svm_cross_validation(&prob,&param,nr_fold,target);
	if(param.svm_type == EPSILON_SVR ||
	   param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}

		retval = total_error/prob.l;
		
		printf("Cross Validation Mean squared error = %g\n",retval);
		printf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;

		retval = 1.0 - ((double)total_correct/prob.l);
		printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
	}

	free(target);
	svm_destroy_param(&param);
	param.weight_label = NULL;
	param.weight = NULL;
	free(prob.y);
	free(prob.x);
	free(x_space);
	free(line);

	return retval;
}

} // namespace