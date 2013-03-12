 /*                                                                 -*- C++ -*-
 * File: svm_train.h
 * 
 * Author: Ellie Ivensky
 * Copyright (c) 2013 Idilia Inc, All rights reserved.
 * 
 * Created on: Mar 4, 2013
 *
 * algorithms/machine_learning/svm/svm_train.h
 *
 * Description:
 *   C++ interface for libSVM 
 *   
 */
 
#ifndef _SVM_H
#define _SVM_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libsvm-3.16/svm.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

namespace SV {

class SVM_Train {

public:
	
	SVM_Train();
	~SVM_Train();

	svm_parameter & parameter() { return param; }
	const svm_parameter & parameter() const { return param; }
	
	double do_cross_validation(const char *filename, int nr_fold);
	void do_train(const char * input_file_name, const char * model_file_name);

private:

	char* readline(FILE *input);
	void read_problem(const char *filename);
	
private:
	svm_parameter param;	
	svm_problem prob;		// set by read_problem
	svm_model *model;
	svm_node *x_space;

	int nr_fold;

	char *line;
	int max_line_len;
};

} // namespace

#endif
