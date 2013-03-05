#ifndef _SVM_H
#define _SVM_H

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

class SVM_Train {

public:
	SVM_Train();
	~SVM_Train();

	double do_cross_validation(const char *filename, int nr_fold);
	void do_train(const char * input_file_name, const char * model_file_name);

private:

	char* readline(FILE *input);
	void read_problem(const char *filename);

public:

	svm_parameter param;	

private:
	svm_problem prob;		// set by read_problem
	svm_model *model;
	svm_node *x_space;

	int nr_fold;

	char *line;
	int max_line_len;
};

} // namespace
#endif