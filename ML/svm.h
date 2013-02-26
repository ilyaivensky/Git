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

template <class T>
svm_node * init(const Matrix<T> & x, const Matrix<T> & y, 
						svm_problem & prob, svm_parameter & param)
{
	param.svm_type = C_SVC;
	param.kernel_type = LINEAR;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 100;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	prob.l = x.row;
	prob.y = Malloc(double, prob.l);
	prob.x = Malloc(struct svm_node *, prob.l);
	struct svm_node * x_space = Malloc(struct svm_node, x.row * x.col + y.row);

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

	return x_space;
}

template <class T>
svm_model * learn(const Matrix<T> & x, const Matrix<T> & y)
{
	struct svm_parameter param;
	struct svm_problem prob;		// set by read_problem
	struct svm_node * x_space = init(x, y, prob, param);
	struct svm_model * model = svm_train(&prob, &param);

	double err_in = 0.0;
	if (model)
	{
		for (signed i = 0; i < prob.l; ++i)
			err_in += (prob.y[i] != svm_predict(model, prob.x[i]));
	}

	err_in /= prob.l;

	cerr << "SVM ein=" << err_in << endl;

	free(prob.y);
	free(prob.x);
	free(x_space);

	return model;
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