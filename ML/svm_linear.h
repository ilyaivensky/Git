#ifndef _SVM_LINEAR_H
#define _SVM_LINEAR_H

#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "liblinear-1.93/linear.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define INF HUGE_VAL

namespace SVM_LIN {

template <class T>
signed label(const T & r)
{
	return r > 0 ? 1 : -1;
}

template <class T>
feature_node * init(const Matrix<T> & x, const Matrix<T> & y, 
						problem & prob, parameter & param)
{
	// default values
	param.solver_type = L2R_L2LOSS_SVC_DUAL;
	param.C = 2000;
	param.eps = INF;// see setting below
	param.p = 0.1;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;

	//prob.bias = -1;
	prob.l = x.row;
	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct feature_node *,prob.l);
	struct feature_node * x_space = Malloc(struct feature_node, x.row * x.col + y.row /*elements+prob.l*/);

	if(param.eps == INF)
	{
		switch(param.solver_type)
		{
			case L2R_LR:
			case L2R_L2LOSS_SVC:
				param.eps = 0.01;
				break;
			case L2R_L2LOSS_SVR:
				param.eps = 0.001;
				break;
			case L2R_L2LOSS_SVC_DUAL:
			case L2R_L1LOSS_SVC_DUAL:
			case MCSVM_CS:
			case L2R_LR_DUAL:
				param.eps = 0.1;
				break;
			case L1R_L2LOSS_SVC:
			case L1R_LR:
				param.eps = 0.01;
				break;
			case L2R_L1LOSS_SVR_DUAL:
			case L2R_L2LOSS_SVR_DUAL:
				param.eps = 0.1;
				break;
		}
	}

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

		if (inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[j++].index = -1;
	}

	if(prob.bias >= 0)
	{
		prob.n=max_index+1;
		for(signed i=1;i<prob.l;i++)
			(prob.x[i]-2)->index = prob.n;
		x_space[j-2].index = prob.n;
	}
	else
		prob.n=max_index;

	return x_space;
}

template <class T>
model * learn(const Matrix<T> & x, const Matrix<T> & y)
{
	struct parameter param;
	struct problem prob;
	struct feature_node * x_space = init(x, y, prob, param);
	struct model* model_ = train(&prob, &param);

	double err_in = 0.0;
	if (model_)
	{
		for (signed i = 0; i < prob.l; ++i)
			err_in += (prob.y[i] != predict(model_, prob.x[i]));
	}

	err_in /= prob.l;

	cerr << "SVM ein=" << err_in << endl;

	free(prob.y);
	free(prob.x);
	free(x_space);

	return model_;
}

template <class T>
T predict(const Matrix<T> & x, const Matrix<T> & y, struct model * model)
{
	struct parameter param;
	struct problem prob;		// set by read_problem
	struct feature_node * x_space = init(x, y, prob, param);

	double eout = 0.0;
	if (model)
	{ 
		for (signed i = 0; i < prob.l; ++i)
		{
			T g = predict(model, prob.x[i]);
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