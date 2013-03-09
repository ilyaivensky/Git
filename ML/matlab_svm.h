#ifndef _SVM_H
#define _SVM_H

#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"

//#pragma comment (lib,"C:/Program Files (x86)/MATLAB/R2012a Student/extern/lib/win32/microsoft/libeng.lib")
//#pragma comment (lib,"C:/Program Files (x86)/MATLAB/R2012a Student/extern/lib/win32/microsoft/libmx.lib")
//#pragma comment (lib,"C:/Program Files (x86)/MATLAB/R2012a Student/extern/lib/win32/microsoft/libut.lib")

namespace MATLAB_SV {

template <class T>
signed label(const T & r)
{
	return r > 0 ? 1 : -1;
}

template <class T> 
T linear_kernel(const vector<T> & v1, const vector<T> & v2)
{
	return dot_product(v1, v2);
}

template <class T> 
T second_order_polynomial_kernel(const vector<T> & v1, const vector<T> & v2)
{
	return pow(1 + dot_product(v1, v2), 2);
}

// x_non_scaled should include x0 = 1
template <class T>
unsigned learn(Engine * ep, const Matrix<T> & x_non_scaled, 
			   const Matrix<T> & y, vector<T> & w, 
			   T (*kernel)(const vector<T> &, const vector<T> &))
{	
	for (unsigned r = 0; r < x_non_scaled.row; ++r)
		if (x_non_scaled[r][0] != 1.0)
			throw exception("Wrong format of input");
	
	const double inf = std::numeric_limits<T>::infinity();
	const double zero = 0.0;
	const double minus_one = -1.0;

	/****************************************************** 
	* Remove x0 and scale
	*******************************************************/

	Matrix<T> x(x_non_scaled.row, x_non_scaled.col - 1);
	for (unsigned r = 0; r < x_non_scaled.row; ++r)
		for (unsigned c = 1; c < x_non_scaled.col; ++c)
			x[r][c - 1] = x_non_scaled[r][c];
	
	/***************************************************************************
	* Create H (quadr matrix)
	* The number of rows and columns in H must equal the number of elements of f
	*****************************************************************************/

	Matrix<T> quadr(x.row, x.row);
	for (unsigned i = 0; i < x.row; ++i)
		for (unsigned j = 0; j < x.row; ++ j)
			quadr[i][j] = y[i][0] * y[j][0] * kernel(x[i], x[j]);

	T * buff_H = (T*)malloc(sizeof(T) * quadr.row * quadr.col);
	
	for (unsigned r = 0; r < quadr.row; ++r)
		for (unsigned c = 0; c < quadr.col; ++c)
			buff_H[(r * quadr.col) + c] = quadr[r][c];

	mxArray * h = mxCreateDoubleMatrix(quadr.row, quadr.col, mxREAL);
	memcpy((void *)mxGetPr(h), (void *)buff_H,  quadr.row * quadr.col * sizeof(T));
	delete buff_H;

	engPutVariable(ep, "H", h);

	/***********************************************************************************
	* Create f 
	* Vector of doubles. Represents the linear term in the expression 1/2*x'*H*x + f'*x.
	* Setting it to '-1'
	* The number of rows and columns in H must equal the number of elements of f
	************************************************************************************/
	T * buff_f = (T*)malloc(sizeof(T) * quadr.row);
	for (unsigned r = 0; r < quadr.row; ++r)
		buff_f[r] = minus_one;
	
	mxArray * f = mxCreateDoubleMatrix(quadr.row, 1, mxREAL);
	memcpy((void *)mxGetPr(f), (void *)buff_f, quadr.row * sizeof(T));
	delete buff_f;

	engPutVariable(ep, "f", f);

	/****************************************************************************************
	* Create A and b (empty)
	* A - matrix of doubles. Represents the linear coefficients in the constraints A*x ≤ b.
	* b - vector of doubles. Represents the constant vector in the constraints A*x ≤ b.
	****************************************************************************************/
	mxArray * a =  mxCreateDoubleMatrix(0, 0, mxREAL);
	mxArray * b =  mxCreateDoubleMatrix(0, 0, mxREAL);

	engPutVariable(ep, "A", a);
	engPutVariable(ep, "b", b);

	/***************************************************************************************
	* Create Aeq 
	* Matrix of doubles. Represents the linear coefficients in the constraints Aeq*x = beq
	* The number of rows and columns in Aeq must be the same as the number of elements of f
	****************************************************************************************/
	T * buff_aeq = (T*)malloc(sizeof(T) * quadr.row * quadr.col);

	for (unsigned r = 0; r < quadr.row; ++r)
		for (unsigned c = 0; c < quadr.col; ++c)
			buff_aeq[(r * quadr.col) + c] = y[r][0];

	mxArray * aeq =  mxCreateDoubleMatrix(quadr.row, quadr.col, mxREAL);
	memcpy((void *)mxGetPr(aeq), (void *)buff_aeq, quadr.row * quadr.col * sizeof(T));
	delete buff_aeq;

	engPutVariable(ep, "Aeq", aeq);
	
	/***********************************************************************************
	* Create beq
	* Vector of doubles. Represents the constant vector in the constraints Aeq*x = beq.
	* The number of rows in Aeq must be the same as the number of elements of beq
	************************************************************************************/
	T * buff_beq = (T*)malloc(sizeof(T) * x.row);	
	for (unsigned c = 0; c < quadr.row; ++c)
		buff_beq[c] = zero;

	mxArray * beq =  mxCreateDoubleMatrix(1, quadr.row, mxREAL);

	memcpy((void *)mxGetPr(beq), (void *)buff_beq, quadr.row * sizeof(T));
	delete buff_beq;

	engPutVariable(ep, "beq", beq);
	
	/**************************************************************************************
	* Create lb and ub
	* Vectors of doubles. Represent the lower and upper bounds elementwise in lb ≤ x ≤ ub.
	***************************************************************************************/
	T * buff_lb = (T*)malloc(sizeof(T) * quadr.col);
	T * buff_ub = (T*)malloc(sizeof(T) * quadr.col);
	
	for (unsigned c = 0; c < quadr.col; ++c)
	{
		buff_lb[c] = zero;
		buff_ub[c] = 2000;//inf;
	}

	mxArray * lb =  mxCreateDoubleMatrix(1, quadr.col, mxREAL);
	memcpy((void *)mxGetPr(lb), (void *)buff_lb, quadr.col * sizeof(T));
	delete buff_lb;
	
	mxArray * ub =  mxCreateDoubleMatrix(1, quadr.col, mxREAL);
	memcpy((void *)mxGetPr(ub), (void *)buff_ub, quadr.col * sizeof(T));
	delete buff_ub;

	engPutVariable(ep, "lb", lb);
	engPutVariable(ep, "ub", ub);

	/*******************************************************************************
	* Creating x0
	* Vector of doubles. Optional. The initial point for some quadprog algorithms.
	*******************************************************************************/

	mxArray * x0 =  mxCreateDoubleMatrix(0, 0, mxREAL);
	engPutVariable(ep, "x0", x0);
	
    /*******************************************************
	* Creating options
	********************************************************/
	//engEvalString(ep, "options=optimset('Algorithm','interior-point-convex')");
	//engPutVariable(ep, "options", "optimset('Algorithm',interior-point-convex)");
	//if (debug)
		//_getch();
	 

	/*****************************************
	* We well read results from this
	******************************************/

	mxArray * result = NULL, * fval = NULL, *exitflag = NULL;

	engPutVariable(ep, "x", result);
	engPutVariable(ep, "fval", fval);
	engPutVariable(ep, "exitflag", exitflag);

	// Evaluate quadprog
	engEvalString(ep, "[x,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,optimset('Algorithm','interior-point-convex'))");
	
	if (debug)
		_getch();
	 
	vector<T> alpha(x.row, 0.0);
	unsigned sv_count = 0;
	if (debug)
		printf("\nRetrieving results...\n");
	if ((result = engGetVariable(ep,"x")) == NULL)
	    printf("Oops! You didn't create a variable 'x' in Matlab workspace.\n\n");
	else {
		if (debug && false)
		{
			printf("'x' is class %s\t\n", mxGetClassName(result));
			double * res = mxGetPr(result);
			cerr << *res << endl;
		}
		T * buff_res = (T*)malloc(sizeof(T) * x.row);	
		memcpy((void *)buff_res, (void *)mxGetPr(result), x.row * sizeof(T));

		for (unsigned i = 0; i < x.row; ++i)
		{
			alpha[i] = buff_res[i];
			if (abs(alpha[i]) > 0.001)
				++sv_count;
		}

		delete buff_res;

		if (debug)
			cerr << "result is:" << endl << alpha << endl;
	}

#if 0
	if ((fval = engGetVariable(ep,"fval")) == NULL)
	    printf("Oops! You didn't create a variable 'fval' in Matlab workspace.\n\n");
	else {
		printf("'fval' is class %s\t\n", mxGetClassName(fval));
		double * res = mxGetPr(fval);
		cerr << *res << endl;
	}
#endif

	mxDestroyArray(result);
	mxDestroyArray(fval);
	mxDestroyArray(exitflag);
	mxDestroyArray(h);
	mxDestroyArray(a);
	mxDestroyArray(b);
	mxDestroyArray(aeq);
	mxDestroyArray(beq);
	mxDestroyArray(lb);
	mxDestroyArray(ub);
	mxDestroyArray(x0);

	vector<T> w_sv(x.col, 0.0);
	unsigned sv = x.row;
	for (unsigned r = 0; r < x.row; ++r)
		for (unsigned c = 0; c < x.col; ++c)
		{
			if (abs(alpha[r]) > 0.001)
			{
				w_sv[c] += alpha[r] * y[r][0] * x[r][c];
				sv = r;
			}
		}

	w[0] = y[sv][0] - dot_product(w_sv, x[sv]);
	for (unsigned c = 0; c < x.col; ++c)
		w[c + 1] = w_sv[c];

	return sv_count;
}

} // namespace

#endif