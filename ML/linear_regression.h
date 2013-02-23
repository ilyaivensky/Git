#ifndef _LINEAR_REGRESSION_H
#define _LINEAR_REGRESSION_H

#include "matrix.h"

template <class T>
Matrix<T> linear_regression(const Matrix<T> & x, const Matrix<T> & y, T lambda = 0.0)
{
	Matrix<T> xtx = x.xtx();
	Matrix<T> reg = Matrix<T>::diag(xtx.row, lambda);
	xtx += reg;
	Matrix<T> inverse = xtx.invert();
	
	Matrix<T> h = inverse.multiply_by_transposed(x);
	Matrix<T> w = h * y;

	return w;
}

#endif