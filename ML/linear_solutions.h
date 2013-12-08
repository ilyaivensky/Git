#ifndef _LINEAR_SOULTIONS_H
#define _LINEAR_SOLUTIONS_H

#include "LA/matrix.h"
#include "LA/linear_algebra.h"

/**********************************************************
* Pseudoinverse solution for linear weights
* (implements left inverse)
***********************************************************/
template <class T>
Matrix<T> linear_pseudoinverse_solution(const Matrix<T> & x, const Matrix<T> & y, T regularization = 0.0)
{
	Matrix<T> xtx = gram(x);
	if (regularization != 0.0)
	{
		Matrix<T> reg = Matrix<T>::diag(xtx.nrow(), regularization);
		xtx += reg;
	}
	Matrix<T> inverted = inv(xtx);
	
	Matrix<T> h = inverted.multiply_by_transposed(x);
	Matrix<T> w = h * y;

	return w;
}

#endif