#ifndef _LINEAR_SOULTIONS_H
#define _LINEAR_SOLUTIONS_H

#include "matrix.h"

/**********************************************************
* Pseudoinverse solution for linear weights
* (implements left inverse)
***********************************************************/
template <class T>
Matrix<T> linear_pseudoinverse_solution(const Matrix<T> & x, const Matrix<T> & y, T lambda = 0.0)
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