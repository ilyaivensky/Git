#ifndef _ERROR_H
#define _ERROR_H

#include "matrix.h"

extern double e;

namespace EF {

template <class T>
T classification(const T & h, const T & y)
{
	return (h * y >= 0.0 ? 0.0 : 1.0);
}

template <class T>
T square(const T & h, const T & y)
{
	return pow((h - y), 2);
}

template <class T>
T cross_entropy(const T & h, const T & y)
{
	return log(1 + pow(e, -y * h));
}

template <class T, class W>
T calc_error(const Matrix<T> & x, const W & w, const Matrix<T> & y, 
		   T (*error_f)(const T &, const T &))
{
	Matrix<T> h = x * w;
	T error = 0.0;
	for (unsigned r = 0; r < h.row; ++r)
		error += error_f(h[r][0], y[r][0]);

	error /= h.row;

	return error;
}

}

#endif