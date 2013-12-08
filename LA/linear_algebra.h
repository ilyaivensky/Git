#ifndef _LINEAR_ALGEBRA_H_
#define _LINEAR_ALGEBRA_H_

/*                                                                 -*- C++ -*-
 * File: linear_algebra.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 7, 2013
 *   
 */

#include "LA/matrix.h"

// Returns trace of matrix
template <class T>
T tr(const Matrix<T> & m);

// Returns determinant of matrix
template <class T>
T det(const Matrix<T> & m);

// Returns m^(-1)
template <class T>
Matrix<T> inv(const Matrix<T> & m);

// Returns Gramian matrix (x^t * x)
template <class T>
Matrix<T> gram(const Matrix<T> & x);

// Returns covariance matrix
template <class T>
Matrix<T> cov(const Matrix<T> & m);

#include "LA/linear_algebra_impl.h"

#endif
