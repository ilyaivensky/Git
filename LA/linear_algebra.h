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
// Implements Laplace expansion 
template <class T>
T det(const Matrix<T> & m);

// Returns m^(-1)
template <class T>
Matrix<T> inv(const Matrix<T> & m);

// Returns Gramian matrix (x^t * x)
template <class T>
Matrix<T> gram(const Matrix<T> & x);

// Returns vector of column medians
template <class T>
vector<T> medc(const Matrix<T> & m);

// Returns deviation scores matrix, where each element of
// column vector is deviation score of corresponding 
// element in x from the mean of its column in x
template <class T>
Matrix<T> dev(const Matrix<T> & x);

// Returns covariance matrix
template <class T>
Matrix<T> cov(const Matrix<T> & m);

template <class T>
bool is_square(const Matrix<T> & m) { return m.nrow() == m.ncol(); }

template <class T>
bool is_singular(const Matrix<T> & m) { return !is_square(m) || det(m) == 0; }

template <class T> 
Matrix<T> linear_solution(const Matrix<T> & a, const Matrix<T> & y) { return inv(a) * y; } 

#include "LA/linear_algebra_impl.h"

#endif
