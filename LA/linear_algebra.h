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
#include <tuple>

// Returns trace of matrix
template <class T>
T tr(const Matrix<T> & m);

// Returns determinant of matrix
// Implements Laplace expansion 
template <class T>
T det(const Matrix<T> & m);

/** 
 Returns m^(-1)
 Implements Gauss-Jordan elimination
 */
template <class T>
Matrix<T> inv(const Matrix<T> & m);

// Returns Gramian matrix (x^t * x)
template <class T>
Matrix<T> gram(const Matrix<T> & x);

// Returns vector of column means
template <class T>
vector<T> mean_col(const Matrix<T> & m);

/** 
 Returns deviation scores matrix, where each element of
 column vector is deviation score of corresponding 
 element in x from the mean of its column in x
 */
template <class T>
Matrix<T> dev(const Matrix<T> & x);

// Returns covariance matrix
template <class T>
Matrix<T> cov(const Matrix<T> & m);

template <class T>
bool is_square(const Matrix<T> & m) { return m.nrow() == m.ncol(); }

template <class T>
bool is_singular(const Matrix<T> & m) { return !is_square(m) || det(m) == 0; }

/**
 Factorizes matrix A into L (lower triangular), U (upper triangular) and P (permutation matrix)
 such as PA = LU
 */
template <class T>
tuple<Matrix<T>, Matrix<T>, Matrix<T>> lu(const Matrix<T> & m);

/**
 Returns row reduced echelon form (RREF) of m
 */
template <class T>
Matrix<T> rref(const Matrix<T> & m);

template <class T> 
Matrix<T> linear_solution(const Matrix<T> & a, const Matrix<T> & y) { return inv(a) * y; } 

/** 
 Solves quadratic equation,
 parameter is a vector of coefficients
 a[0] * x^2 + a[1] * x + a[2]
 */
template <class T>
vector<T> solve_quadratic(const vector<T> & polynom);

/**
 Returns coefficients of characteristic polynomial
 (e.g., vector 'a' is coefficients of polynomial 
  a[0] * x^3 + a[1] * x^2 + a[2] * x + a[3] = 0)
 Supports only 2x2 and 3x3 matrices
 */
template <class T>
vector<T> characteristic_polynomial(const Matrix<T> & m);

/** Returns eigenvalues of 2x2 matrix */
template <class T>
vector<T> eigenvalues_2x2(const Matrix<T> & m);

/**
 Returns vector of pairs <eigenval, eigenvec> for given 2x2 matrix
 */
template <class T>
vector<pair<T, vector<T>>> eigen_2x2(const Matrix<T> & m);

#include "LA/linear_algebra_impl.h"

#endif
