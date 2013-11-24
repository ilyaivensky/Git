/*                                                                 -*- C++ -*-
 * File: matrix.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Mar 1, 2013
 *
 * Description:
 *   Implementation of matrix and some operations
 *   
 */

#ifndef _MATRIX_H_
#define	_MATRIX_H_

#include <vector>
#include <limits>

using namespace std;

template <class T>
struct Matrix : public vector<vector<T> > 
{
	typedef vector<T> Row;
			
	// Creates matrix row * col and initializes with 0
	Matrix(unsigned row, unsigned col) : 
		vector<Row>(row, Row(col, 0)), row(row), col(col) {}

	// Creates square matrix n * n and initializes with 0.0
 	Matrix(unsigned n) : 
		vector<Row>(n, Row(n, 0)), row(n), col(n) {}

	// Creates outer product of v1 and v2
	Matrix(const vector<T> & v1, const vector<T> & v2);

	// Creates empty matrix
	Matrix() : row(0), col(0) {}

	// Random initialization 
	// Include "ML/random.h" before including this file (matrix.h)
	// Include "ML/random.h" is not needed if those functions are not used
    void random_init();
	void random_init_0();

	void interactive_init();

	bool empty() const { return row == 0; } 

	void resize(unsigned resized_row, unsigned resized_col)
	{
		vector<Row>::resize(resized_row, Row(col));
		row = resized_row;
		for (vector<Row>::iterator it = begin(), itEnd = end(); it != itEnd; ++it)
			it->resize(resized_col);
		col = resized_col;
	}

	void add_row(const Row & x)
	{
		push_back(x);
		++row;
	}
	
	// throws exception if not square 
	Matrix invert() const;

	Matrix xtx() const;
	T determinant() const;
	T trace() const;

	// The same as above, but we transpose other
	Matrix multiply_by_transposed(const Matrix & other) const;

	// Creates transformed variant of self according to transform function
	Matrix get_transformed(vector<T> (*t)(const vector<T> &)) const;

	// Transforms according to transform function
	void transform_self(vector<T> (*t)(const vector<T> &));

	// Creates binary matrix according to the label function 
	// (e.g. 1:0 or 1:-1) for positive and negative values from 'm'
	static Matrix binary(const Matrix & m, signed (*label)(const T &));

	void scale(T lb, T ub);

	// Creates diagonal matrix 
	static Matrix<T> diag(unsigned size, T value); 

	bool is_square() const { return row == col; }

	unsigned row;
	unsigned col;
};

template <class T>
ostream & operator<<(ostream & os, const Matrix<T> & m);

template <class T>
Matrix<T> operator *(const Matrix<T> & m, const vector<T> & v);

template <class T>
Matrix<T> operator *(const Matrix<T> & m, const Matrix<T> & v);

template <class T>
Matrix<T> & operator+=(Matrix<T> & m1, const Matrix<T> & m2);

template <class T>
bool operator ==(const Matrix<T> & m1, const Matrix<T> & m2);

template <class T>
bool operator !=(const Matrix<T> & m1, const Matrix<T> & m2);

#include "matrix_impl.h"

#endif