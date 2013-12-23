/*                                                                 -*- C++ -*-
 * File: matrix.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Mar 1, 2013
 *
 * Description:
 *   Implementation of matrix
 *   
 */

#ifndef _MATRIX_H_
#define	_MATRIX_H_

#include <vector>
#include <initializer_list>

using namespace std;

template <class T>
class Matrix : public vector<vector<T> > 
{
public:
	typedef vector<T> Row;

	typedef vector<Row> Super;

	Matrix() {};

	Matrix(const initializer_list<Row> & list) : 
		Super(list) {};
	
	// Creates matrix row * col and initializes with val
	Matrix(unsigned row, unsigned col, T val = T()) : 
		vector<Row>(row, Row(col, val)) {}

	// Creates square matrix n * n and initializes with default value of T
 	Matrix(unsigned n) : 
		vector<Row>(n, Row(n)) {}

	// Creates outer product of v1 and v2
	Matrix(const vector<T> & v1, const vector<T> & v2);

	// Creates matrix with identical rows 
	Matrix(const vector<T> & row, unsigned nrow) :
		vector<Row>(nrow, row) {}

	// Creates matrix with identical columns
	Matrix(unsigned ncol, const vector<T> & col);

	unsigned nrow() const { return size(); }
	unsigned ncol() const { return size() ? front().size() : 0; }

	// Random initialization 
	// Include "ML/random.h" before including this file (matrix.h)
	// Include "ML/random.h" is not needed if those functions are not used
    void random_init();
	void random_init_0();

	void interactive_init();

	void resize(unsigned row, unsigned col)
	{
		vector<Row>::resize(row, Row(col));
		for (auto it = begin(), itEnd = end(); it != itEnd; ++it)
			it->resize(col);
	}

	void add_row(const Row & x)
	{
		if (!empty() && front().size() != x.size())
			throw exception("Matrix: inconsistent num of columns");

		push_back(x);
	}
	
	// Returns (*this) * (other)^t
	Matrix multiply_by_transposed(const Matrix & other) const;

	// Creates transformed variant of self according to transform function
	Matrix get_transformed(vector<T> (*t)(const vector<T> &)) const;

	// Transforms according to transform function
	void transform_self(vector<T> (*t)(const vector<T> &));

	// Creates binary matrix of type B according 
	// to the label function (e.g. 1:0 or 1:-1 or true:false) 
	template <class B>
	Matrix<B> get_binary(B (*label)(const T &)) const
	{
		Matrix<B> res(nrow(), ncol());
			
		const_iterator itRow = begin(), itRowEnd = end();
		auto resRowIt = res.begin();

		for (; itRow != itRowEnd; ++itRow, ++resRowIt)
		{
			auto itCol = itRow->begin(), itColEnd = itRow->end();
			auto resColIt = resRowIt->begin();

			for (; itCol != itColEnd; ++itCol, ++resColIt)
				*resColIt = label(*itCol);
		}

		return res;
	}

	// Scales the content of matrix 
	// into range [lower_bound, upper_bound]
	void scale(T lower_bound, T upper_bound);

	// Creates diagonal matrix 
	static Matrix diag(unsigned size, const T & value); 

	// Returns minor version of m (removed row r and col c)
	Matrix minor(unsigned r, unsigned c) const;

	Matrix operator * (const Matrix & m) const;
	Matrix operator * (const vector<T> & v) const;

	// Multiply each element by scalar
	Matrix & operator *= (const T & t);
	Matrix & operator /= (const T & t);

	Matrix operator + (const Matrix & m) const;
	Matrix & operator += (const Matrix & m);

	Matrix operator - (const Matrix & m) const;
	Matrix & operator -= (const Matrix & m);
};

template <class T>
ostream & operator << (ostream & os, const Matrix<T> & m);

// Creates bi-vector
template <class T>
Matrix<T> operator ^ (const vector<T> & v1, const vector<T> & v2);

template <class T>
Matrix<T> operator ^ (const Matrix<T> & m, const vector<T> & v);

template <class T>
Matrix<T> operator ^ (const vector<T> & v, const Matrix<T> & m);

template <class T>
Matrix<T> operator ^ (const Matrix<T> & m1, const Matrix<T> & m2);


#include "matrix_impl.h"

#endif
