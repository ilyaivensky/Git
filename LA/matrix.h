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
#include <limits>

using namespace std;

template <class T>
class Matrix : public vector<vector<T> > 
{
public:
	typedef vector<T> Row;

	// Creates matrix row * col and initializes with val
	Matrix(unsigned row, unsigned col, T val = T()) : 
		vector<Row>(row, Row(col, val)), row_(row), col_(col) {}

	// Creates square matrix n * n and initializes with default value of T
 	Matrix(unsigned n) : 
		vector<Row>(n, Row(n)), row_(n), col_(n) {}

	// Creates outer product of v1 and v2
	Matrix(const vector<T> & v1, const vector<T> & v2);

	// Creates empty matrix
	Matrix() : row_(0), col_(0) {}

	unsigned nrow() const { return row_; }
	unsigned ncol() const { return col_; }

	// Random initialization 
	// Include "ML/random.h" before including this file (matrix.h)
	// Include "ML/random.h" is not needed if those functions are not used
    void random_init();
	void random_init_0();

	void interactive_init();

	bool empty() const { return row_ == 0; } 

	void resize(unsigned row, unsigned col)
	{
		vector<Row>::resize(row, Row(col));
		row_ = row;
		for (vector<Row>::iterator it = begin(), itEnd = end(); it != itEnd; ++it)
			it->resize(col);
		col_ = col;
	}

	void add_row(const Row & x)
	{
		if (!row_) col_ = x.size();
		else if (x.size() != col_)
			throw exception("Matrix: inconsistent num of columns");

		push_back(x);
		++row_;
	}
	
	// The same as above, but we transpose other
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
		Matrix<B> res(row_, col_);
			
		const_iterator itRow = begin(), itRowEnd = end();
		Matrix<B>::iterator resRowIt = res.begin();

		for (; itRow != itRowEnd; ++itRow, ++resRowIt)
		{
			Row::const_iterator itCol = itRow->begin(), itColEnd = itRow->end();
			Matrix<B>::Row::iterator resColIt = resRowIt->begin();

			for (; itCol != itColEnd; ++itCol, ++resColIt)
				*resColIt = label(*itCol);
		}

		return res;
	}

	// Scales the content of matrix 
	// into range [lower_bound, upper_bound]
	void scale(T lower_bound, T upper_bound);

	// Creates diagonal matrix 
	static Matrix diag(unsigned size, T value); 

	bool is_square() const { return row_ == col_; }

	Matrix operator * (const Matrix & m) const;
	Matrix operator * (const vector<T> & v) const;

	// Multiply each element by scalar
	Matrix & operator *= (const T & t);
	Matrix & operator /= (const T & t);

	Matrix operator + (const Matrix & m) const;
	Matrix & operator += (const Matrix & m);

	Matrix operator - (const Matrix & m) const;
	Matrix & operator -= (const Matrix & m);

private:
	unsigned row_;
	unsigned col_;
};

template <class T>
ostream & operator << (ostream & os, const Matrix<T> & m);

#include "matrix_impl.h"

#endif
