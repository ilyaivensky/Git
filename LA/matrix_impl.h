/*                                                                 -*- C++ -*-
 * File: matrix_impl.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Mar 18, 2013
 */

#include "vector_utils.h"

template <class T>
Matrix<T>::Matrix(const vector<T> & v1, const vector<T> & v2) :
	vector<Row>(v1.size(), Row(v2.size(), 0.0)), row(row), col(col)
{
	if (v1.size() != v2.size())
		throw exception("outer_product: v1.size() != v2.size()");

	for (unsigned i = 0; i < v1.size(); ++i)
	{
		Row & r = (*this)[i];
		for (unsigned j = 0; j < v2.size(); ++j)
			r[j] = v1[i] * v2[j];
	}
}

template <class T>
Matrix<T> Matrix<T>::diag(unsigned n, T lambda)
{
	Matrix<T> r(n);
	for (unsigned i = 0; i < n; ++i)
		r[i][i] = lambda;

	return r;
}

template <class T>
void Matrix<T>::scale(T lb, T ub)
{
	if (ub <= lb)
		throw exception("upper bound is not greater than lower bound");

	vector<T> feature_max(col_, numeric_limits<T>::min());
	vector<T> feature_min(col_, numeric_limits<T>::max());

	for (unsigned r = 0; r < row; ++r)
		for (unsigned c = 0; c < col; ++c)
		{
			const Row & row = (*this)[r];
			feature_max[c] = std::max(row[c], feature_max[c]);
			feature_min[c] = std::min(row[c], feature_min[c]);
		}

	for (unsigned r = 0; r < row; ++r)
	{
		const Row & row = (*this)[r];
		for (unsigned c = 0; c < col; ++c)
		{
			if (feature_max[c] == feature_min[c])
				continue;
			if (row[c] == feature_min[c])
				row[c] = lb;
			else if (row[c] == feature_max[c])
				row[c] = ub;
			else
				row[c] = lb + (ub - lb) * (row[c]  - feature_min[c]) / (feature_max[c] - feature_min[c]);
		}
	}
}

template <class T>
void Matrix<T>::random_init() 
{
	for (iterator it = begin(), itEnd = end(); it != itEnd; ++it)
		(*it) = random_example<T>(col_);
}

template <class T>
void Matrix<T>::random_init_0() 
{
	for (iterator it = begin(), itEnd = end(); it != itEnd; ++it)
		(*it) = random_example_0<T>(col_);
}

template <class T>
Matrix<T> Matrix<T>::get_transformed(vector<T> (*t)(const vector<T> &)) const
{
	if (*t == 0)
		return *this;

	Matrix transformed(this->row, 0);
	iterator itTr = transformed.begin();
	for (const_iterator it = begin(), itEnd = end(); it != itEnd; ++it, ++itTr)
		*itTr = t(*it);

	transformed.col_ = transformed[0].size();

	return transformed;
}

template <class T>
void Matrix<T>::transform_self(vector<T> (*t)(const vector<T> &))
{
	if (*t == 0)
		return;

	for (vector<vector<T> >::iterator it = this->begin(), itEnd = this->end(); it != itEnd; ++it)
		*it = t(*it);

	col_ = (*this)[0].size();
}

// X(m,n), Y(m,n) - treated as transposed
template <class T>
Matrix<T> Matrix<T>::multiply_by_transposed(const Matrix<T> & other) const
{
	Matrix<T> res(this->nrow(), other.nrow());
	for (unsigned m = 0; m < this->nrow(); ++m)
	{
		for (unsigned n = 0; n < other.nrow(); ++n) // we transpose other
		{
			T elem = 0;
			for (unsigned k = 0; k < other.ncol(); ++k)
			{
				elem += (*this)[m][k] * other[n][k]; 
			}
			res[m][n] = elem;
		}
	}

#if 0
	cerr << "multiply_by_transposed is:" << endl << res << endl;
#endif

	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const Matrix<T> & m2) const
{
	const Matrix<T> & m1 = *this;

	Matrix<T> res(m1.nrow(), m2.ncol());
	for (unsigned m = 0; m < m1.nrow(); ++m)
	{
		for (unsigned n = 0; n < m2.ncol(); ++n) 
		{
			T & elem = res[m][n];
			for (unsigned k = 0; k < m2.nrow(); ++k)
				elem += m1[m][k] * m2[k][n]; 
		}
	}

	return res;
}

template <class T>
Matrix<T> Matrix<T>::operator * (const vector<T> & v) const
{
	Matrix<T> res(this->nrow(), 1);
	for (unsigned r = 0; r < this->nrow(); ++r)
		res[r][0] = inner_product((*this)[r], v);

	return res;
}

template <class T>
Matrix<T> & Matrix<T>::operator *= (const T & t)
{
	for (Matrix<T>::iterator itRow = this->begin(), itRowEnd = this->end(); itRow != itRowEnd; ++itRow)
		for (Matrix<T>::Row::iterator itCol = itRow->begin(), itColEnd = itRow->end(); itCol != itColEnd; ++itCol)
			(*itCol) *= t;

	return *this;
}

template <class T>
Matrix<T> & Matrix<T>::operator /= (const T & t)
{
	return operator *= (static_cast<T>(1) / t);
}

template <class T>
Matrix<T> Matrix<T>::operator + (const Matrix<T> & m2) const
{
	Matrix<T> m(*this);
	return m += m2;
}

template <class T>
Matrix<T> & Matrix<T>::operator += (const Matrix<T> & m2)
{
	Matrix & m1 = *this;

	if (m1.nrow() != m2.nrow() || m1.ncol() != m2.ncol())
		throw exception("not compatible for operator '+='");

	Matrix<T>::iterator rit1 = m1.begin(), rit1End = m1.end();
	Matrix<T>::const_iterator rit2 = m2.begin();

	for (; rit1 != rit1End; ++rit1, ++rit2)
	{
		Matrix<T>::Row::iterator cit1 = rit1->begin(), cit1End = rit1->end();
		Matrix<T>::Row::const_iterator cit2 = rit2->begin();

		for (; cit1 != cit1End; ++cit1, ++cit2)
			(*cit1) += (*cit2);
	}

	return m1;
}


template <class T>
Matrix<T> Matrix<T>::operator - (const Matrix<T> & m2) const
{
	Matrix<T> m(*this);
	return m -= m2; 
}

template <class T>
Matrix<T> & Matrix<T>::operator -= (const Matrix<T> & m2)
{
	Matrix<T> & m1 = *this;

	if (m1.nrow() != m2.nrow() || m1.ncol() != m2.ncol())
		throw exception("not compatible for operator '-='");

	Matrix<T>::iterator rit1 = m1.begin(), rit1End = m1.end();
	Matrix<T>::const_iterator rit2 = m2.begin();

	for (; rit1 != rit1End; ++rit1, ++rit2)
	{
		Matrix<T>::Row::iterator cit1 = rit1->begin(), cit1End = rit1->end();
		Matrix<T>::Row::const_iterator cit2 = rit2->begin();

		for (; cit1 != cit1End; ++cit1, ++cit2)
			(*cit1) -= (*cit2);
	}

	return m1;
}

template <class T> 
ostream & operator<<(ostream & os, const Matrix<T> & m)
{
	for (Matrix<T>::const_iterator it = m.begin(), itEnd = m.end(); it != itEnd; ++it)
		os << *it << endl;
     
	return os;
}
