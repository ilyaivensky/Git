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
Matrix<T> Matrix<T>::minor(unsigned r, unsigned c) const
{
	if (r >= nrow() || c >= ncol())
		throw exception("minor: out of range");

	Matrix<T> m(nrow() - 1, ncol() - 1);
	for (unsigned i = 0, mi = 0; i < nrow(); ++i)
	{
		if (i == r) continue;

		const Row & row = (*this)[i];
		Row & mrow = m[mi];

		for (unsigned j = 0, mj = 0; j < ncol(); ++j)
		{
			if (j == c) continue;

			mrow[mj] = row[j];
			++mj;
		}
		++mi;
	}

	return m;
}

template <class T>
Matrix<T>::Matrix(unsigned ncol, const vector<T> & col)
{
	reserve(ncol);
	for (const auto & c : col)
		push_back(Row(ncol, c));
}

template <class T>
void Matrix<T>::resize(unsigned row, unsigned col)
{
	vector<Row>::resize(row, Row(col));
	for (auto & row : *this)
		row.resize(col);
}

template <class T>
void Matrix<T>::add_row(const Row & x)
{
	if (!empty() && front().size() != x.size())
		throw exception("Matrix: inconsistent num of columns");

	push_back(x);
}

template <class T>
Matrix<T> Matrix<T>::diag(unsigned n, const T & val)
{
	Matrix<T> r(n);
	for (unsigned i = 0; i < n; ++i)
		r[i][i] = val;

	return r;
}

template <class T>
void Matrix<T>::scale(const T & lb, const T & ub)
{
	if (ub <= lb)
		throw exception("upper bound is not greater than lower bound");

	vector<T> feature_max(ncol(), numeric_limits<T>::min());
	vector<T> feature_min(ncol(), numeric_limits<T>::max());

	for (const auto & row : *this)
		for (unsigned c = 0; c < ncol(); ++c)
		{
			feature_max[c] = std::max(row[c], feature_max[c]);
			feature_min[c] = std::min(row[c], feature_min[c]);
		}

	for (auto & row : *this)
		for (unsigned c = 0; c < ncol(); ++c)
		{
			if (feature_max[c] == feature_min[c])
				continue;
			if (row[c] == feature_min[c])
				row[c] = lb;
			else if (row[c] == feature_max[c])
				row[c] = ub;
			else
				row[c] = lb + (ub - lb) * (row[c] - feature_min[c]) / (feature_max[c] - feature_min[c]);
		}
}

template <class T>
void Matrix<T>::random_init() 
{
	for (auto & row : *this)
		row = random_example<T>(ncol());
}

template <class T>
void Matrix<T>::random_init_0() 
{
	for (auto & row : *this)
		row = random_example_0<T>(ncol());
}

template <class T>
Matrix<T> Matrix<T>::get_transformed(vector<T> (*t)(const vector<T> &)) const
{
	if (*t == 0)
		return *this;

	Matrix transformed(nrow(), 0);
	auto itTr = transformed.begin();
	for (auto it = begin(), itEnd = end(); it != itEnd; ++it, ++itTr)
		*itTr = t(*it);

	return transformed;
}

template <class T>
void Matrix<T>::transform_self(vector<T> (*t)(const vector<T> &))
{
	if (*t == 0)
		return;

	for (auto & row : *this)
		row = t(row);
}

// X(m,n), Y(m,n) - treated as transposed
template <class T>
Matrix<T> Matrix<T>::multiply_by_transposed(const Matrix<T> & other) const
{
	Matrix<T> res(nrow(), other.nrow());
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
	for (auto & row : *this)
		row *= t; 

	return *this;
}

template <class T>
Matrix<T> & Matrix<T>::operator /= (const T & t)
{
	for (auto & row : *this)
		row /= t;
		
	return *this;
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

	auto rit1 = m1.begin(), rit1End = m1.end();
	auto rit2 = m2.begin();

	for (; rit1 != rit1End; ++rit1, ++rit2)
		*rit1 += *rit2;

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

	auto rit1 = m1.begin(), rit1End = m1.end();
	auto rit2 = m2.begin();

	for (; rit1 != rit1End; ++rit1, ++rit2)
		*rit1 -= *rit2;

	return m1;
}

template <class T> 
ostream & operator << (ostream & os, const Matrix<T> & m)
{
	for (const auto & row : m)
		os << row << endl;

	return os;
}

template <class T>
Matrix<T> operator ^ (const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("Incompatible vectors"); 

	Matrix<T> res;

	for (auto itV1 = v1.begin(), itV2 = v2.begin(), 
		itVEnd = v1.end(); itV1 != itVEnd; ++itV1, ++itV2)
	{
		res.add_row({ *itV1, *itV2 });
	}

	return res;
}

template <class T>
Matrix<T> operator ^ (const Matrix<T> & m, const vector<T> & v)
{
	if (m.size() != v.size())
		throw exception("Incompatible sizes"); 

	Matrix<T> res(m);

	auto itV = v.begin();
	for (auto itM = m.begin(), itMEnd = m.end(); itM != itMEnd; ++itM, ++itV)
		itM->push_back(*itV);

	return res;
}

template <class T>
Matrix<T> operator ^ (const vector<T> & v, const Matrix<T> & m)
{
	if (m.size() != v.size())
		throw exception("Incompatible sizes"); 

	Matrix<T> res(v.size(), 0);
	auto itV = v.begin();
	auto itM = m.begin();
	for (auto itRes = res.begin(), itResEnd = res.end(); itRes != itResEnd; ++itRes, ++itV, ++itM)
	{
		itRes->reserve(m.ncol() + 1);
		itRes->push_back(*itV);
		itRes->insert(itRes->end(), itM->begin(), itM->end()); 
	}

	return res;
}
