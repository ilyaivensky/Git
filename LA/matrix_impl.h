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
		for (unsigned j = 0; j < v2.size(); ++j)
			(*this)[i][j] = v1[i] * v2[j];
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
Matrix<T> Matrix<T>::binary(const Matrix<T> & m, signed (*label)(const T &))
{
	Matrix<T> res(m.row, m.col);
	for (unsigned r = 0; r < m.row; ++r)
		for (unsigned c = 0; c < m.col; ++c)
			res[r][c] = label(m[r][c]);

	return res;
}

template <class T>
void Matrix<T>::scale(T lb, T ub)
{
	if (ub <= lb)
		throw exception("upper bound is not greater than lower bound");

	vector<T> feature_max(col, numeric_limits<T>::min());
	vector<T> feature_min(col, numeric_limits<T>::max());

	for (unsigned r = 0; r < row; ++r)
		for (unsigned c = 0; c < col; ++c)
		{
			feature_max[c] = std::max((*this)[r][c], feature_max[c]);
			feature_min[c] = std::min((*this)[r][c], feature_min[c]);
		}

	for (unsigned r = 0; r < row; ++r)
		for (unsigned c = 0; c < col; ++c)
		{
			if (feature_max[c] == feature_min[c])
				continue;
			if ((*this)[r][c] == feature_min[c])
				(*this)[r][c] = lb;
			else if ((*this)[r][c] == feature_max[c])
				(*this)[r][c] = ub;
			else
				(*this)[r][c] = lb + (ub - lb) * ((*this)[r][c]  - feature_min[c]) / (feature_max[c] - feature_min[c]);
		}
}

template <class T>
void Matrix<T>::random_init() 
{
	unsigned TRAINING_SET = row;
	unsigned DIM = col;

	for (unsigned m = 0; m < TRAINING_SET; ++m)
		(*this)[m] = random_example<T>(DIM);

#if 0
	cerr << "Training set is initialized as:" << endl << *this << endl;
#endif
}

template <class T>
void Matrix<T>::random_init_0() 
{
	unsigned TRAINING_SET = row;
	unsigned DIM = col;

	for (unsigned m = 0; m < TRAINING_SET; ++m)
		(*this)[m] = random_example_0<T>(DIM);

#if 0
	cerr << "Training set is initialized as:" << endl << *this << endl;
#endif
}

template <class T>
T Matrix<T>::determinant() const
{
	if (!is_square())
		throw exception("determinant is not defined for non-square matrix"); 

	T det = 0;
	for (unsigned i = 0; i < row; ++i)
	{
		for (unsigned r = 0, c = i; r < row; ++r, ++c)
		{
			if (c == row) c = 0;
			det += (*this)[r][c];
		}

		for (unsigned r = 0; c = row - i; r < row; ++r, --c)
		{
			if (c == 0) c = row;
			det -= (*this)[r][c - 1];
		}
	}
	return det;
}

template <class T>
T Matrix<T>::trace() const
{
	if (!is_square())
		throw exception("trace is not defined for non-square matrix");

	T tr = 0;
	for (unsigned i = 0; i < row; ++i)
		tr += (*this)[i][i];

	return tr;
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

	transformed.col = transformed[0].size();

	return transformed;
}

template <class T>
void Matrix<T>::transform_self(vector<T> (*t)(const vector<T> &))
{
	if (*t == 0)
		return;

	for (vector<vector<T> >::iterator it = this->begin(), itEnd = this->end(); it != itEnd; ++it)
		*it = t(*it);

	col = (*this)[0].size();
}

template <class T>
Matrix<T> Matrix<T>::invert() const
{
	if (!is_square())
		throw exception("cannot inverse non-square matrix\n");

	Matrix<T> matrix(*this);
	unsigned n = matrix.size();
	Matrix<T> inverse(n);

	// Init inverse as ident matrix
	for (unsigned i = 0; i < n; ++i)
		inverse[i][i] = 1.0;

	for (unsigned i = 0; i < n; ++i)
	{
		for (unsigned j = 0; j < n; ++j)
		{
			if (i != j)
			{
				if (matrix[i][i] == 0.0)
					throw exception("matrix cannot be inverted\n");
				T ratio = matrix[j][i]/matrix[i][i];
				for (unsigned k = 0; k < n; ++k)
				{
					matrix[j][k] -= ratio * matrix[i][k];
					inverse[j][k] -= ratio * inverse[i][k];
				}
			}
		}
	}

	for (unsigned i = 0; i < n; ++i)
	{
		T a = matrix[i][i];
		if (a == 0.0)
			throw exception("matrix cannot be inverted\n");
		for (unsigned j = 0; j < n; ++j)
		{
			matrix[i][j] /= a;
			inverse[i][j] /= a;
		}
	}

#if 0
	cerr << "The inverse matrix is: " << endl << inverse << endl;
#endif
	return inverse;
}

template <class T>
Matrix<T> Matrix<T>::xtx() const
{
	unsigned DIM = col;

	Matrix<T> res(DIM);
	for (unsigned n = 0; n < DIM; ++n)
	{  
		Row & res_row = res[n];
		// We start from col n because xtx is symmetric (i.e. res[n][q] == res[q][n])
		for (unsigned q = n; q < DIM; ++q)
		{
			T & res_n_q = res_row[q];
			for (const_iterator it = begin(), itEnd = end(); it != itEnd; ++it)
			{
				const Row & src_row = *it; 
				res_n_q += src_row[n] * src_row[q];
			}
			// Copy res[n][q] to res[q][n]
			res[q][n] = res_n_q;
		}
	}

#if 0
	cerr << "XTX is:" << endl << res << endl;
#endif

	return res;
}

template <class T>
Matrix<T> & operator+=(Matrix<T> & m1, const Matrix<T> & m2)
{
	if (m1.row != m2.row || m1.col != m2.col)
		throw exception("not compatible for operator +=");

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
Matrix<T> operator *(const Matrix<T> & m1, const Matrix<T> & m2)
{
	Matrix<T> res(m1.row, m2.col);
	for (unsigned m = 0; m < m1.row; ++m)
	{
		for (unsigned n = 0; n < m2.col; ++n) 
		{
			T & elem = res[m][n];
			for (unsigned k = 0; k < m2.row; ++k)
				elem += m1[m][k] * m2[k][n]; 
		}
	}

	return res;
}

template <class T>
Matrix<T> operator *(const Matrix<T> & m, const vector<T> & v)
{
	Matrix<T> res(m.row, 1);
	for (unsigned r = 0; r < m.row; ++r)
		res[r][0] = inner_product(m[r], v);

	return res;
}

template <class T>
bool operator ==(const Matrix<T> & m1, const Matrix<T> & m2)
{
	if (m1.row != m2.row || m1.col != m2.col)
		return false;

	for (unsigned r = 0; r < m1.row; ++r)
		if (m1[r] != m2[r])
			return false;

	return true;
}

template <class T>
bool operator !=(const Matrix<T> & m1, const Matrix<T> & m2)
{
	return !operator==(m1, m2);
}
// X(m,n), Y(m,n) - treated as transposed
template <class T>
Matrix<T> Matrix<T>::multiply_by_transposed(const Matrix<T> & other) const
{
	Matrix<T> res(this->row, other.row);
	for (unsigned m = 0; m < this->row; ++m)
	{
		for (unsigned n = 0; n < other.row; ++n) // we transpose other
		{
			T elem = 0.0;
			for (unsigned k = 0; k < other.col; ++k)
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
ostream & operator<<(ostream & os, const Matrix<T> & m)
{
	for (unsigned i = 0; i < m.row; ++i)
		os << m[i] << endl;
     
	return os;
}
