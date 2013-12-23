/*                                                                 -*- C++ -*-
 * File: linear_algebra_impl.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 7, 2013
 *   
 */

#include <limits>

template <class T>
T tr(const Matrix<T> & m)
{
	if (!is_square(m))
		throw exception("trace is not defined for non-square matrix");

	T tr = 0;
	for (unsigned i = 0; i < m.nrow(); ++i)
		tr += m[i][i];

	return tr;
}

template <class T>
T det(const Matrix<T> & m)
{
	if (!is_square(m))
		throw exception("determinant is not defined for non-square matrix"); 

	switch (m.nrow())
	{
	case 0:
		return 0;
	case 1:
		return m[0][0];
	case 2:
		return m[0][0] * m[1][1] - m[0][1] * m[1][0];
	default:
		break;
	}

	T d = 0;
	for (unsigned j = 0; j < m.nrow(); ++j)
		d += static_cast<T>(pow(-1, j + 2)) * m[0][j] * det(m.minor(0, j));

	return d;
}
		
template <class T>
Matrix<T> gram(const Matrix<T> & m)
{
	unsigned DIM = m.ncol();

	Matrix<T> res(DIM, DIM, 0);	
	for (unsigned i = 0; i < DIM; ++i)
	{
		Matrix<T>::Row & ri = res[i];
		for (unsigned j = i; j < DIM; ++j)
		{
			T & rij = ri[j]; 

			for (auto row : m)
				rij += row[i] * row[j];

			res[j][i] = rij;
		}
	}

	return res;
}

template <class T>
vector<T> mean_col(const Matrix<T> & m)
{
	Matrix<T> k(1, m.nrow(), 1);
	Matrix<T> res(k * m);
	res /= static_cast<T>(m.nrow());
	return res[0];
}

template <class T>
Matrix<T> dev(const Matrix<T> & m)
{
	return m - Matrix<T>(mean_col(m), m.nrow());
}

template <class T>
Matrix<T> cov(const Matrix<T> & m)
{
	Matrix<T> deviation = dev(m);
	Matrix<T> cov = gram(deviation);
	cov /= static_cast<T>(m.nrow());
	return cov;
}

template <class T>
Matrix<T> inv(const Matrix<T> & m)
{
	if (is_singular(m))
		throw exception("cannot inverse singular matrix\n");

	Matrix<T> matrix(m);
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

	return inverse;
}






