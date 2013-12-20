/*                                                                 -*- C++ -*-
 * File: linear_algebra_impl.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 7, 2013
 *   
 */

template <class T>
T tr(const Matrix<T> & m)
{
	if (!m.is_square())
		throw exception("trace is not defined for non-square matrix");

	T tr = 0;
	for (unsigned i = 0; i < m.nrow(); ++i)
		tr += m[i][i];

	return tr;
}

template <class T>
T det(const Matrix<T> & m)
{
	if (!m.is_square())
		throw exception("determinant is not defined for non-square matrix"); 

	T det = 0;
	for (unsigned i = 0; i < m.nrow(); ++i)
	{
		for (unsigned r = 0, c = i; r < m.nrow(); ++r, ++c)
		{
			if (c == m.nrow()) c = 0;
			det += m[r][c];
		}

		for (unsigned r = 0, c = m.nrow() - i; r < m.nrow(); ++r, --c)
		{
			if (c == 0) c = m.nrow();
			det -= m[r][c - 1];
		}
	}
	return det;
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
			for (Matrix<T>::const_iterator it = m.begin(), itEnd = m.end(); it != itEnd; ++it)
				rij += (*it)[i] * (*it)[j];

			res[j][i] = rij;
		}
	}

	return res;
}

template <class T>
vector<T> medc(const Matrix<T> & m)
{
	Matrix<T> k(1, m.nrow(), 1);
	Matrix<T> res(k * m);
	res /= static_cast<T>(m.nrow());
	return res[0];
}

template <class T>
Matrix<T> dev(const Matrix<T> & m)
{
	return m - Matrix<T>(medc(m), m.nrow());
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
	if (!m.is_square())
		throw exception("cannot inverse non-square matrix\n");

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






