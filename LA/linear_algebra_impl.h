/*                                                                 -*- C++ -*-
 * File: linear_algebra_impl.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 7, 2013
 *   
 */

#include <limits>
#include <cassert>

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

			for (const auto & row : m)
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
	auto deviation = dev(m);
	auto cov = gram(deviation);
	cov /= static_cast<T>(m.nrow());
	return cov;
}

template <class T>
Matrix<T> inv(const Matrix<T> & m)
{
	if (!is_square(m))
		throw exception("cannot inverse non-square matrix\n");

	unsigned n = m.size();
	Matrix<T> upper = m, inverse = Matrix<T>::diag(n, 1);
	
	for (unsigned pivot = 0; pivot < n; ++pivot)
	{
		for (unsigned r = 0; r < n; ++r)
		{
			if (r != pivot)
			{
				if (upper[pivot][pivot] == 0)
					throw exception("matrix cannot be inverted\n");
				
				T ratio = upper[r][pivot]/upper[pivot][pivot];

				upper[r] -= upper[pivot] * ratio;
				inverse[r] -= inverse[pivot] * ratio;
			}
		}
	}

	for (unsigned i = 0; i < n; ++i)
	{
		T a = upper[i][i];
		if (a == 0)
			throw exception("matrix cannot be inverted\n");

		inverse[i] /= a;
	}

	return inverse;
}

template <class T>
vector<T> characteristic_polynomial(const Matrix<T> & m)
{
	if (!is_square(m))
		throw exception("characteristic_polynomial(): Matrix is not square");

	switch (m.size())
	{
	case 2:
		return vector<T>({ 1, (-1) * tr(m), det(m) });
	case 3:
	{
		T c = (pow(tr(m), 2) - tr(m * m)) / (-2);
		return vector<T>({ (-1), tr(m), c, det(m) });
	}
	default:
		break;
	}

	throw exception("characteristic_polynomial(): matrices of this size are not supported");
	return vector<T>();
}

template <class T>
vector<T> solve_quadratic(const vector<T> & p)
{
	if (p.size() != 3)
		throw exception("solve_quadratic(): polynom length is not 3");

	T d = sqrt(pow(p[1], 2) - (4 * p[0] * p[2]));
	T a2 = 2 * p[0];
	T b = (-1) * p[1];

	if (d == 0)
		return vector<T>(1, (b / a2));

	T root1 = (b - d) / a2;
	T root2 = (b + d) / a2;

	return vector<T>( { root1, root2 } );
}

template <class T>
vector<T> eigenvalues_2x2(const Matrix<T> & m)
{
	if (m.nrow() != 2 || m.ncol() != 2)
		throw exception("eigen_values_2x2(): matrices of this size are not supported");

	auto cp = characteristic_polynomial(m);

	return solve_quadratic(cp);
}

template <class T>
tuple <Matrix<T>, Matrix<T>, Matrix<T>> lu(const Matrix<T> & m)
{
	if (!is_square(m))
		throw exception("lu: m is not square");

	Matrix<T> lower(Matrix<T>::diag(m.nrow(), 1));
	Matrix<T> upper(m);
	Matrix<T> permutation(Matrix<T>::diag(m.nrow(), 1));

	unsigned n = m.nrow();

	// Find the optimal permutation matrix (P)
	for (unsigned pivot = 0; pivot < n; ++pivot)
	{
		T max_pivot = 0;
		unsigned max_pivot_row = m.nrow();
		for (unsigned r = pivot; r < n; ++r)
		{
			if (abs(upper[r][pivot]) > abs(max_pivot))
			{
				max_pivot_row = r;
				max_pivot = upper[r][pivot];
			}
		}

		if (max_pivot == 0)
			throw exception("lu: cannot decompose m because it is singular");

		if (max_pivot_row != pivot)
		{
			swap(upper[max_pivot_row], upper[pivot]);
			swap(permutation[max_pivot_row], permutation[pivot]);
		}
	}

	// Generate lower (L) and upper (U) matrices
	for (unsigned pivot = 0; pivot < n; ++pivot)
	{
		for (unsigned r = pivot + 1; r < n; ++r)
		{
			T ratio = upper[r][pivot] / upper[pivot][pivot];
			upper[r] -= upper[pivot] * ratio;
			lower[r][pivot] = ratio;
		}
	}

	return make_tuple(lower, upper, permutation);
}

template <class T>
Matrix<T> rref_int(const Matrix<T> & m)
{
	Matrix<T> a(m);
	for (unsigned pivot_row = 0, pivot_col = 0; pivot_col < a.ncol(); ++pivot_col)
	{
		// Find the best pivot row for this column
		T pivot = 0;
		unsigned max_pivot_row = m.nrow();
		for (unsigned r = pivot_row; r < a.nrow(); ++r)
		{
			auto & row = a[r];
			if (row[pivot_col] != 0)
			{
				pivot = row[pivot_col];
				max_pivot_row = r;
				break;
			}
		}

		if (pivot == 0)
			continue;

		auto k = a[max_pivot_row][pivot_col];
		if (k != 1) a[max_pivot_row] /= k;

		if (max_pivot_row != pivot_row)
			swap(a[max_pivot_row], a[pivot_row]);

		auto & prow = a[pivot_row];
		// Substruct pivot row (multiplied by proper k) from each non-pivot row
		for (unsigned r = 0; r < a.nrow(); ++r)
		{
			if (r == pivot_row)
				continue;

			auto & row = a[r];
			auto k = row[pivot_col];
			for (unsigned c = pivot_col; c < a.ncol(); ++c)
			{
				if (row[c] == row[pivot_col]) row[c] = 0;
				else row[c] -= prow[c] * k;
			}
		}

		++pivot_row;
	}

	// Eliminate zero rows
	// If there are any, they should be at the bottom of matrix
	while (!a.empty())
	{
		const auto & row = a.back();
		if (is_zero(row))
			a.pop_back();
		else
			break;
	}

	return a;
}

Matrix<int> rref(const Matrix<int> & m)
{
	return rref_int(m);
}

Matrix<long> rref(const Matrix<long> & m)
{
	return rref_int(m);
}

template <class T>
Matrix<T> rref(const Matrix<T> & m)
{
	Matrix<T> a(m);
	for (unsigned pivot_row = 0, pivot_col = 0; pivot_col < a.ncol(); ++pivot_col)
	{
		// Find the best pivot row for this column
		T max_pivot = 0;
		unsigned max_pivot_row = m.nrow();
		for (unsigned r = pivot_row; r < a.nrow(); ++r)
		{
			auto & row = a[r];
			if (abs(row[pivot_col]) > abs(max_pivot))
			{
				max_pivot = row[pivot_col];
				max_pivot_row = r;
			}
		}

		if (max_pivot == 0)
			continue;

		auto k = a[max_pivot_row][pivot_col];
		if (k != 1) a[max_pivot_row] /= k;

		if (max_pivot_row != pivot_row)
			swap(a[max_pivot_row], a[pivot_row]);

		auto & prow = a[pivot_row];
		// Substruct pivot row (multiplied by proper k) from each non-pivot row
		for (unsigned r = 0; r < a.nrow(); ++r)
		{
			if (r == pivot_row)
				continue;

			auto & row = a[r];
			auto k = row[pivot_col];
			for (unsigned c = pivot_col; c < a.ncol(); ++c)
			{
				if (row[c] == row[pivot_col]) row[c] = 0;
				else
				{
					row[c] -= prow[c] * k;
					auto rounded = round(row[c]);
					if (fabs(row[c] - rounded) < EPSILON) row[c] = rounded;
				}
			}
		}

		++pivot_row;
	}

	// Eliminate zero rows
	// If there are any, they should be at the bottom of matrix
	while (!a.empty())
	{
		const auto & row = a.back();
		if (is_zero(row))
			a.pop_back();
		else
			break;
	}

	return a;
}

template <class T>
vector<pair<T, vector<T>>> eigen_2x2(const Matrix<T> & m)
{
	auto eigenvalues = eigenvalues_2x2(m);

	static Matrix<T> zeros(m.nrow(), 1);

	vector<std::pair<T, vector<T>>> eigens;

	for (const auto & eigenval : eigenvalues)
	{
		auto a = m - Matrix<T>::diag(m.nrow(), eigenval);
		auto constraints = rref(a);

		// Because det(a) has to be 0
		assert(constraints.nrow() == 1);
		
		/*
		   constraints[0][0] has to be always 1 
		   (because it is echelon form of a)
		   so if we set eigenvec[0] to be 1,
		   then eigenvec[1] = -1 / constraints[0][1]
		   because a * eigenvec = 0
		 */
		vector<T> eigenvec = { 1, (-1) / constraints[0][1] };
		eigens.push_back(make_pair(eigenval, eigenvec));
	}
	
	return eigens;
}




