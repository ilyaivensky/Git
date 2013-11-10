#ifndef _VECTOR_UTILS_H
#define _VECTOR_UTILS_H

#include <vector>
#include <cstdlib>

#include "matrix.h"

using namespace std;

template <class T> 
T inner_product(const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("inner_product: v1.size() != v2.size()");

	T retval = 0.0;

	for (unsigned i = 0; i < v1.size(); ++i)
		retval += (v1[i] * v2[i]);

	return retval;
}

template <class T> 
Matrix<T> outer_product(const vector<T> & v1, const vector<T> & v2)
{
	return Matrix<T>(v1, v2);
}

template <class T> 
T square_dist(const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("dist: v1.size() != v2.size()");

	T result = 0.0;
	for (vector<T>::const_iterator it1 = v1.begin(), it1End = v1.end(), 
		it2 = v2.begin(); it1 != it1End; ++it1, ++it2)
		result += pow(*it1 - *it2, 2);

	return result;
}
	
template <class T> 
T euclidian_dist(const vector<T> & v1, const vector<T> & v2)
{
	return sqrt(square_dist(v1, v2));
}

template <class T>
T norm(const vector<T> & v, T p = 2.0)
{
	if (p < 1.0)
		throw exception("norm of vector is not defined for p < 1");

	T sum = 0.0;
	for (vector<T>::const_iterator it = v.begin(), itEnd = v.end(); it != itEnd; ++it)
		sum += pow(*it, p);

	return pow(sum, 1/p);
}

template <class T>
ostream & operator<<(ostream & os, const vector<T> & v)
{
	for (unsigned i = 0; i < v.size(); ++i)
		os << v[i] << "\t";

	return os;
}


#endif