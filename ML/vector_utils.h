#ifndef _VECTOR_UTILS_H
#define _VECTOR_UTILS_H

#include <vector>
#include <cstdlib>

using namespace std;

template <class T> 
T dot_product(const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("dot_product: v1.size() != v2.size()");

	T retval = 0.0;

	for (unsigned i = 0; i < v1.size(); ++i)
		retval += (v1[i] * v2[i]);

	return retval;
}

template <class T> 
T dist(const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("dist: v1.size() != v2.size()");

	T result = 0.0;
	for (unsigned i = 0; i < v1.size(); ++i)
		result += pow(v1[i] - v2[i], 2);

	return sqrt(result);
}

template <class T>
ostream & operator<<(ostream & os, const vector<T> & v)
{
	for (unsigned i = 0; i < v.size(); ++i)
		cerr << v[i] << "\t";

	return os;
}


#endif