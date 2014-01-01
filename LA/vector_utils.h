#ifndef _VECTOR_UTILS_H_
#define _VECTOR_UTILS_H_

#include <vector>

using namespace std;

#define EPSILON 1e-5

template <class T>
bool is_zero(const vector<T> & v);

template <class T> 
T inner_product(const vector<T> & v1, const vector<T> & v2);

template <class T> 
Matrix<T> outer_product(const vector<T> & v1, const vector<T> & v2);

template <class T> 
T square_dist(const vector<T> & v1, const vector<T> & v2);
	
template <class T> 
T euclidian_dist(const vector<T> & v1, const vector<T> & v2);

double norm(const vector<int> & v, unsigned p);

template <class T>
T norm(const vector<T> & v, unsigned p);

template <class T>
T cosine(const vector<T> & v1, const vector<T> & v2);

template <class T>
void make_vector_set(vector<T> & v);

template <class T>
T operator * (const vector<T> & v1, const vector<T> & v2);

template <class T>
vector<T> & operator *= (vector<T> & v, const T & scalar);

template <class T>
vector<T> operator * (const vector<T> & v, const T & scalar);

template <class T>
vector<T> & operator /= (vector<T> & v, const T & scalar);

template <class T>
vector<T> operator / (const vector<T> & v, const T & scalar);

template <class T>
vector<T> & operator += (vector<T> & v1, const vector<T> & v2);

template <class T>
vector<T> & operator -= (vector<T> & v1, const vector<T> & v2);

template <class T>
ostream & operator << (ostream & os, const vector<T> & v);

#include "LA/vector_utils_impl.h"

#endif
