#ifndef _RANDOM_H
#define _RANDOM_H

#include <vector>
#include <conio.h>

using namespace std;

// Returns random in range [-1.0, 1.0]
template<class T>
inline T random(T min = -1.0, T max = 1.0)
{
	T delta = max - min;
	unsigned m = unsigned(100 * delta);
	T r = ((T(rand() % m)) / 100) + min;
	return r;
}


template <class T>
inline vector<T> random_example(unsigned dim)
{
	vector<T> example(dim);
	example[0] = 1.0;
	for (unsigned n = 1; n < dim; ++n)
		example[n] = random<T>();

	return example;
}

template<class T>
inline vector<T> random_example_0(unsigned dim)
{
	vector<T> example(dim);
	for (unsigned n = 0; n < dim; ++n)
		example[n] = random<T>();

	return example;
}

template <class T>
inline vector<T> random_linear_func_3()
{
	vector<T> w(3);
    T x1 = random<T>();
	T y1 = random<T>();
	T x2 = random<T>();
	T y2 = random<T>();

	if (x1 > x2)
	{
		swap(x1, x2);
		swap(y1, y2);
	}

	T dx = x2 - x1;
	T dy = y2 - y1;

	w[1] = dy / dx;
	w[0] = y1 - (w[1] * x1);
	w[2] = -1 * ((w[1] * x1) + w[0]) / y1;

	return w;
}

void random_permutation(vector<unsigned> & v)
{
	for (unsigned i = v.size(); i > 0; --i)
	{
		unsigned r = rand() % i;
		swap(v[r], v[i-1]);
	}
}

#endif