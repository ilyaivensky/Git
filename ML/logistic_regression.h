#ifndef _LOGISTIC_REGRESSION_H
#define _LOGISTIC_REGRESSION_H

#include <cstdlib>
#include <algorithm>
#include "matrix.h"
#include "vector_utils.h"

extern double e;
extern bool debug;

template <class T> 
unsigned logistic_regression(const Matrix<T> & x, const Matrix<T> & y, vector<T> & w)
{
	const T learning_rate = 0.01;

	unsigned epoch = 1;

	while (true)
	{
		// Implements SGD
		vector<unsigned> examples(x.row);
		for (unsigned i = 0; i < x.row; ++i)
			examples[i] = i;

		std::random_shuffle(examples.begin(), examples.end());

		vector<T> prev_w(w);
		for (unsigned i = 0; i < examples.size(); ++i)
		{
			unsigned m = examples[i];
			vector<T> next_w(w);
		
			T sig = dot_product(x[m], w);
			T agreement = sig * y[m][0];
			T k = (1 + pow(e, agreement));
			if (debug)
				cerr << "example " << m << " agreement " << agreement 
				<< " y[" << m << "]=" << y[m][0] << " dot_product=" << sig
				<< " k=" << k << endl;
		
			if (k == 0.0)
				throw exception("div is 0");
			for (unsigned n = 0; n < x.col; ++n)
				next_w[n] += learning_rate * ((x[m][n] * y[m][0]) /  k);

			swap(w, next_w);
		}

		T d = dist(w, prev_w);

		if (debug)
			cerr << "Iter " << epoch << " distance " << d << endl;

		if (d < 0.01)
		{
			if (debug)
			{
				cerr << "w(" << epoch-1 << ") " << prev_w << endl;
				cerr << "w(" << epoch << ") " << w << endl;
			}
			break;
		}

		if (epoch % 100 == 0)
			cerr << "epoch " << epoch << endl;

		++epoch;
	}

	return epoch;
}

#endif