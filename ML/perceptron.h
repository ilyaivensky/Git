#ifndef _PERCEPTRON_H
#define _PERCEPTRON_H

#include "matrix.h"

namespace PLA
{
template <class T>
signed label(const T & r)
{
	return r > 0 ? 1 : 0;
}

template <class T>
unsigned train_perceptron(
	const Matrix<T> & x, 
	const Matrix<T> & y,
	vector<T> & w)
{
	unsigned count = 0;

	vector<signed> target(y.row);
	for (unsigned r = 0; r < y.row; ++r)
		target[r] = label(y[r][0]); 

	//cerr << "Perceptron input w=[" << w << "]" << endl;
	//cerr << "Data:" << endl << data << endl;
	//cerr << "Desired target: " << desired_target << endl;

	while (true)
	{
		// Index of misclassified data + error sign (+ or -)
		vector<pair<unsigned, signed> > misclassified; 

		// Test and get a list of misclassified
		for (unsigned m = 0; m < x.row; ++m)
		{
			signed g = label(dot_product(w, x[m]));
			signed error = target[m] - g; // set error -1 or 0 or 1
			if (error)
			{
				misclassified.push_back(make_pair(m, error)); 
			}
		}

		++count;

		if (misclassified.empty())
		{
			//cerr << "DONE! Precision 1.0" << endl;
		    break;
		}
		else 
		{
			// Select random misclassified 
			int toFix = rand() % misclassified.size();
			signed error = misclassified[toFix].second;
			unsigned m = misclassified[toFix].first;
			
			// Update w using selected random misclassified
			vector<T> next_w(w.size());
			for (unsigned n = 0; n < x.col; ++n)
				next_w[n] = w[n] + error * x[m][n];

			if (w == next_w)
				throw exception("Failed to ajust weights\n");

			swap(w, next_w);
		}
	}
	
	return count;
}

} // namespace
#endif