#ifndef _PERCEPTRON_H
#define _PERCEPTRON_H

#include "matrix.h"

template <class T>
unsigned train_perceptron(
	const Matrix<T> & data, 
	const vector<signed> & desired_target,
	vector<T> & w)
{
	unsigned count = 0;

	//cerr << "Perceptron input w=[" << w << "]" << endl;
	//cerr << "Data:" << endl << data << endl;
	//cerr << "Desired target: " << desired_target << endl;

	while (true)
	{
		// Index of misclassified data + error sign (+ or -)
		vector<pair<unsigned, signed> > misclassified; 

		// Test and get a list of misclassified
		for (unsigned m = 0; m < 10/*data.row*/; ++m)
		{
			float res = 0.0;
			for (unsigned n = 0; n < data.col; ++n)
				res += (w[n] * data[m][n]);

			signed g = (res > 0 ? 1 : 0);
			signed error = desired_target[m] - g;
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
			vector<float> next_w(w.size());
			for (unsigned n = 0; n < data.col; ++n)
				next_w[n] = w[n] + error * data[m][n];

			if (w == next_w)
				throw exception("Failed to ajust weights\n");

			swap(w, next_w);
		}
	}
	
	return count;
}

#endif