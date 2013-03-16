/*                                                                 -*- C++ -*-
 * File: rbf.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Mar 16, 2013
 *
 * Description:
 *   RBF classifier 
 *   
 */

#ifndef _RBF_H
#define _RBF_H

#include "linear_solutions.h"
#include "clustering.h"
#include "matrix.h"

using namespace std;

template <class T>
class RBF_Classifier
{
public:

	// Creates and configures classifier
	RBF_Classifier(T gamma, unsigned k) : gamma(gamma), k(k){}
	
	// Trains using dataset x with labels y. 
	void train(const Matrix<T> & x, const Matrix<T> & y);
	
	// Predict single example
	T predict(const vector<T> & x) const;

	// Predict data set. Returns matrix 1*n
	Matrix<T> predict(const Matrix<T> & x) const;

private:
	
	// Parameters
	T gamma; 
	unsigned k;
    
	// Trained 
	Matrix<T> w;
	Matrix<T> mu;
};

template <class T>
void RBF_Classifier<T>::train(const Matrix<T> & x, const Matrix<T> & y)
{
	//Matrix<T> mu;
	while (mu.empty())
		mu = k_means_clusterung(x, k);

	Matrix<T> phi(x.row, k + 1);
	for (unsigned n = 0; n < phi.row; ++n)
	{
		phi[n][0] = 1.0; // adding for calculating bias (aka w[0])
		for (unsigned k = 1; k < phi.col; ++k)
			phi[n][k] = exp((-1) * gamma * square_dist(x[n], mu[k - 1])); 
	}

	w = linear_pseudoinverse_solution(phi, y);
}

template <class T>
T RBF_Classifier<T>::predict(const vector<T> & x) const
{
	T retval = w[0][0]; // bias
	for (unsigned r = 0; r < mu.row; ++r)
		retval += w[r + 1][0] * exp((-1) * gamma * square_dist(x, mu[r])); 

	return retval /= k;
}

template <class T>
Matrix<T> RBF_Classifier<T>::predict(const Matrix<T> & x) const
{
	Matrix<T> h;

	// predict each example
	for (unsigned r = 0; r < x.row; ++r)
		h.add_row(vector<T>(1, predict(x[r])));
	
	return h;
}


#endif
