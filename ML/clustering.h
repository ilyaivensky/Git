/*                                                                 -*- C++ -*-
 * File: clustering.h
 * 
 * Author: Ilya Ivensky
 * Created on: Mar 12, 2013
 *
 * Description: implements clustering algorithms   
 */

#ifndef _CLUSTERING_H
#define _CLUSTERING_H

#include "LA/matrix.h"

using namespace std;

/************************************************************
* Implements Lloyd's algorithm
* Returns matrix of k cluster centers
* (or empty matrix if fails to create k non-empty clusters)
************************************************************/
template <class T>
Matrix<T> k_means_clusterung(const Matrix<T> & x, unsigned k)
{
	// Select initial k points
	Matrix<T> mu(k, x.col);
	mu.random_init_0();

	unsigned count = 0;

	while (true)
	{
		Matrix<T> mu_next(mu.row, mu.col);
		// Each cluster contains IDs of examples
		typedef vector<unsigned> Cluster;
		vector<Cluster> clusters(k);

		/******************************************
		* Create clusters for fixed centers
		*******************************************/
		for (unsigned n = 0; n < x.row; ++n)
		{
			unsigned nearest = mu.size();
			T min_dist = numeric_limits<T>::max();
			for (unsigned c = 0; c < mu.size(); ++c)
			{
				T curr_dist = euclidian_dist(x[n], mu[c]);
				if (curr_dist < min_dist)
				{
					min_dist = curr_dist;
					nearest = c;
				}
			}
			clusters[nearest].push_back(n);
		}

		/*******************************************
		* Update centers for fixed clusters 
		********************************************/
		for (unsigned c = 0; c < clusters.size(); ++c)
		{
			const Cluster & cluster = clusters[c];

			if (cluster.empty())
			{
				//cerr << "k_means_clustering: cluster is empty" << endl;
				return Matrix<T>();
			}

			vector<T> & mu_c = mu_next[c]; 
			
			for (vector<unsigned>::const_iterator itC = cluster.begin(), itCEnd = cluster.end(); 
				itC != itCEnd; ++itC)
			{
				const vector<T> & example = x[*itC];
				vector<T>::iterator itMu = mu_c.begin(); 
				for (vector<T>::const_iterator itEx = example.begin(), 
					itExEnd = example.end(); itEx != itExEnd; ++itEx, ++itMu)
					*itMu += *itEx;
			}

			// Average each feature
			for (vector<T>::iterator itMu = mu_c.begin(), itMuEnd = mu_c.end(); itMu != itMuEnd; ++itMu)
				*itMu /= cluster.size();
		}

		swap(mu, mu_next);

		if (mu == mu_next)
			break;

		++count;
	}

	//cerr << "Lloyd's done, count=" << count << endl;
	return mu;
}

#endif