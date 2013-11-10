/*                                                                 -*- C++ -*-
 * File: knn.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 10, 2013
 *
 * Description:
 *   Implementation of KNN classifier
 *   
 */

#ifndef _KNN_H
#define _KNN_H

#include "LA/Matrix.h"
#include <ostream>
#include <map>

class KNN 
{
public:
	KNN(unsigned k) : neighbours(k) {}

	template <class DataT> 
	Matrix<unsigned> classify(const Matrix<DataT> & trainingData, const Matrix<unsigned> & trainingLabels,
							  const Matrix<DataT> & testingData, const Matrix<unsigned> & testingLabels,
							  unsigned numClasses, ostream & os) const;
private:
	unsigned neighbours;
};

template <class DataT> 
Matrix<unsigned> KNN::classify(const Matrix<DataT> & trainingData, const Matrix<unsigned> & trainingLabels,
							  const Matrix<DataT> & testingData, const Matrix<unsigned> & testingLabels,
							  unsigned numClasses, ostream & os) const
{
	Matrix<unsigned> confusion(numClasses, numClasses);

	// Foe each example ...
	for (unsigned test = 0; test < testingData.row; ++test)
	{
		// Step 1: find nearest members and their labels

		vector<DataT> min_dist(this->neighbours, std::numeric_limits<DataT>::max());
		vector<unsigned> nn(this->neighbours);

		for (unsigned train = 0; train < trainingData.row; ++train)
		{
			// Actually, we have to calculate Euclidian distance
			// Instead, for efficiency reason, we calculate only a sum of quadrats, and skip a sqrt 
			float dist = square_dist(trainingData[train], testingData[test]);

			// Find nearest neighbours
			for (unsigned i = 0; i < nn.size(); ++i)
			{
				if (dist < min_dist[i])
				{
					if (i < nn.size() - 1)
					{
						min_dist[i+1] = min_dist[i];
						nn[i+1] = nn[0];
					}
					min_dist[i] = dist;
					nn[i] = trainingLabels[train][0];
				}
			}
		}

		// Step 2: Find label by looking at nearest neighbours

		// Count all predicted labels
		map<unsigned, unsigned> predicted;
		for (unsigned i = 0; i <  this->neighbours; ++i)
			predicted[nn[i]] += 1;

		// Find the most popular prediction
		unsigned maxPredictions = 0;
		unsigned label = nn[0]; // initialize by the nearest neighbour
		for (map<unsigned, unsigned>::const_iterator it = predicted.begin(), 
			itEnd = predicted.end(); it != itEnd; ++it)
		{
			if (it->second > maxPredictions)
				maxPredictions = it->second; label = it->first;
		}

		// Step 3: Update confusion matrix
		confusion[label][testingLabels[test][0]] += 1;
		
		// Report
		os << "ID: " << test << " Pred: " << label << ", corr: " << testingLabels[test][0] 
			<< " NNs (sorted): " << nn << endl; 
	}

	return confusion;
}

#endif