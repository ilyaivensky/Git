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
#include <set>

template <class T>
class KNN 
{
public:
	// Lazy training - just memorize
	void train(const Matrix<T> & trainingData, const Matrix<unsigned> & trainingLabels)
	{
		trainedData_ = trainingData;
		trainedLabels_ = trainingLabels;
	}

	Matrix<unsigned> classify(const Matrix<T> & data, unsigned numNN, 
		ostream & os, 
		const Matrix<unsigned> & refLabels = Matrix<unsigned>()) const;

private:
	Matrix<T> trainedData_;
	Matrix<unsigned> trainedLabels_;
};

template <class T> 
Matrix<unsigned> KNN<T>::classify(const Matrix<T> & data, unsigned numNN,
								  ostream & os, 
								  const Matrix<unsigned> & refLabels) const
{
	Matrix<unsigned> predictedLabels(data.row, 1);

	// For each example ...
	for (unsigned m = 0; m < data.row; ++m)
	{
		// Step 1: find nearest members and their labels

		vector<T> min_dist(numNN, std::numeric_limits<T>::max());
		vector<unsigned> nn(numNN);

		for (unsigned train = 0; train < trainedData_.row; ++train)
		{
			// Actually, we have to calculate Euclidian distance
			// Instead, for efficiency reason, we calculate only a sum of quadrats, and skip a sqrt 
			float dist = square_dist(trainedData_[train], data[m]);

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