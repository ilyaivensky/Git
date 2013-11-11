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
	// Lazy training - just memorize the data and labels
	void train(const Matrix<T> & trainingData, const Matrix<unsigned> & trainingLabels)
	{
		trainedData_ = trainingData;
		trainedLabels_ = trainingLabels;
	}

	Matrix<unsigned> classify(const Matrix<T> & data, 
		unsigned numNN, 
		ostream & os, 
		const Matrix<unsigned> & refLabels = Matrix<unsigned>()) const;

private:
	Matrix<T> trainedData_;
	Matrix<unsigned> trainedLabels_;
};

template <class T> 
Matrix<unsigned> KNN<T>::classify(const Matrix<T> & data, 
								  unsigned numNN,
								  ostream & os, 
								  const Matrix<unsigned> & refLabels) const
{
	// Column vector of predictions
	Matrix<unsigned> predictedLabels(data.row, 1);

	// For each example ...
	for (unsigned m = 0; m < data.row; ++m)
	{
		// Step 1: find nearest members and their labels

		vector<T> min_dist(numNN, std::numeric_limits<T>::max());
		vector<unsigned> nn(numNN);

		for (unsigned tr = 0; tr < trainedData_.row; ++tr)
		{
			// Actually, we have to calculate Euclidian distance
			// Instead, for efficiency reason, we calculate only a sum of quadrats and skip a sqrt
			// This optimization is transparent for results because sqrt is monotone ascending function 
			float dist = square_dist(trainedData_[tr], data[m]);

			// Find nearest neighbours using priority queue
			for (unsigned i = 0; i < nn.size(); ++i)
			{
				if (dist >= min_dist[i])
					continue;

				// Current element is closer than any other element in the queue
				// Shift all of them 1 position up 
				for (unsigned j = i + 1; j < nn.size(); ++j)
				{
					min_dist[j] = min_dist[j-1];
					nn[j] = nn[j-1];
				}
					
				// ...and insert current element
				min_dist[i] = dist;
				nn[i] = trainedLabels_[tr][0];
			}
		}

		// Step 2: Find label by looking at nearest neighbours

		// Count all predicted labels
		map<unsigned, unsigned> predicted;
		for (unsigned i = 0; i < numNN; ++i)
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

		predictedLabels[m][0] = label;

		if (!refLabels.empty())
			os << "ID: " << m << " Pred: " << label << " Corr: " << refLabels[m] << " NNs: " << nn << endl;
		else
			os << "ID: " << m << " Pred: " << label << " NNs: " << nn << endl;
	}
	
	return predictedLabels;
}

#endif