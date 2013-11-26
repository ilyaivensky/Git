#ifndef _SVM_FEATURES_FILE_H_
#define _SVM_FEATURES_FILE_H_

#include <fstream>
#include <iostream>

#include "LA/matrix.h"

template <class D, class L>
void write_SVM_feature_file(const std::string & fileName, Matrix<D> & features, vector<L> & labels)
{
	std::cerr << "Writing " << fileName << std::endl;
	// Open output file
	ofstream of;
	stringstream str;
	str << fileName;
	of.open(str.str(), 'w');
	
	unsigned i = 0;
	vector<L>::const_iterator itL = labels.begin();
	for (Matrix<D>::const_iterator itD = features.begin(), itDEnd = features.end(); 
		itD != itDEnd; ++itD, ++itL, ++i)
	{
		of << *itL << " ";
		unsigned id = 1;
		for (vector<D>::const_iterator itF = itD->begin(), itFEnd = itD->end(); itF != itFEnd; ++itF, ++id)
				of << id << ":" << *itF << " ";
		of << "\n";
		if (i % 1000 == 0)
			cerr << i << endl;
	}

	of.close();
	cerr << "Done!" << endl;
}

#endif