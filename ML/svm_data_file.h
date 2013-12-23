/*                                                                 -*- C++ -*-
 * File: svm_data_file.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 26, 2013
 *
 * Description:
 *   Utility to dump features into SVM supported data file
 *   
 */

#ifndef _SVM_DATA_FILE_H_
#define _SVM_DATA_FILE_H_

#include <fstream>
#include <iostream>

#include "LA/matrix.h"

template <class D, class L>
void write_SVM_data_file(const std::string & fileName, Matrix<D> & features, vector<L> & labels)
{
	std::cerr << "Writing " << fileName << std::endl;
	// Open output file
	ofstream of;
	stringstream str;
	str << fileName;
	of.open(str.str(), 'w');
	
	unsigned i = 0;
	auto itL = labels.begin();
	for (auto itD = features.begin(), itDEnd = features.end(); 
		itD != itDEnd; ++itD, ++itL, ++i)
	{
		of << *itL << " ";
		unsigned id = 1;
		for (auto itF = itD->begin(), itFEnd = itD->end(); itF != itFEnd; ++itF, ++id)
				of << id << ":" << *itF << " ";
		of << "\n";
		if (i % 1000 == 0)
			cerr << i << endl;
	}

	of.close();
	cerr << "Done!" << endl;
}

#endif