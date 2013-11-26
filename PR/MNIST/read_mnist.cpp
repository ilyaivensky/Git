/*                                                                 -*- C++ -*-
 * File: read_mnist.cpp
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 26, 2013
 *
 * Description:
 *   Utility file to extract data from a pair (image, label) raw MNIST data files
 * 
 * NOTE: this implementation works only on little endian machines
 * (change function read32() to use this code on big endian machines 
 *   
 */

#include "PR/MNIST/read_mnist.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace MNIST {

static uint32_t read32(fstream & fs)
{
	uint8_t buff[4];
	fs.read(reinterpret_cast<char*>(&buff), sizeof(buff));
	uint32_t retval = 0;
	retval = buff[3];  
	retval |= buff[2] << 8;
	retval |= buff[1] << 16; 
	retval |= buff[0] << 24;

	return retval;
}

static inline uint8_t read8(fstream & fs)
{
	uint8_t val;
	fs.read(reinterpret_cast<char*>(&val), sizeof(val));
	return val;
}

int read_mnist(const std::string & imagesFile, const std::string & labelsFile, 
			 vector<Image> & images, vector<int> & labels, 
			 const std::string & outputFile, unsigned maxImgs)
{
	fstream ifs, lfs;
	ifs.open(imagesFile, std::ios::in | std::ios::binary);
	lfs.open(labelsFile, std::ios::in | std::ios::binary);
	
	if (!ifs.is_open())
	{
		cerr << "No " << imagesFile << endl;
		throw exception("file cannot be open");
	}

	if (!lfs.is_open())
	{
		cerr << "No " << labelsFile << endl;
		throw exception("file cannot be open");
	}
	
	// Read headers
	uint32_t magicFS = read32(ifs);
	uint32_t numImgs = read32(ifs);
	uint32_t numRow = read32(ifs);
	uint32_t numCol = read32(ifs);

	uint32_t magicLS = read32(lfs);
	uint32_t numLabels = read32(lfs);

	if (numImgs != numLabels)
		throw exception("numImgs != numLabels");

#if 1
	cerr << "image file, magic = "
		<< (int) magicFS << ", numImgs = " << (int) numImgs 
		<< ", numRow = " << (int) numRow << ", numCol = " << numCol << endl;
	
	cerr << "labels file, magic = "
		<< (int) magicLS << ", numLabels = " << (int) numLabels << endl;
#endif
/*
	// Open output files
	vector<ofstream> ofs(10);
	for (unsigned i = 0; i < ofs.size(); ++i)
	{
		stringstream str;
		str << outputFile << i << ".txt";
		ofs[i].open(str.str(), 'w');
	}
	*/
	vector<unsigned> counters(10, 0);

	images.resize(std::min(numImgs, maxImgs), Image(numRow, numCol));
	labels.resize(std::min(numImgs, maxImgs), -1);

	unsigned i = 1;
	vector<int>::iterator itLbls = labels.begin();
	for (vector<Image>::iterator itImgs = images.begin(), itEnd = images.end(); itImgs != itEnd; ++itImgs, ++itLbls, ++i)
	{
		Image & img = *itImgs;
		for (unsigned r = 0; r < numRow; ++r)
		{
			for (unsigned c = 0; c < numCol; ++c)
			{
				uint8_t val = read8(ifs);
				img[r][c] = val ? 1 : 0;
			}
		}

		int label = read8(lfs);
		*itLbls = label;

		//ofs[label] << img << endl;
		counters[label] += 1;
		//cerr << label << endl;

		if ((i % 1000) == 0)
			cerr << i << endl;
	}
/*
	for (unsigned i = 0; i < ofs.size(); ++i)
	{
		ofs[i].close();
		cerr << "Total " << i << ": " << counters[i] << endl;
	}
	*/
	return 0;
}

} // namespace