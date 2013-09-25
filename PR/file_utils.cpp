#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "file_utils.h"

using namespace std;

void writeFile(const std::string & fname, const Matrix & m)
{
	ofstream myWriteFile;
	myWriteFile.open(fname, 'w');

	if (!myWriteFile.is_open())
		throw exception("cannot open file for writing");

	// Write the file line by line
	for (Matrix::const_iterator itR = m.begin(), itREnd = m.end(); itR != itREnd; ++itR)
	{
		std::ostream_iterator<unsigned> output_iterator(myWriteFile, "\t");
		std::copy(itR->begin(), itR->end(), output_iterator);
		myWriteFile << "\n";
	}
	
	myWriteFile.close();
}

Matrix readFile(const std::string & fname)
{
	ifstream myReadFile;
	myReadFile.open(fname);

	Matrix bits;

	if (!myReadFile.is_open()) 
		throw exception("cannot open file for reading");

	signed line_size = -1;

	// Read the file line by line
	while (!myReadFile.eof()) 
	{
		string line; //buffer to read into
		while (getline(myReadFile,line))
		{
			std::stringstream strstr(line);
			std::istream_iterator<int> it(strstr), end;
			vector<unsigned> line_bits(it, end);

			// Sanity check
			if (line_size < 0)
				line_size = line_bits.size();
			else if (line_size != line_bits.size())
				throw exception("line size is not stable");

			bits.push_back(line_bits);
		}
 	}

	bits.row = bits.size();
	if (bits.row)
		bits.col = bits.front().size();
			
	myReadFile.close();
	return bits;
}
