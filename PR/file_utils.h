#ifndef _FILE_UTILS_H_
#define _FILE_UTILS_H_

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "LA/matrix.h"

template <class T>
void writeFile(const std::string & fname, const Matrix<T> & m);

template <class T>
Matrix<T> readFile(const std::string & fname);

//////////////////////////////////////////////////////////////
//
// Implementations
//
//////////////////////////////////////////////////////////////

template <class T>
void writeFile(const std::string & fname, const Matrix<T> & m)
{
	std::ofstream myWriteFile;
	myWriteFile.open(fname, 'w');

	if (!myWriteFile.is_open())
		throw exception("cannot open file for writing");

	// Write the file line by line
	for (const auto & row : m)
	{
		std::ostream_iterator<T> output_iterator(myWriteFile, "\t");
		std::copy(row.begin(), row.end(), output_iterator);
		myWriteFile << "\n";
	}

	myWriteFile.close();
}

template <class T>
Matrix<T> readFile(const std::string & fname)
{
	ifstream myReadFile;
	myReadFile.open(fname);
	Matrix<T> data;

	if (!myReadFile.is_open())
	{
		cerr << "No " << fname << endl;
		return data;
	}

	while (!myReadFile.eof())
	{
		string line; //buffer to read into
		while (getline(myReadFile, line))
		{
			std::stringstream strstr(line);

			// use stream iterators to copy the stream to the vector as whitespace separated strings
			std::istream_iterator<T> it(strstr), end;
			std::vector<T> tokens(it, end);

			data.add_row(tokens);
		}
	}

	myReadFile.close();
	return data;
}

#endif