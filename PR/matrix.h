#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <vector>
#include <string>

using namespace std;

// Some infrastructure workout
struct Matrix : public vector<vector<unsigned> > 
{
	typedef vector<unsigned> Row;
			
	// Creates matrix row * col and initializes with 0
	Matrix(unsigned row, unsigned col) : 
		vector<Row>(row, Row(col, 0)), row(row), col(col) {}

	// Creates square matrix n * n and initializes with 0
 	Matrix(unsigned n) : 
		vector<Row>(n, Row(n, 0)), row(n), col(n) {}

	// Creates empty matrix
	Matrix() : row(0), col(0) {}

	unsigned row;
	unsigned col;
};

ostream & operator<<(ostream & os, const vector<unsigned> & v);
ostream & operator<<(ostream & os, const Matrix & m);

#endif