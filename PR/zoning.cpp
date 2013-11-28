/*                                                                 -*- C++ -*-
 * File: zoning.cpp
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 27, 2013
 *
 * Description:
 *   Utility to break image into zones
 *   
 */


#include "PR/zoning.h"

#include <iostream>

Matrix<Zone> zoning(unsigned row, unsigned col, unsigned num)
{
	if (num > row || num > col || num == 0)
		throw exception("Illegal zoning parameter");

	Matrix<Zone> zones(num, num);

	double zone_row = static_cast<double>(row) / num;
	double zone_col = static_cast<double>(col) / num;

	// Initialize the right most zones (by vertical) 
	// end column to be the end column of image  
	for (unsigned i = 0; i < num; ++i)
		zones[i][num - 1].colsEnd() = col;

	// Initialize the bottom zones 
	// row to be the end row of image
	for (unsigned j = 0; j < num; ++j)
		zones[num - 1][j].rowsEnd() = row;

	// Initialize start row and start col for all zones
	for (unsigned i = 0; i < num; ++i)
		for (unsigned j = 0; j < num; ++j)
		{
			Zone & zone = zones[i][j];

			// Add 0.5 to double to force rounding to nearest integer
			zone.rowsBegin() = static_cast<unsigned>(i * zone_row + 0.5);
			zone.colsBegin() = static_cast<unsigned>(j * zone_col + 0.5); 
		}

	// Initialize end rows for all zones except the most bottom zones 
	// (bottom zones were initialized earlier)
	for (unsigned i = num - 1; i > 0; --i)
		for (unsigned j = 0; j < num; ++j)
			zones[i - 1][j].rowsEnd() = zones[i][j].rowsBegin();

	// Initialize end columns for all zones except the right most zones 
	// (right most zones were initialized earlier)
	for (unsigned i = 0; i < num; ++i)
		for (unsigned j = num - 1; j > 0; --j)
			zones[i][j - 1].colsEnd() = zones[i][j].colsBegin();

	//cerr << "Zoning: " << zones << endl;
	return zones;
}

ostream & operator<<(ostream & os, const Zone & z)
{
	os << "((" << z.first.first << "," << z.first.second << "),(" << z.second.first << "," << z.second.second << "))";
	return os;
}