/*                                                                 -*- C++ -*-
 * File: zoning.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 27, 2013
 *
 * Description:
 *   Utility to break image into zones
 *   
 */

#ifndef _ZONING_H_
#define _ZONING_H_

#include "LA/matrix.h"
#include <iostream>

using namespace std;

// Defines zone on the image (2D range)
struct Zone : public pair<pair<unsigned, unsigned>, pair<unsigned, unsigned> > 
{
	typedef public pair<pair<unsigned, unsigned>, pair<unsigned, unsigned> > Super;
		
	Zone() {}
	Zone(unsigned rb, unsigned re, unsigned cb, unsigned ce)
		: Super(make_pair(rb, re), make_pair(cb, ce)) {}

	// Const methods
	const unsigned & rowsBegin() const { return first.first; }
	const unsigned & rowsEnd() const { return first.second; }
	const unsigned & colsBegin() const { return second.first; }
	const unsigned & colsEnd() const { return second.second; }

	// Non-const methods
	unsigned & rowsBegin() { return first.first; }
	unsigned & rowsEnd() { return first.second; }
	unsigned & colsBegin() { return second.first; }
	unsigned & colsEnd() { return second.second; }
};

ostream & operator<<(ostream & os, const Zone & z);

// Breaks image into n^2 zones approximately equal zones
// Returns nxn matrix of zone ranges
Matrix<Zone> zoning(unsigned row, unsigned col, unsigned num);

#endif