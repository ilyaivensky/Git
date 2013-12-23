/*                                                                 -*- C++ -*-
 * File: contour_tracing.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Oct, 2013
 *
 * Description:
 *   Implements various algorithms of extracting contours
 *   
 */


#ifndef _CONTOUR_TRACING_H_
#define _CONTOUR_TRACING_H_

#include "LA/matrix.h"

template <class Image>
Image scan_contours(const Image & img);

template <class T>
Matrix<T> trace_contours(const Matrix<T> & img);

/******************************************************************************
* 
* IMPLEMENTATIONS
*
*******************************************************************************/

template <class Image>
Image scan_contours(const Image & m)
{
	Image contour(m.nrow(), m.ncol());
	
	// Scan left-right each row
	for (unsigned i = 0; i < m.nrow(); ++i)
	{	
		bool state = false;
		for (unsigned j = 0; j < m.ncol(); ++j)
		{
			if (bool(m[i][j]) != state)
			{
				// Record new state
				state = static_cast<bool>(m[i][j]);
				if (m[i][j])
				{
					// Add contour for the current field
					contour[i][j] = '*';
				}	
				else
				{
					// Add contour for the field on the prev col, same row
					contour[i][j - 1] = '*';
				}
			}
		}
	}

	// Scan up-down each col
	for (unsigned j = 0; j < m.ncol(); ++j)
	{	
		bool state = false;
		for (unsigned i = 0; i < m.nrow(); ++i)
		{
			if (bool(m[i][j]) != state)
			{
				// Record new state
				state = bool(m[i][j]);
				if (m[i][j])
				{
					// Add contour for the current field
					contour[i][j] = '*';
				}	
				else
				{
					// Add contour for the field on the prev row, same col
					contour[i - 1][j] = '*';
				}
			}
		}
	}

	return contour;
}

template <class T>
Matrix<T> bug_walk(const Matrix<T> & m, signed start_i, signed start_j)
{
	Matrix<T> contour(m.nrow(), m.ncol());

	enum Directions {
		NORTH = 0,
		EAST,
		SOUTH,
		WEST,
		LAST
	};

	signed i = start_i, j = start_j;

	if (!m[i][j])
		throw exception("invalid start");
 
	contour[i][j] = m[i][j];

	unsigned dir = WEST;

	while (true)
	{
		switch (dir)
		{
		case NORTH:
			--i;
			break;
		case EAST:
			++j;
			break;
		case SOUTH:
			++i;
			break;
		case WEST:
			--j;
			break;
		default:
			throw exception("unknown direction");
		}

		// Check whether it is a black pixel or white. Treat out-of-range pixels as white
		if (i >= 0 && i < (signed)m.nrow() && j >= 0 && j < (signed)m.ncol() && m[i][j])
		{
			// This is a black pixel. Record it
			contour[i][j] = m[i][j];
			// ... and make a left turn
			if (dir == NORTH) dir = LAST;
			--dir;
		}
		else
		{
			// This is a white pixel. Make a right turn
			++dir; 
			if (dir == LAST) dir = NORTH;
		}

		// Check for stop condition: 
		// 1. We should be on the starting pixel
		// 2. We should be going in the same direction as when we started
		if (i == start_i && j == start_j && dir == WEST)
			break;
	}

	return contour;
}

template <class T>
Matrix<T> trace_contours(const Matrix<T> & m)
{
	// Scan from bottom to top and from left to right to find a black pixel to start bug_walk from it
	signed i = m.nrow(), j = 0;
	for (; i > 0; --i)
		for (j = 0; j < (signed)m.ncol(); ++j)
			if (m[i - 1][j])
				return bug_walk(m, --i, j);

	// No black pixel was found. Return empty contour
	return Matrix<T>(m.nrow(), m.ncol());
}

#endif
