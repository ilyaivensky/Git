#ifndef _CONTOUR_TRACING_H_
#define _CONTOUR_TRACING_H_

#include "LA/matrix.h"

template <class T>
Matrix<T> scan_contours(const Matrix<T> & img);

template <class T>
Matrix<T> trace_contours(const Matrix<T> & img);

template <class T>
Matrix<T> scan_contours(const Matrix<T> & m)
{
	Matrix<T> contour(m.row, m.col);
	
	// Scan left-right each row
	for (unsigned i = 0; i < m.row; ++i)
	{	
		unsigned state = 0;
		for (unsigned j = 0; j < m.col; ++j)
		{
			if (m[i][j] != state)
			{
				// Record new state
				state = m[i][j];
				if (m[i][j])
				{
					// Add contour for the current field
					contour[i][j] = m[i][j];
				}	
				else
				{
					// Add contour for the field on the prev col, same row
					contour[i][j - 1] = m[i][j - 1];
				}
			}
		}
	}

	// Scan up-down each col
	for (unsigned j = 0; j < m.col; ++j)
	{	
		unsigned state = 0;
		for (unsigned i = 0; i < m.row; ++i)
		{
			if (m[i][j] != state)
			{
				// Record new state
				state = m[i][j];
				if (m[i][j])
				{
					// Add contour for the current field
					contour[i][j] = m[i][j];
				}	
				else
				{
					// Add contour for the field on the prev row, same col
					contour[i - 1][j] = m[i - 1][j];
				}
			}
		}
	}

	return contour;
}

template <class T>
Matrix<T> bug_walk(const Matrix<T> & m, signed start_i, signed start_j)
{
	Matrix<T> contour(m.row, m.col);

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
		if (i >= 0 && i < (signed)m.row && j >= 0 && j < (signed)m.col && m[i][j])
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
	signed i = m.row, j = 0;
	for (; i > 0; --i)
		for (j = 0; j < (signed)m.col; ++j)
			if (m[i - 1][j])
				return bug_walk(m, --i, j, fname);

	// No black pixel was found. Return empty contour
	return Matrix<T>(m.row, m.col);
}

#endif