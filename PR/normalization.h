#ifndef _NORMALIZATION_H_
#define _NORMALIZATION_H_

#include "LA/matrix.h"

template <class T>
Matrix<T> linear_normalization(const Matrix<T> & image, unsigned size);

// Aspect Ratio Adaptive Normalization
template <class T>
Matrix<T> linear_aran_normalization(const Matrix<T> & image, unsigned size);

/******************************************************************************
* 
* IMPLEMENTATIONS
*
*******************************************************************************/

template <class T>
Matrix<T> linear_aran_normalization(const Matrix<T> & image, unsigned size)
{
	double k = static_cast<double>(size) / std::max(image.row, image.col);

	double dsize = static_cast<double>(size);
	double drow = static_cast<double>(image.row);
	double dcol = static_cast<double>(image.col);

    // Adding 0.5 in order to correctly round contious value (double) into discrete value (integer) 
	unsigned vert_shift = static_cast<unsigned>((dsize - drow * k) / 2 + 0.5);
	unsigned horiz_shift = static_cast<unsigned>((dsize - dcol * k) / 2 + 0.5);

	Matrix<T> tmp(size);

	// Read the image into standard pane using forward mapping
	for (unsigned i = 0; i < image.row; ++i)
		for (unsigned j = 0; j < image.col; ++j)
		{
			unsigned i1 = static_cast<unsigned>(i * k + 0.5) + vert_shift;
			unsigned j1 = static_cast<unsigned>(j * k + 0.5) + horiz_shift;

			if (i1 < size && j1 < size)
				tmp[i1][j1] = image[i][j]; 
		}

	return tmp;
}

template <class T>
Matrix<T> linear_normalization(const Matrix<T> & image, unsigned size)
{
	Matrix tmp(size);

	double alpha = (static_cast<double>(size)) / image.col;
	double beta = (static_cast<double>(size)) / image.row;

	// Read the image into standard pane using forward mapping
	for (unsigned i = 0; i < image.row; ++i)
		for (unsigned j = 0; j < image.col; ++j)
		{
			unsigned i1 = static_cast<unsigned>(i * beta + 0.5);
			unsigned j1 = static_cast<unsigned>(j * alpha + 0.5);

			if (i1 < size && j1 < size)
				tmp[i1][j1] = image[i][j]; 
		}

	return tmp;
}

#endif