#ifndef _SLANT_CORRECTION_H_
#define _SLANT_CORRECTION_H_

#include "LA/matrix.h"

template <class T>
double moment(const Matrix<T> & m, unsigned p, unsigned q)
{
	double result = 0.0;

	for (unsigned y = 0;  y < m.row; ++y)
		for (unsigned x = 0; x < m.col; ++x)
			if (m[y][x] != 0)
				result += pow(x, p) * pow(y, q);
			
	return result;
}

template <class T>
double central_moment(const Matrix<T> & m, unsigned xc, unsigned yc, unsigned p, unsigned q)
{
	double result = 0.0;

	for (unsigned y = 0;  y < m.row; ++y)
		for (unsigned x = 0; x < m.col; ++x)
			if (m[y][x] != 0)
				result += pow(x - xc, p) * pow(y - yc, q);
			
	return result;
}

template <class T>
Matrix<T> slant_correction(const Matrix<T> & m)
{
	//cerr << "slantCorrection: m.row=" << m.row << " m.col=" << m.col << endl;
	double m00 = moment(m, 0, 0);
	double m10 = moment(m, 1, 0);
	double m01 = moment(m, 0, 1);

	//cerr << "m00=" << m00 << " m10=" << m10 << " m01=" << m01 << endl;

	unsigned xc = static_cast<unsigned>(m10 / m00 + 0.5);
	unsigned yc = static_cast<unsigned>(m01 / m00 + 0.5);
	
	//cerr << "xc=" << xc << " yc=" << yc << endl;
	
	double mu11 = central_moment(m, xc, yc, 1, 1);
	double mu02 = central_moment(m, xc, yc, 0, 2);
	
	//cerr << "mu11=" << mu11 << " m02=" << mu02 << endl;
	
	double tan = -(mu11 / mu02); 
	
	Matrix<T> result(m.row, m.col);
		
	for (unsigned y = 0; y < m.row; ++y)
	{
		for (unsigned x = 0; x < m.col; ++x)
		{
			signed y_tmp = y - yc;
			signed x1 = x - static_cast<signed>(y_tmp * tan + 0.5);

			// If are we out of the range due to rotation?
			if (x1 < 0 || x1 >= static_cast<signed>(result.col))
				continue;
			
			result[y][x1] = m[y][x];
		}
	}

	return result;
}

#endif