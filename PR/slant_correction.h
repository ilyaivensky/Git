/*                                                                 -*- C++ -*-
 * File: slant_correction.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 5, 2013
 *
 * Description:
 *   Implements various algorithms of correnting slant
 *   
 */


#ifndef _SLANT_CORRECTION_H_
#define _SLANT_CORRECTION_H_

#include "LA/matrix.h"

template <class T>
double moment(const Matrix<T> & m, unsigned p, unsigned q)
{
	double result = 0.0;

	for (unsigned y = 0;  y < m.nrow(); ++y)
		for (unsigned x = 0; x < m.ncol(); ++x)
			if (m[y][x] != 0)
				result += pow(x, p) * pow(y, q);
			
	return result;
}

template <class T>
double central_moment(const Matrix<T> & m, unsigned xc, unsigned yc, unsigned p, unsigned q)
{
	double result = 0.0;

	for (unsigned y = 0;  y < m.nrow(); ++y)
		for (unsigned x = 0; x < m.ncol(); ++x)
			if (m[y][x] != 0)
				result += pow(x - xc, p) * pow(y - yc, q);
			
	return result;
}

template <class T>
Matrix<T> moment_based_normalization(const Matrix<T> & m)
{
	//cerr << "moment_based_normalization: m.row=" << m.row << " m.col=" << m.col << endl;
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
	
	Matrix<T> result(m.nrow(), m.ncol());
		
	for (unsigned y = 0; y < m.nrow(); ++y)
	{
		for (unsigned x = 0; x < m.ncol(); ++x)
		{
			signed y_tmp = y - yc;
			signed x1 = x - static_cast<signed>(y_tmp * tan + 0.5);

			// If are we out of the range due to rotation?
			if (x1 < 0 || x1 >= static_cast<signed>(result.ncol()))
				continue;
			
			result[y][x1] = m[y][x];
		}
	}

	return result;
}

template <class Image>
Image slant_correction(const Image & m)
{
	// Step 1. Calculate bounding diagonals

	unsigned b1 = 0;
	bool found_b1 = false;
	for (; b1 < m.nrow() + m.ncol() - 1 && !found_b1; ++b1)
	{
		unsigned r = (b1 > m.nrow() ? m.nrow() : b1 + 1);
		unsigned c = (b1 < m.nrow() ? 0 : b1 - m.nrow() + 1);

		// scan right_up-left_down diagonal
		for (; r > 0 && c < m.ncol() && !found_b1; --r, ++c)
		{
			if (m[r - 1][c])
				found_b1 = true;
		}
	}

	unsigned b2 = m.nrow() + m.ncol() - 1;
	bool found_b2 = false;
	for (; b2 > 0 && !found_b2; --b2)
	{
		unsigned r = (b2 > m.nrow() ? m.nrow() : b2);
		unsigned c = (b2 > m.nrow() ? b2 - m.nrow() : 0);

		// scan right_up-left_down diagonal
		for (; r > 0 && c < m.ncol() && !found_b2; --r, ++c)
		{
			if (m[r - 1][c])
				found_b2 = true;
		}
	}

	unsigned b3 = 0;
	bool found_b3 = false;
	for (; b3 < m.nrow() + m.ncol() - 1 && !found_b3; ++b3)
	{
		unsigned r = (b3 < m.ncol() ? 0 : b3 - m.ncol() + 1);
		unsigned c = (b3 < m.ncol() ? m.ncol() - b3 - 1 : 0);
		
		// scan left_down-right_up diagonal
		for (; r < m.nrow() && c < m.ncol() && !found_b3; ++r, ++c)
		{
			if (m[r][c])
				found_b3 = true;
		}
	}

	unsigned b4 = m.nrow() + m.ncol() - 1;
	bool found_b4 = false;
	for (; b4 > 0 && !found_b4; --b4)
	{
		unsigned r = (b4 > m.ncol() ? b4 - m.ncol() : 0);
		unsigned c = (b4 > m.ncol() ? 0 : m.ncol() - b4);

		// scan left_down-right_up diagonal
		for (; r < m.nrow() && c < m.ncol() && !found_b4; ++r, ++c)
		{
			if (m[r][c])
				found_b4 = true;
		}
	}

	// Step 2: Calculate angle of slant

	// Calculate A
	unsigned min_r = m.nrow(), max_r = 0;

	for (unsigned r = 0; r < m.nrow(); ++r)
		for (unsigned c = 0; c < m.ncol(); ++c)
			if (m[r][c])
			{
				if (r > max_r) max_r = r;
				if (r < min_r) min_r = r;
			}

	double a = max_r - min_r;

	// Calculate B
	signed b = b4 + b1 - b3 - b2;
	double b_double = static_cast<double>(b) / 2;

	double theta = b_double / a;

	// Step 3. Apply correction

	Image sl_corrected(m.nrow(), m.ncol());

	for (unsigned r = 0; r < m.nrow(); ++r)
		for (unsigned c = 0; c < m.ncol(); ++c)
			if (m[r][c])
			{
				signed c_res = static_cast<signed>(static_cast<double>(c) + r * theta + 0.5);
				if (c_res >= 0 && c_res < static_cast<signed>(m.ncol()))
					sl_corrected[r][c_res] = m[r][c];
			}
	
	// Step 4. Align with the center of frame

	unsigned r_total = 0, c_total = 0, total = 0;
	
	for (unsigned r = 0; r < sl_corrected.nrow(); ++r)
		for (unsigned c = 0; c < sl_corrected.ncol(); ++c)
			if (sl_corrected[r][c])
				r_total += r, c_total += c, ++total;

	double r_center = static_cast<double>(r_total);
	r_center /= total;

	double c_center = static_cast<double>(c_total);
	c_center /= total;

	signed r_shift = static_cast<signed>(r_center - m.nrow() / 2);
	signed c_shift = static_cast<signed>(c_center - m.ncol() / 2);

	Image result(m.nrow(), m.ncol());
	for (unsigned r = 0; r < sl_corrected.nrow(); ++r)
		for (unsigned c = 0; c < sl_corrected.ncol(); ++c)
			if (sl_corrected[r][c])
			{
				signed res_r = r - r_shift;
				signed res_c = c - c_shift;

				if (res_r >= 0 && res_r < static_cast<signed>(result.nrow()) &&
					res_c >= 0 && res_c < static_cast<signed>(result.ncol()))
					result[res_r][res_c] = sl_corrected[r][c];
			}

	return result;
}

#endif
