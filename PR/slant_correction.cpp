#include "slant_correction.h"

double moment(const Matrix & m, unsigned p, unsigned q)
{
	double result = 0.0;

	for (unsigned y = 0;  y < m.row; ++y)
		for (unsigned x = 0; y < m.col; ++x)
			result += pow(x, p) * pow(y, q) * m[y][x];
			
	return result;
}

double central_moment(const Matrix & m, unsigned xc, unsigned yc, unsigned p, unsigned q)
{
	double result = 0.0;

	for (unsigned y = 0;  y < m.row; ++y)
		for (unsigned x = 0; x < m.col; ++x)
			result += pow(x - xc, p) * pow(y - yc, q) * m[y][x];
			
	return result;
}

Matrix slant_correction(const Matrix & m)
{
	cerr << "slantCorrection: m.row=" << m.row << " m.col=" << m.col << endl;
	double m00 = moment(m, 0, 0);
	double m10 = moment(m, 1, 0);
	double m01 = moment(m, 0, 1);

	cerr << "m00=" << m00 << " m10=" << m10 << " m01=" << m01 << endl;

	unsigned xc = unsigned(m10 / m00);
	unsigned yc = unsigned(m01 / m00);
	
	cerr << "xc=" << xc << " yc=" << yc << endl;
	
	double mu11 = central_moment(m, xc, yc, 1, 1);
	double mu02 = central_moment(m, xc, yc, 0, 2);
	
	cerr << "mu11=" << mu11 << " m02=" << mu02 << endl;
	
	double tan = -(mu11 / mu02); 
	
	Matrix result(m.row, m.col);
		
	for (unsigned y = 0; y < m.row; ++y)
	{
		for (unsigned x = 0; x < m.col; ++x)
		{
			signed y_tmp = y - yc;
			int x1 = x - (int)(y_tmp * tan);

			// If are we out of the range due to rotation?
			if (x1 < 0 || x1 >= (signed)result.col)
				continue;
			
			result[y][x1] = m[y][x];
		}
	}
	return result;
}
