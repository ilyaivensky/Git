/*                                                                 -*- C++ -*-
 * File: statistical_features.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 26, 2013
 *
 * Description:
 *   Implements statistical features extraction from binary image
 *   
 */


#ifndef _STATISTICAL_FEATURES_H_
#define _STATISTICAL_FEATURES_H_

#include "LA/matrix.h"
#include "LA/linear_algebra.h"
#include "PR/zoning.h"
#include "PR/utils.h"

template <class Image, class Feature>
vector<Feature> raw_binary_pixels(const Image & image)
{
	vector<Feature> features;

	for (Image::const_iterator itR = image.begin(), itREnd = image.end(); itR != itREnd; ++itR)
		features.insert(features.end(), itR->begin(), itR->end());

	return features;
}

template <class Image, class Feature>
vector<Feature> histograms(const Image & image, const Matrix<Zone> & zones)
{
	vector<Feature> features;

	// Calculate features separately for each zone
	for (Matrix<Zone>::const_iterator itRow = zones.begin(), itRowEnd = zones.end(); itRow != itRowEnd; ++itRow)
	{
		for (vector<Zone>::const_iterator itCol = itRow->begin(), itColEnd = itRow->end(); itCol != itColEnd; ++itCol)
		{
			const Zone & z  = *itCol;

			unsigned rows = z.rowsEnd() - z.rowsBegin();
			unsigned cols = z.colsEnd() - z.colsBegin();

			// Extracts vertical and horizontal histograms
			vector<unsigned> vp(cols), hp(rows);

			for (unsigned r = z.rowsBegin(), r_loc = 0; r < z.rowsEnd(); ++r, ++r_loc)
				for (unsigned c = z.colsBegin(), c_loc = 0; c < z.colsEnd(); ++c, ++c_loc)
					if (image[r][c])
						++hp[r_loc], ++vp[c_loc]; 

			features.insert(features.end(), hp.begin(), hp.end());
			features.insert(features.end(), vp.begin(), vp.end());
#if 1
			// Extracts left-right diagonal histograms
			vector<unsigned> lrd(rows + cols - 1);
			for (unsigned k = 0; k < lrd.size(); ++k)
			{
				// Set the break point in intioalization of row idx and col idx
				// Before break point the driver (i.e. the index defining diagonal)
				// is row idx, after - col idx
				// Therefore the break point is defined by rows
				unsigned bp = rows;

				unsigned r = k < bp ? z.rowsEnd() - k - 1 : z.rowsBegin();
				unsigned c = k < bp ? z.colsBegin() : z.colsBegin() + k - bp + 1;

				for (; r < z.rowsEnd() && c < z.colsEnd(); ++r, ++c)
				{
					if (image[r][c])
						++lrd[k];
				}
			}
			features.insert(features.end(), lrd.begin(), lrd.end());

			// Extracts right-left diagonal histograms
			vector<unsigned> rld(rows + cols - 1);
			for (unsigned k = 0; k < rld.size(); ++k)
			{
				// Set the break point in intioalization of row idx and col idx
				// Before break point the driver (i.e. the index defining diagonal)
				// is col idx, after - row idx
				// Therefore the break point is defined by cols
				unsigned bp = cols;

				unsigned r = k < bp ? z.rowsBegin() : z.rowsBegin() + k - bp + 1;
				unsigned c = k < bp ? z.colsBegin() + k + 1 : z.colsEnd();

				for (; r < z.rowsEnd() && c > z.colsBegin(); ++r, --c)
				{
					if (image[r][c - 1])
						++lrd[k];
				}
			}

			features.insert(features.end(), rld.begin(), rld.end());
#endif
		}
	}
	
	return features;
}

template <class Image, class Feature>
vector<Feature> radial_histograms(const Image & img)
{
	vector<Feature> features(32);

	unsigned total_r = 0, total_c = 0, total = 0;

	for (unsigned r = 0; r < img.nrow(); ++r)
		for (unsigned c = 0; c < img.ncol(); ++c)
		{
			if (img[r][c])
			{
				++total;
				total_r += r;
				total_c += c;
			}
		}

	float fmean_r = static_cast<float>(total_r) / total;
	float fmean_c = static_cast<float>(total_c) / total;

	unsigned granularity = 8;

	float step = (float) 1 / granularity;

	for (unsigned k = 0; k < granularity; ++k)
	{
		float k_step = k * step;

		for (float r = fmean_r, c = fmean_c; 
			static_cast<unsigned>(r + 0.5) < img.nrow() && static_cast<unsigned>(c + 0.5) < img.ncol(); 
			r = r + sqrt(1 - k_step), c = c + sqrt(k_step))
		{
			if (img[static_cast<unsigned>(r + 0.5)][static_cast<unsigned>(c + 0.5)])
				features[k] += 1;
		}
	}

	for (unsigned k = 0; k < granularity; ++k)
	{
		float k_step = k * step;

		for (float r = fmean_r, c = fmean_c; 
			static_cast<unsigned>(r + 0.5) < img.nrow() && static_cast<signed>(c + 0.5) >= 0; 
			r = r + sqrt(1 - k_step), c = c - sqrt(k_step))
		{
			if (img[static_cast<unsigned>(r + 0.5)][static_cast<unsigned>(c + 0.5)])
				features[k + granularity] += 1;
		}
	}

	for (unsigned k = 0; k < granularity; ++k)
	{
		float k_step = k * step;

		for (float r = fmean_r, c = fmean_c; 
			static_cast<signed>(r + 0.5) >= 0 && static_cast<unsigned>(c + 0.5) < img.ncol(); 
			r = r - sqrt(1 - k_step), c = c + sqrt(k_step))
		{
			if (img[static_cast<unsigned>(r + 0.5)][static_cast<unsigned>(c + 0.5)])
				features[k + (granularity * 2)] += 1;
		}
	}

	for (unsigned k = 0; k < granularity; ++k)
	{
		float k_step = k * step;

		for (float r = fmean_r, c = fmean_c; 
			static_cast<signed>(r + 0.5) >= 0 && static_cast<signed>(c + 0.5) >= 0; 
			r = r - sqrt(1 - k_step), c = c - sqrt(k_step))
		{
			if (img[static_cast<unsigned>(r + 0.5)][static_cast<unsigned>(c + 0.5)])
				features[k + (granularity * 3)] += 1;
		}
	}

	return features;
}

template <class Image, class Feature>
vector<Feature> _chain_codes(const Image & contour)
{
	vector<Feature> codes(8);

	for (unsigned r = 1; r < contour.nrow() - 1; ++r)
		for (unsigned c = 1; c < contour.ncol() - 1; ++c)
		{
			if (!contour[r][c])
				continue;
			
			if (!contour[r][c + 1])
			{
				if (contour[r - 1][c + 1])
					codes[1] += 1;
				else if (contour[r - 1][c])
					codes[2] += 1;
			}
			if (!contour[r - 1][c])
			{
				if (contour[r - 1][c - 1])
					codes[3] += 1;
				else if (contour[r][c - 1])
					codes[4] += 1;
			}
			if (!contour[r][c - 1])
			{
				if (contour[r + 1][c - 1])
					codes[5] += 1;
				else if (contour[r + 1][c])
					codes[6] += 1;
			}
			if (!contour[r + 1][c])
			{
				if (contour[r + 1][c + 1])
					codes[7] += 1;
				else if (contour[r][c + 1])
					codes[0] += 1;
			}
		}

	return codes;
}

template <class Image, class Feature>
vector<Feature> chain_codes(const Image & img)
{
	bool frame = true;
	for (unsigned c = 0; frame && c < img.ncol(); ++c)
		if (img[0][c] || img[img.nrow() - 1][c])
			frame = false;

	for (unsigned r = 0; frame && r < img.nrow(); ++r)
		if (img[r][0] || img[r][img.ncol() - 1])
			frame = false;

	if (frame)
		return _chain_codes<Image, Feature>(img);

	Image framed = make_frame(img);
	
	return _chain_codes<Image, Feature>(framed);
}

// For each zone in zones, calculates distance from global centroid to zone's centroid 
// Each distance to local centroid is accompanied by the number of black pixels in that zone
template <class Image, class Feature>
vector<Feature> fourier_centroid_distances(const Image & contour, const Matrix<Zone> & zones)
{
	vector<Feature> retval(zones.nrow() * zones.ncol() * 2);

	unsigned x_acc = 0, y_acc = 0, total = 0;
	
	for (unsigned x = 0; x < contour.nrow(); ++x)
		for (unsigned y = 0; y < contour.ncol(); ++y)
			if (contour[x][y])
				x_acc += x, y_acc += y, ++total;
				

	if (!total)
		return retval;

	float x_centr = static_cast<float>(x_acc) / total;
	float y_centr = static_cast<float>(y_acc) / total;

	vector<float> centr(2);
	centr[0] = x_centr, centr[1] = y_centr;

	// offset in the retval
	unsigned zone_offset = 0;

	for (unsigned zx = 0; zx < zones.nrow(); ++zx)
	{
		for (unsigned zy = 0; zy < zones.ncol(); ++zy, zone_offset += 2)
		{
			const Zone & z = zones[zx][zy];

			unsigned zx_acc = 0, zy_acc = 0, total_z = 0;

			for (unsigned x = z.rowsBegin(); x < z.rowsEnd(); ++x)
				for (unsigned y = z.colsBegin(); y < z.colsEnd(); ++y)
					if (contour[x][y])
						zx_acc += x, zy_acc += y, ++total_z;

			// This zone is empty
			if (!total_z) continue;

			float zx_centr = static_cast<float>(zx_acc) / total_z;
			float zy_centr = static_cast<float>(zy_acc) / total_z;

			vector<float> centr_z(2);
			centr_z[0] = zx_centr, centr_z[1] = zy_centr;
	
			float dist_z = square_dist(centr, centr_z);
			retval[zone_offset] = dist_z;
			retval[zone_offset + 1] = static_cast<float>(total_z);
		}
	}

	return retval;
}

template <class T>
float label(const T & r)
{
	return static_cast<float>(r > 0 ? 1 : 0);
}

template <class Image, class Feature>
vector<Feature> covarience(const Image & img)
{
	vector<Feature> features(3);

	// Create matrix of pixel coordinates 
	Matrix<Feature> pixels;
	for (unsigned r = 0; r < img.nrow(); ++r)
		for (unsigned c = 0; c < img.ncol(); ++c)
			if (img[r][c])
			{
				vector<Feature> pixel(2);
				pixel[0] = static_cast<Feature>(r), pixel[1] = static_cast<Feature>(c);
				pixels.add_row(pixel);
			}

	// Get covariance (2x2 matrix)
	Matrix<Feature> c = cov(pixels);

	// Put down to the feature vector (only 3 features since c is symmetric)
	features[0] = c[0][0], features[1] = c[0][1], features[2] = c[1][1];

	return features;
}

#endif
