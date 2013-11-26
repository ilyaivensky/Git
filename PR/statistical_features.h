/*                                                                 -*- C++ -*-
 * File: statistical_features.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 26, 2013
 *
 * Description:
 *   Implements statistical features extraction
 *   
 */


#ifndef _STATISTICAL_FEATURES_H_
#define _STATISTICAL_FEATURES_H_

#include "LA/matrix.h"

template <class Image, class Feature>
vector<Feature> raw_binary_pixels(const Image & image)
{
	vector<Feature> features;

	for (Image::const_iterator itR = image.begin(), itREnd = image.end(); itR != itREnd; ++itR)
		features.insert(features.end(), itR->begin(), itR->end());

	return features;
}

template <class Image, class Feature>
vector<Feature> histograms(const Image & image)
{
	vector<Feature> features;

	// Extracts vertical and horizontal histograms
	vector<unsigned> vp(image.col), hp(image.row);

	for (unsigned r = 0; r < image.row; ++r)
		for (unsigned c = 0; c < image.col; ++c)
			if (image[r][c])
				++hp[r], ++vp[c]; 

	features.insert(features.end(), hp.begin(), hp.end());
	features.insert(features.end(), vp.begin(), vp.end());

	// Extracts left-right diagonal histograms
	vector<unsigned> lrd(image.row + image.col - 1);
	for (unsigned k = 0; k < lrd.size(); ++k)
	{
		for (unsigned r = (k < image.row ? k : 0), c = (k < image.row ? 0 : k - image.row + 1); 
			r < image.row && c < image.col; ++r, ++c)
		{
			if (image[r][c])
				++lrd[k];
		}
	}
	features.insert(features.end(), lrd.begin(), lrd.end());

	// Extracts right-left diagonal histograms
	vector<unsigned> rld(image.row + image.col - 1);
	for (unsigned k = 0; k < rld.size(); ++k)
	{
		for (unsigned r = (k < image.row ? k : 0), c = (k < image.row ? image.col : k - image.col - 1); 
			r < image.row && c > 0; ++r, --c)
		{
			if (image[r][c - 1])
				++lrd[k];
		}
	}
	features.insert(features.end(), rld.begin(), rld.end());

	return features;
}

#endif