/*                                                                 -*- C++ -*-
 * File: thinning.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 2, 2013
 *
 * Description:
 *   Implementation of thinning algorithms
 *   
 */

#ifndef _THINNING_H_
#define _THINNING_H_

#include "PR/utils.h"

template <class Image>
static Image _zhang_suen_thinning (const Image & img)
{
	Image skeleton(img);
	bool changes = true;
	while (changes) 
	{
		changes = false;
		{
			// Step 1
			Image contour(img.nrow(), img.ncol());
			for (unsigned r = 1; r < skeleton.nrow() - 1; ++r)
				for (unsigned c = 1; c < skeleton.ncol() - 1; ++c)
				{
					if (!skeleton[r][c])
						continue;

					unsigned b = b_score(skeleton, r, c);
					if (b < 2 || b > 6)
						continue;

					unsigned a = a_score(skeleton, r, c);
					if (a != 1)
						continue;

					if (skeleton[r - 1][c] && skeleton[r][c + 1] && skeleton[r + 1][c])
						continue;

					if (skeleton[r][c + 1] && skeleton[r + 1][c] && skeleton[r][c - 1])
						continue;

					if (skeleton[r][c])
						contour[r][c] = '*';
				}
			
			for (unsigned r = 1; r < contour.nrow() - 1; ++r)
				for (unsigned c = 1; c < contour.ncol() - 1; ++c)
				{
					if (!contour[r][c])
						continue;

					skeleton[r][c] = 0;
					changes = true;
				}
		}

		{
			// Step 2
			Image contour(img.nrow(), img.ncol());
			for (unsigned r = 1; r < skeleton.nrow() - 1; ++r)
				for (unsigned c = 1; c < skeleton.ncol() - 1; ++c)
				{
					if (!skeleton[r][c])
						continue;

					unsigned b = b_score(skeleton, r, c);
					if (b < 2 || b > 6)
						continue;

					unsigned a = a_score(skeleton, r, c);
					if (a != 1)
						continue;

					if (skeleton[r - 1][c] && skeleton[r][c + 1] && skeleton[r][c - 1])
						continue;

					if (skeleton[r - 1][c] && skeleton[r + 1][c] && skeleton[r][c - 1])
						continue;

					if (skeleton[r][c])
						contour[r][c] = '*';
				}
			
			for (unsigned r = 1; r < contour.nrow() - 1; ++r)
				for (unsigned c = 1; c < contour.ncol() - 1; ++c)
				{
					if (!contour[r][c])
						continue;

					skeleton[r][c] = 0;
					changes = true;
				}
		}
	}

	// Make contour binary
	for (unsigned r = 0; r < skeleton.nrow(); ++r)
		for (unsigned c = 0; c < skeleton.ncol(); ++c)
			if (skeleton[r][c])
				skeleton[r][c] = '*';

	return skeleton;
}

template <class Image>
Image zhang_suen_thinning (const Image & img)
{
	// Detect whether there is a white frame. 
	// It is needed for optimization only
	// If there is a frame, then we can scan image from 1 to last -1,
	// and check neighbors without checking boundaries

	bool frame = true;
	for (unsigned c = 0; frame && c < img.ncol(); ++c)
		if (img[0][c] || img[img.nrow() - 1][c])
			frame = false;

	for (unsigned r = 0; frame && r < img.nrow(); ++r)
		if (img[r][0] || img[r][img.ncol() - 1])
			frame = false;

	if (frame)
		return _zhang_suen_thinning(img);

	// There is no frame - we will create one temporarily 
	
	Image framed = make_frame(img);
	Image framed_skeleton = _zhang_suen_thinning(framed);
	Image skeleton(img.nrow(), img.ncol());
	for (unsigned r = 1; r < framed_skeleton.nrow() - 1; ++r)
		for (unsigned c = 1; c < framed_skeleton.ncol() - 1; ++c)
			skeleton[r - 1][c - 1] = framed_skeleton[r][c];

	return skeleton;
}

#endif
