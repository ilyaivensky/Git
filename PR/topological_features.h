/*                                                                 -*- C++ -*-
 * File: topological_features.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 3, 2013
 *
 * Description:
 *   Utilities to extract endpoints, branches and crossings
 *	 Parameter zones allows collecting features separately per zones 
 *   
 */


#ifndef _TOPOLOGICAL_FEATURES_H_
#define _TOPOLOGICAL_FEATURES_H_

#include "PR/zoning.h"
#include "PR/utils.h"

template <class Image, class Feature>
vector<Feature> feature_points(const Image & img, const Matrix<Zone> & zones)
{
	vector<Feature> features(zones.nrow() * zones.ncol() * 3);
	Image framed = make_frame(img);

	unsigned zone_offset = 0;

	// Calculate features separately for each zone
	for (Matrix<Zone>::const_iterator itRow = zones.begin(), itRowEnd = zones.end(); itRow != itRowEnd; ++itRow)
	{
		for (vector<Zone>::const_iterator itCol = itRow->begin(), itColEnd = itRow->end(); itCol != itColEnd; ++itCol, zone_offset += 3)
		{
			const Zone & z  = *itCol;

			unsigned rows = z.rowsEnd() - z.rowsBegin();
			unsigned cols = z.colsEnd() - z.colsBegin();

			unsigned eps = 0, branches = 0, crosses = 0;

			for (unsigned r = z.rowsBegin(); r < z.rowsEnd(); ++r)
				for (unsigned c = z.colsBegin(); c < z.colsEnd(); ++c)
				{
					if (!img[r][c])
						continue;

					unsigned b = b_score(framed, r + 1, c + 1); 
					if (b < 2) ++eps;
					if (b == 4) ++crosses;
					if (b == 3)
					{
						unsigned a = a_score(framed, r + 1, c + 1);
						if (a == 3) ++branches;
					}
				}

			features[zone_offset] = static_cast<float>(eps);
			features[zone_offset + 1] = static_cast<float>(branches);
			features[zone_offset + 2] = static_cast<float>(crosses);
		}
	}

	return features;
}
#endif
