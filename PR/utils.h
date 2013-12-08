/*                                                                 -*- C++ -*-
 * File: utils.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 3, 2013
 *
 * Description:
 *   Utilities for image pricessing
 *   
 */


#ifndef _UTILS_H_
#define _UTILS_H_

template <class Image>
Image make_frame(const Image & img)
{
	Image framed(img.nrow() + 2, img.ncol() + 2);
	for (unsigned r = 0; r < img.nrow(); ++r)
		for (unsigned c = 0; c < img.ncol(); ++c)
			framed[r + 1][c + 1] = img[r][c];

	return framed;
}

// Calculates B score of the pixel iterating counter the clock
template <class Image>
unsigned b_score(const Image & img, unsigned r, unsigned c)
{
	if (r == 0 || r == img.nrow() - 1 || 
		c == 0 || c == img.ncol() - 1)
		throw exception("First/last row/column");

	unsigned b = 0;
	if (img[r][c - 1])
		++b;
	if (img[r + 1][c - 1])
		++b;
	if (img[r + 1][c])
		++b;
	if (img[r + 1][c + 1])
		++b;
	if (img[r][c + 1])
		++b;
	if (img[r - 1][c + 1])
		++b;
	if (img[r - 1][c])
		++b;
	if (img[r - 1][c - 1])
		++b;

	return b;
}

// Calculates A score of the pixel iterating clockwise
template <class Image>
unsigned a_score(const Image & img, unsigned r, unsigned c)
{
	if (r == 0 || r == img.nrow() - 1 || 
		c == 0 || c == img.ncol() - 1)
		throw exception("First/last row/column");

	enum {
		WHITE = 0,
		BLACK = 1
	};

	bool state = WHITE;
	unsigned a = 0;
	if (img[r - 1][c])
		state = BLACK;
	if (img[r - 1][c + 1])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r][c + 1])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r + 1][c + 1])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r + 1][c])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r + 1][c - 1])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r][c - 1])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r - 1][c - 1])
	{
		if (state == WHITE) ++a;
		state = BLACK;
	}
	else
		state = WHITE;
	if (img[r - 1][c])
		if (state == WHITE) ++a;

	return a;
}

#endif
