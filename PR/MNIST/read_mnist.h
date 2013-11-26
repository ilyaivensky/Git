/*                                                                 -*- C++ -*-
 * File: read_mnist.h
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Nov 26, 2013
 *
 * Description:
 *   Utility file to extract data from a pair (image, label) raw MNIST data files
 *   
 */

#ifndef _READ_MNIST_H_
#define _READ_MNIST_H_

#include <cstdint>
#include <limits>
#include "LA/matrix.h"

namespace MNIST {

typedef uint8_t Pixel;
typedef Matrix<Pixel> Image;

int read_mnist(const std::string & imagesFile, const std::string & labelsFile, 
			 vector<Image> & images, vector<int> & labels, const std::string & outputFile, 
			 unsigned maxImgs = std::numeric_limits<unsigned>::max());

} // namespace

#endif