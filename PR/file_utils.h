#ifndef _FILE_UTILS_H_
#define _FILE_UTILS_H_

#include <string>
#include "matrix.h"

void writeFile(const std::string & fname, const Matrix & m);
Matrix readFile(const std::string & fname);

#endif