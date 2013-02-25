#ifndef _SVM_H
#define _SVM_H

#include "matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"
#define  BUFSIZE 256

//#pragma comment (lib,"C:/Program Files (x86)/MATLAB/R2012a Student/extern/lib/win32/microsoft/libeng.lib")
//#pragma comment (lib,"C:/Program Files (x86)/MATLAB/R2012a Student/extern/lib/win32/microsoft/libmx.lib")
//#pragma comment (lib,"C:/Program Files (x86)/MATLAB/R2012a Student/extern/lib/win32/microsoft/libut.lib")

namespace SVM {

template <class T>
void learn(const Matrix<T> & x, const Matrix<T> & y)
{
	Engine *ep;
	mxArray *T = NULL, *result = NULL;
//	char buffer[BUFSIZE+1];
	double time[10] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };

	/*
	 * Call engOpen with a NULL string. This starts a MATLAB process 
     * on the current host using the command "matlab".
	 */
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return;// EXIT_FAILURE;
	}
	else
		fprintf(stderr, "\nStarted MATLAB engine\n");
}

} // namespace
#endif