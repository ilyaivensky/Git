#include "matrix.h"

using namespace std;

ostream & operator<<(ostream & os, const vector<unsigned> & v)
{
	for (vector<unsigned>::const_iterator it = v.begin(), itEnd = v.end(); it != itEnd; ++it)
		cerr << *it << "\t";

	return os;
}

ostream & operator<<(ostream & os, const Matrix & m)
{
	for (Matrix::const_iterator it = m.begin(), itEnd = m.end(); it != itEnd; ++it)
		os << *it << endl;
     
	return os;
}