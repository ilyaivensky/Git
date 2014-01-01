#include "LA/vector_utils.h"

vector<int> & operator /= (vector<int> & v, const int & scalar)
{
	for (auto & el : v)
		el /= scalar;

	return v;
}