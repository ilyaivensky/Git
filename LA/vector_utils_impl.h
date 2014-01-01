#include <cstdlib>

template <class T>
bool is_zero(const vector<T> & v)
{
	static const T zero = T();
	for (const auto & el : v)
	{
		if (el != zero)
			return false;
	}

	return true;
}

template <class T>
T inner_product(const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("inner_product: v1.size() != v2.size()");

	T retval = 0.0;

	for (unsigned i = 0; i < v1.size(); ++i)
		retval += (v1[i] * v2[i]);

	return retval;
}

template <class T>
Matrix<T> outer_product(const vector<T> & v1, const vector<T> & v2)
{
	return Matrix<T>(v1, v2);
}

template <class T>
T square_dist(const vector<T> & v1, const vector<T> & v2)
{
	if (v1.size() != v2.size())
		throw exception("dist: v1.size() != v2.size()");

	T result = 0;
	for (auto it1 = v1.begin(), it1End = v1.end(),
		it2 = v2.begin(); it1 != it1End; ++it1, ++it2)
		result += static_cast<T>(pow(*it1 - *it2, 2));

	return result;
}

template <class T>
T euclidian_dist(const vector<T> & v1, const vector<T> & v2)
{
	return sqrt(square_dist(v1, v2));
}

double norm(const vector<int> & v, unsigned p)
{
	if (p < 1)
		throw exception("norm of vector is not defined for p < 1");

	double sum = 0;
	for (const auto & val : v)
		sum += pow(val, p);

	return pow(sum, static_cast<double>(1) / p);
}

template <class T>
T norm(const vector<T> & v, unsigned p)
{
	if (p < 1)
		throw exception("norm of vector is not defined for p < 1");

	T sum = 0;
	for (const auto & val : v)
		sum += static_cast<T>(pow(val, p));

	return pow(sum, static_cast<T>(1) / p);
}

template <class T>
T cosine(const vector<T> & v1, const vector<T> & v2)
{
	if (!is_zero(v1) && v1 == v2)
		return 1;

	return (v1 * v2) / (norm(v1) * norm(v2));
}

template <class T>
void make_vector_set(vector<T> & v)
{
	sort(v.begin(), v.end());
	auto itEnd = unique(v.begin(), v.end());
	v.erase(itEnd, v.end());
}

template <class T>
T operator * (const vector<T> & v1, const vector<T> & v2)
{
	return inner_product(v1, v2);
}

template <class T>
vector<T> & operator *= (vector<T> & v, const T & scalar)
{
	for (auto & el : v)
		el *= scalar;

	return v;
}

template <class T>
vector<T> operator * (const vector<T> & v, const T & scalar)
{
	vector<T> tmp = v;
	tmp *= scalar;

	return tmp;
}

template <>
vector<long> & operator /= (vector<long> & v, const long & scalar)
{
	for (auto & el : v)
	{
		el /= scalar;
	}

	return v;
}

template <>
vector<int> & operator /= (vector<int> & v, const int & scalar)
{
	for (auto & el : v)
	{
		el /= scalar;
	}

	return v;
}

template <class T>
vector<T> & operator /= (vector<T> & v, const T & scalar)
{
	for (auto & el : v)
	{
		el /= scalar;
		auto rounded = round(el);
		if (fabs(el - rounded) < EPSILON) el = rounded;
	}

	return v;
}

template <class T>
vector<T> operator / (const vector<T> & v, const T & scalar)
{
	vector<T> tmp = v;
	tmp /= scalar;

	return tmp;
}

template <class T>
vector<T> & operator += (vector<T> & v1, const vector<T> & v2)
{
	auto it2 = v2.begin();
	for (auto it1 = v1.begin(), it1End = v1.end(); it1 != it1End; ++it1, ++it2)
		*it1 += *it2;

	return v1;
}

template <class T>
vector<T> & operator -= (vector<T> & v1, const vector<T> & v2)
{
	auto it2 = v2.begin();
	for (auto it1 = v1.begin(), it1End = v1.end(); it1 != it1End; ++it1, ++it2)
		*it1 -= *it2;

	return v1;
}
template <class T>
ostream & operator << (ostream & os, const vector<T> & v)
{
	std::copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
	return os;
}
