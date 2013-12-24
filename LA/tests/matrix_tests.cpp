/*                                                                 -*- C++ -*-
 * File: matrix_tests.cpp
 * 
 * Author: Ilya Ivensky
 * 
 * Created on: Dec 19, 2013
 *
 * Description:
 *   Boost unit tests for matrix
 *   
 */

#include <boost/assign/list_of.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "LA/matrix.h"
#include "LA/linear_algebra.h"

using namespace std;
using boost::unit_test_framework::test_suite;

void gramian_test()
{
	Matrix<int> m = { 
		{  24,   0,  30 },
		{  24,  30, -30 },
		{  -6,   0,   0 },
		{  -6,   0,  30 },
		{ -36, -30, -30 }
	};

	Matrix<int> g = {
		{ 2520, 1800,  900 },
		{ 1800, 1800,    0 },
		{  900,    0, 3600 }
	};
	
	BOOST_REQUIRE(gram(m) == g); 
}

void covariance_test()
{
	Matrix<int> m = { 
		{ 90, 60, 90 },
		{ 90, 90, 30 },
		{ 60, 60, 60 },
		{ 60, 60, 90 },
		{ 30, 30, 30 }
	};

	Matrix<int> c = { 
		{ 504, 360, 180 },
		{ 360, 360,   0 },
		{ 180,   0, 720 }
	};

	BOOST_REQUIRE(cov(m) == c); 
}

void inverse_test()
{
	Matrix<float> a;
	a.add_row({ 1, 0 });
	a.add_row({ 2, 2 });

	Matrix<float> a1;
	a1.add_row(boost::assign::list_of<float>( 1)(  0));
	a1.add_row(boost::assign::list_of<float>(-1)(0.5));

	BOOST_REQUIRE(inv(a) == a1);
}

void inverse_test2()
{
	Matrix<float> a(2, 2);
	a[0][0] = 1, a[0][1] = 0;
	a[1][0] = 2, a[1][1] = 2;

	Matrix<float> a1 = inv(a);

	BOOST_REQUIRE(a * a1 == Matrix<float>::diag(2, 1.0));
}

void determinant3_test()
{
	Matrix<float> m = {
		{ -2,  2, -3 },
		{ -1,  1,  3 },
		{  2,  0, -1 }
	};

	BOOST_REQUIRE(det(m) == 18);
}

void determinant4_test()
{
	Matrix<int> m = {
		{  1,  2,  3,  4 },
		{  5,  6,  7,  8 },
		{  9, 10, 11, 12 },
		{ 13, 14, 15, 16 }
	};

	BOOST_REQUIRE(det(m) == 0);
}

void determinant5_test()
{
	Matrix<double> m(
		{
			{ 1,  8, -9,  7,  5 },
			{ 0,  1,  0,  4,  4 },
			{ 0,  0,  1,  2,  5 },
			{ 0,  0,  0,  1, -5 },
			{ 0,  0,  0,  0,  1 }
		}
	);

	BOOST_REQUIRE(det(m) == 1);
}

void linear_solution_test()
{
	Matrix<float> a = {
		{ 1, 2, 2 },
		{ 2, 2, 2 },
		{ 2, 2, 1 }
	};

	Matrix<float> y = {
		{ 1 },
		{ 2 },
		{ 3 }
	};

	Matrix<float> x = {
		{ 1 },
		{ 1 },
		{ -1 }
	};

	BOOST_REQUIRE(linear_solution(a, y) == x);
}

template <class T>
void eigenvalue_2x2_test(const Matrix<T> & m)
{
	vector<float> eigenvalues = eigenvalues_2x2(m);

	Matrix<float> l0 = {
		{ eigenvalues[0], 0 },
		{ 0, eigenvalues[0] }
	};

	Matrix<float> l1 = {
		{ eigenvalues[1], 0 },
		{ 0, eigenvalues[1] }
	};

	BOOST_REQUIRE(det(l0 - m) == 0);
	BOOST_REQUIRE(det(l1 - m) == 0);
}

void eigenvalue_2x2_test1()
{
	Matrix<float> m = { 
		{ 1, 2 },
		{ 4, 3 }
	};

	eigenvalue_2x2_test(m);
}

void eigenvalue_2x2_test2()
{
	Matrix<float> m = {
		{ 1, -4 },
		{ 4, -7 }
	};

	eigenvalue_2x2_test(m);
}

void characteristic_polynomial_3x3_test()
{
	Matrix<double> m = {
		{ 3, 2, 4 },
		{ 2, 0, 2 },
		{ 4, 2, 3 }
	};

	vector<double> cp = { -1, 6, 15, 8 };
	BOOST_REQUIRE(characteristic_polynomial(m) == cp);
}

void square_dist_test()
{
	vector<int> v1 = { 1, 2, 3 };
	vector<int> v2 = { 0, 4, 1 };

	BOOST_REQUIRE(square_dist(v1, v2) == 9);
}

void norm_test()
{
	vector<double> v = { 3, 4 };
	BOOST_REQUIRE(norm(v) == 5);
}

boost::unit_test_framework::test_suite * init_unit_test_suite(int argc, char *argv[])
{
    test_suite* test = BOOST_TEST_SUITE("Matrix test suite");

    test->add(BOOST_TEST_CASE(&gramian_test));
	test->add(BOOST_TEST_CASE(&covariance_test));
	test->add(BOOST_TEST_CASE(&inverse_test));
	test->add(BOOST_TEST_CASE(&inverse_test2));
	test->add(BOOST_TEST_CASE(&determinant3_test));
	test->add(BOOST_TEST_CASE(&determinant4_test));
	test->add(BOOST_TEST_CASE(&determinant5_test));
	test->add(BOOST_TEST_CASE(&linear_solution_test));
	test->add(BOOST_TEST_CASE(&square_dist_test));
	test->add(BOOST_TEST_CASE(&norm_test));
	test->add(BOOST_TEST_CASE(&eigenvalue_2x2_test1));
	test->add(BOOST_TEST_CASE(&eigenvalue_2x2_test2));
	test->add(BOOST_TEST_CASE(&characteristic_polynomial_3x3_test));

    return test;
}

int run_test(int argc, char* argv[])
{
  boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
  return ::boost::unit_test::unit_test_main(init_func, argc, argv );
}
