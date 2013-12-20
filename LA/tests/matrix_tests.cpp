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
	Matrix<int> m(5, 3);
	m[0][0] =  24, m[0][1] =   0, m[0][2] =  30;
	m[1][0] =  24, m[1][1] =  30, m[1][2] = -30;
	m[2][0] =  -6, m[2][1] =   0, m[2][2] =   0;
	m[3][0] =  -6, m[3][1] =   0, m[3][2] =  30;
	m[4][0] = -36, m[4][1] = -30, m[4][2] = -30;

	Matrix<int> g(3, 3);
	g[0][0] = 2520, g[0][1] = 1800, g[0][2] =  900;
	g[1][0] = 1800, g[1][1] = 1800, g[1][2] =    0;
	g[2][0] =  900, g[2][1] =    0, g[2][2] = 3600;

	BOOST_REQUIRE(gram(m) == g); 
}

void covariance_test()
{
	Matrix<int> m(5, 3);
	m[0][0] = 90, m[0][1] = 60, m[0][2] = 90;
	m[1][0] = 90, m[1][1] = 90, m[1][2] = 30;
	m[2][0] = 60, m[2][1] = 60, m[2][2] = 60;
	m[3][0] = 60, m[3][1] = 60, m[3][2] = 90;
	m[4][0] = 30, m[4][1] = 30, m[4][2] = 30;

	Matrix<int> c(3, 3);
	c[0][0] = 504, c[0][1] = 360, c[0][2] = 180;
	c[1][0] = 360, c[1][1] = 360, c[1][2] =   0;
	c[2][0] = 180, c[2][1] =   0, c[2][2] = 720;

	BOOST_REQUIRE(cov(m) == c); 
}

void inverse_test()
{
	Matrix<float> a(2, 2);
	a[0][0] = 1, a[0][1] = 0;
	a[1][0] = 2, a[1][1] = 2;

	Matrix<float> a1(2, 2);
	a1[0][0] =  1, a1[0][1] =   0;
	a1[1][0] = -1, a1[1][1] = 0.5;

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

boost::unit_test_framework::test_suite * init_unit_test_suite(int argc, char *argv[])
{
    test_suite* test = BOOST_TEST_SUITE("Matrix test suite");

    test->add(BOOST_TEST_CASE(&gramian_test));
	test->add(BOOST_TEST_CASE(&covariance_test));
	test->add(BOOST_TEST_CASE(&inverse_test));
	test->add(BOOST_TEST_CASE(&inverse_test2));

    return test;
}

int run_test(int argc, char* argv[])
{
  boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
  return ::boost::unit_test::unit_test_main(init_func, argc, argv );
}

