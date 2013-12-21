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
	Matrix<int> m;
	m.add_row(boost::assign::list_of( 24)(  0)( 30));
	m.add_row(boost::assign::list_of( 24)( 30)(-30));
	m.add_row(boost::assign::list_of( -6)(  0)(  0));
	m.add_row(boost::assign::list_of( -6)(  0)( 30));
	m.add_row(boost::assign::list_of(-36)(-30)(-30));

	Matrix<int> g;
	g.add_row(boost::assign::list_of(2520)(1800)( 900));
	g.add_row(boost::assign::list_of(1800)(1800)(   0));
	g.add_row(boost::assign::list_of( 900)(   0)(3600));

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

void determinant3_test()
{
	Matrix<float> m;

	m.add_row(boost::assign::list_of<float>(-2)( 2)(-3));
	m.add_row(boost::assign::list_of<float>(-1)( 1)( 3));
	m.add_row(boost::assign::list_of<float>( 2)( 0)(-1));

	BOOST_REQUIRE(det(m) == 18);
}

void determinant4_test()
{
	Matrix<int> m;

	m.add_row(boost::assign::list_of( 1)( 2)( 3)( 4));
	m.add_row(boost::assign::list_of( 5)( 6)( 7)( 8));
	m.add_row(boost::assign::list_of( 9)(10)(11)(12));
	m.add_row(boost::assign::list_of(13)(14)(15)(16));
	
	BOOST_REQUIRE(det(m) == 0);
}

void determinant5_test()
{
	Matrix<double> m;

	m.add_row(boost::assign::list_of( 1)( 8)(-9)( 7)( 5));
	m.add_row(boost::assign::list_of( 0)( 1)( 0)( 4)( 4));
	m.add_row(boost::assign::list_of( 0)( 0)( 1)( 2)( 5));
	m.add_row(boost::assign::list_of( 0)( 0)( 0)( 1)(-5));
	m.add_row(boost::assign::list_of( 0)( 0)( 0)( 0)( 1));

	BOOST_REQUIRE(det(m) == 1);
}

void linear_solution_test()
{
	Matrix<float> a;
	a.add_row(boost::assign::list_of<float>( 1)( 2)( 2));
	a.add_row(boost::assign::list_of<float>( 2)( 2)( 2));
	a.add_row(boost::assign::list_of<float>( 2)( 2)( 1));

	Matrix<float> y(3, 1);
	y[0][0] = 1;
	y[1][0] = 2;
	y[2][0] = 3;

	Matrix<float> x(3, 1);
	x[0][0] =  1;
	x[1][0] =  1;
	x[2][0] = -1;
	
	BOOST_REQUIRE(linear_solution(a, y) == x);
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

    return test;
}

int run_test(int argc, char* argv[])
{
  boost::unit_test::init_unit_test_func init_func = &init_unit_test_suite;
  return ::boost::unit_test::unit_test_main(init_func, argc, argv );
}

