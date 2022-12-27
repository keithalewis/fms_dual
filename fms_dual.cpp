// fms_dual.cpp - test dual numbers
#include <cassert>
#include "fms_dual.h"

using namespace fms;

template<class X>
int test_dual()
{
	{
		dual<> d;
		assert(d == 0);
		dual d2{ d };
		assert(d2 == d);
		d = d2;
		assert(!(d != d2));
		assert(d == -d);
	}
	{
		dual<> d(1);
		assert(d == 1);
	}
	{
		dual<> d(2, 3);
		auto i = d * inv(d);
		assert(i == 1);
	}

	return 0;
}
int test_dual_double = test_dual<double>();

template<class X>
dual<X> sq(const dual<X>& x)
{
	return x * x;
}

template<class X>
int test_derivative()
{
	for (X x = -2; x < 2; x += .1) {
		auto sqx = sq(dual(x,1));
		assert(sqx._1 == 2 * x);
		//X dsqx = _1(sq, x);
	}

	return 0;
}
int test_derivative_double = test_derivative<double>();

int main()
{
	return 0;
}
