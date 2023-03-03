// fms_dual.cpp - test dual numbers
#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>
#include "fms_dual.h"

using namespace fms;

template<class X>
constexpr int test_dual()
{
	{
		dual<X> d;
		assert(d == 0);
		dual d2{ d };
		assert(d2 == d);
		d = d2;
		assert(!(d != d2));
		assert(d == -d);
	}
	{
		dual<X> d(1);
		assert(d == 1);
	}
	{
		dual<X> d(2, 3);
		auto i_ = d * inv(d);
		assert(i_ == 1);
		auto _i = inv(d) * d;
		assert(_i == 1);
	}
	{
		dual<X> e(0, 1);
		assert(e != 0);
		assert(e * e == 0);
	}
	{
		assert(-dual<X>(1, 2) == dual<X>(-1, -2));

		assert(dual<X>(1, 2) + dual<X>(3,4) == dual<X>(4, 6));
		assert(X(1) + dual<X>(3, 4) == dual<X>(4, 4));
		assert(dual<X>(1, 2) + X(3) == dual<X>(4, 2));

		assert(dual<X>(1, 2) - dual<X>(3, 4) == dual<X>(-2, -2));
		assert(X(1) - dual<X>(3, 4) == dual<X>(-2, -4));
		assert(dual<X>(1, 2) - X(3) == dual<X>(-2, 2));

		assert(dual<X>(1, 2) * dual<X>(3, 4) == dual<X>(3, 10));
		assert(X(1) * dual<X>(3, 4) == dual<X>(3, 4));
		assert(dual<X>(1, 2) * X(3) == dual<X>(3, 6));

		auto i = inv(dual<X>(3, 4));
		assert(dual<X>(1, 2) / dual<X>(3, 4) == dual<X>(1, 2) * i);
		assert(X(1) / dual<X>(3, 4) == i);
		assert(dual<X>(1, 2) / X(3) == dual<X>(1./3, 2/3.));
	}

	return 0;
}
int test_dual_double = test_dual<double>();

template<class X>
X sq(X x)
{
	return x * x;
}
template<class X>
X dsq(X x)
{
	return 2 * x;
}

template<class X>
int test_derivative()
{
	{
		// derivative of sq
		auto Dsq = D(sq<fms::dual<X>>);

		for (X x = -2; x < 2; x += X(.1)) {
			assert(Dsq(x) == 2 * x);
		}
	}
	{
		auto _sq = _(sq<X>, dsq<X>);

		for (X x = -2; x < 2; x += X(.1)) {
			auto sqx = sq(dual(x, X(1)));
			assert(sqx._0 == x * x);
			assert(sqx._1 == 2 * x);
			auto _sqrx = _sq(fms::dual<X>(x, X(1)));
			assert(_sqrx._0 == sq(x));
			assert(_sqrx._1 == dsq(x));
		}
	}
	{
		auto _sq = _(sq<X>, dsq<X>);
		auto D_sq = D(_sq);
		assert(2*3 == D_sq(X(3)));
		auto D2_sq = D(D_sq);
		X x;
		x = D2_sq(X(3));
		x = x;
	}

	return 0;
}
int test_derivative_double = test_derivative<double>();
int test_derivative_float = test_derivative<float>();

template<class X>
inline X cu(const X& x)
{
	return x * x * x;
}
template<class X>
inline X dcu(const X& x)
{
	return X(3) * x * x;
}

// second derivative
template<class X>
int test_derivative2()
{
	{
		auto _cu = _(cu<X>, dcu<X>);
		dual<X> _x(1, 1);
		auto _cux = _cu(_x); // [1, 3; 0, 1]
		assert(_cux._0 == 1);
		assert(_cux._1 == 3);

		auto Dcu = D(cu<dual<X>>);
		for (const X& x : { X(-1), X(0), X(.1) }) {
			assert(Dcu(x) == dcu(x));
		}

		//auto D2cu = _(cu<dual<X>>, D(cu<dual<X>>));
		//dual<dual<X>> __x(_x, dual<X>(1));
		//auto __D2cu = D2cu(__x);
	}

	return 0;
}
int test_derivative2_double = test_derivative2<double>();

#define M_SQRT2PI 2.5066282746310005024157652848110452530069867406099383166299235763

// standard normal cumulative distribtion
template<class X>
inline auto N_0 = [](X x) { return X((1 + erf(x / M_SQRT2)) / 2); };

// stanard normal density
template<class X>
inline auto N_1 = [](X x) { return X(exp(-x * x / 2) / M_SQRT2PI); };

// standard normal cumulative for dual numbers
template<class X>
inline auto _N = _(N_0<X>, N_1<X>);

template<class X>
int test_log()
{
	return 0;
}
int test_log_double = test_log<double>();

template<class X>
int test_black()
{
	// F = f exp(s Z - s^2/2), Z standard normal
	// put value E[(k - F)^+ = k P(F <= k) - f P_s(F <= k), P_s(Z <= z) = P(Z <= z - s)

	X f = 100;
	X s = X(0.1);
	X k = 100;

	// F <= k iff Z <= log(k/f)/s + s/2
	auto moneyness = [&](dual<X> f) { return _log<X>(k / f) / s + s / X(2); };
	auto put = [&](dual<X> f) {
		auto z = moneyness(f);
		auto N = _N<X>(z);
		auto Ns = _N<X>(z - s);

		return k * N - f * Ns;
	};

	dual<X> v = put(dual<X>(f, X(1)));

	// (d/df) E[(k - F)^+] = E[1(F <= k) F/f] = P_s(F <= k)
	X z = log(k / f) / 2 + s / 2;
	X delta = -N_0<X>(z - s); // closed form

	X err = v._1 - delta;
	assert(abs(err) <= std::numeric_limits<X>::epsilon());
	
	return 0;
}
int test_black_double = test_black<double>();
int test_black_float = test_black<float>();

int main()
{
	return 0;
}
