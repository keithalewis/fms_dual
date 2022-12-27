// fms_dual.h - dual numbers  x0 + x1 e  where e != 0 but e^2 = 0
#pragma once
#include <concepts>
#include <functional>

namespace fms {

	// x = _0 + _1 e
	template<class X = double>
	struct dual {
		X _0, _1;

		constexpr dual(const X& _0 = 0, const X& _1 = 0)
			: _0(_0), _1(_1)
		{ }
		template<typename Y, typename Z>
		constexpr dual(const Y& _0 = 0, const Z& _1 = 0)
			: dual(X(_0), X(_1))
		{ }
		dual(const dual&) = default;
		dual& operator=(const dual&) = default;
		~dual()
		{ }

		constexpr bool operator==(const dual& x) const
		{
			return _0 == x._0 and _1 == x._1;
		}
		constexpr bool operator!=(const dual& x) const
		{
			return !operator==(x);
		}

	};

} // namespace fms

template<typename X>
inline constexpr fms::dual<X> operator+(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 + y._0, x._1 + y._1);
}
// negative of a dual number
template<typename X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x)
{
	return fms::dual(-x._0, -x._1);
}
template<typename X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return a + -b;
}
template<typename X>
inline constexpr fms::dual<X> operator*(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return fms::dual(a._0 * b._0, a._0 * b._1 + a._1 * b._0);
}

// inverse of a dual number, a*inv(a) == dual(1, 0) == inv(a)*a
template<typename X>
inline constexpr fms::dual<X> inv(const fms::dual<X>& a)
{
	return fms::dual(1 / a._0, -a._1 / (a._0 * a._0));
}

template<typename X>
inline constexpr fms::dual<X> operator/(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return a * inv(b);
}

// derivative of f at x, f(x + e) = f(x) + f'(x) e
template<class F, typename X>
inline constexpr X _1(const F& f, const X& x)
{
	return f(fms::dual<X>(x, 1))._1;
}
