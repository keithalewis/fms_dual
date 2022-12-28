// fms_dual.h - dual numbers  x0 + x1 e  where e != 0 but e^2 = 0
#pragma once
#include <concepts>
#include <functional>

namespace fms {

	// x = _0 + _1 e
	template<class X>
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

	// promote to dual function
	template<class F, class dF>
	class _ {
		F f;
		dF df;
	public:
		_(F f, dF df)
			: f{ f }, df{ df }
		{ }
		template<class X>
		fms::dual<X> operator()(fms::dual<X> x)
		{
			return fms::dual<X>(f(x._0), df(x._0)*x._1);
		}
	};

	// F: dual<X> -> dual<X>
	template<class F>
	class D {
		F f;
	public:
		D(F f)
			: f{ f }
		{ }
		template<class X>
		X operator()(X x) const
		{
			auto df = [&](X x) { return f(dual<X>(x, X(1))); };

			return df(x)._1;
		}
	};

} // namespace fms

template<typename X>
inline constexpr fms::dual<X> operator+(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 + y._0, x._1 + y._1);
}
template<typename X>
inline constexpr fms::dual<X> operator+(const X& x, const fms::dual<X>& y)
{
	return fms::dual(x + y._0, y._1);
}
template<typename X>
inline constexpr fms::dual<X> operator+(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0 + y, x._1);
}

// negative of a dual number
template<typename X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x)
{
	return fms::dual(-x._0, -x._1);
}

template<typename X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 - y._0, x._1 - y._1);
}
template<typename X>
inline constexpr fms::dual<X> operator-(const X& x, const fms::dual<X>& y)
{
	return fms::dual(x - y._0, -y._1);
}
template<typename X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0 - y, x._1);
}

template<typename X>
inline constexpr fms::dual<X> operator*(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 * y._0, x._0 * y._1 + x._1 * y._0);
}
template<typename X>
inline constexpr fms::dual<X> operator*(const X& x, const fms::dual<X>& y)
{
	return fms::dual(x * y._0, x * y._1);
}
template<typename X>
inline constexpr fms::dual<X> operator*(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0 * y, x._1 * y);
}

// inverse of a dual number, x*inv(x) == dual(1, 0) == inv(x)*x
template<typename X>
inline constexpr fms::dual<X> inv(const fms::dual<X>& x)
{
	return fms::dual(1 / x._0, -x._1 / (x._0 * x._0));
}

template<typename X>
inline constexpr fms::dual<X> operator/(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return x * inv(y);
}
template<typename X>
inline constexpr fms::dual<X> operator/(const X& x, const fms::dual<X>& y)
{
	return x * inv(y);
}
template<typename X>
inline constexpr fms::dual<X> operator/(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0/y, x._1/y);
}


// derivative of f at x, f(x + e) = f(x) + f'(x) e
template<class F, typename X>
inline constexpr X _1(F f, const X& x)
{
	auto _x = fms::dual<X>(x, 1);

	return f<fms::dual<X>(_x)._1;
}
