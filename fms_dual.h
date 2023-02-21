// fms_dual.h - dual numbers  x0 + x1 e  where e != 0 but e^2 = 0
#pragma once
#include <cmath>
#include <concepts>

namespace fms {

	// x = _0 + _1 e where e != 0 and e*e = 0
	template<class X>
	struct dual {
		X _0, _1;
		constexpr dual(const X& _0 = 0, const X& _1 = 0) noexcept
			: _0(_0), _1(_1)
		{ }
		template<class Y, class Z>
		constexpr dual(const Y& _0, const Z& _1)
			: dual(X(_0), X(_1))
		{ }
		dual(const dual&) = default;
		dual& operator=(const dual&) = default;
		~dual()
		{ }

		constexpr bool operator==(const dual& x) const = default;

		X norm2() const
		{
			return _0 * _0 + _1 * _1;
		}
	};

	template<class X>
	inline X norm2(const dual<X>& x)
	{
		return x.norm2();
	}

	// operator norm || [x_0 x_1; 0 x_0] ||
	template<class X>
	inline X norm(const dual<X>& x)
	{
		X c = 2 * x._0 * x._0 + x._1 * x._1;
		X d = sqrt(2 * x._0 * x._0 + c);

		return std::max(sqrt(-x._1 * d + c), sqrt(x._1 * d) + c) / sqrt(X(2));
	}

	// promote to dual function given function and first derivative
	template<class F, class dF>
	class _ {
		F f;
		dF df;
	public:
		_() { }
		_(F f, dF df)
			: f{ f }, df{ df }
		{ }
		// f(x0 + x1 e) = f(x0) + f'(x0) x1 e
		template<class X>
		fms::dual<X> operator()(fms::dual<X> x)
		{
			return fms::dual<X>(f(x._0), df(x._0)*x._1);
		}
	};

	// derivative of f at x, f(x + e) = f(x) + f'(x) e
	// f:dual<X> -> dual<X>
	template<class F>
	class D {
		F f;
	public:
		D(F f)
			: f{ f }
		{ }
		template<class X>
		X operator()(X x)
		{
			return f(dual<X>(x, X(1)))._1;
		}
	};

} // namespace fms

template<std::floating_point X>
inline constexpr fms::dual<X> operator+(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 + y._0, x._1 + y._1);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator+(const X& x, const fms::dual<X>& y)
{
	return fms::dual(x + y._0, y._1);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator+(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0 + y, x._1);
}

// negative of a dual number
template<std::floating_point X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x)
{
	return fms::dual(-x._0, -x._1);
}

template<std::floating_point X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 - y._0, x._1 - y._1);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator-(const X& x, const fms::dual<X>& y)
{
	return fms::dual(x - y._0, -y._1);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator-(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0 - y, x._1);
}

template<std::floating_point X>
inline constexpr fms::dual<X> operator*(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return fms::dual(x._0 * y._0, x._0 * y._1 + x._1 * y._0);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator*(const X& x, const fms::dual<X>& y)
{
	return fms::dual(x * y._0, x * y._1);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator*(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0 * y, x._1 * y);
}

// inverse of a dual number, x*inv(x) == dual(1, 0) == inv(x)*x
template<std::floating_point X>
inline constexpr fms::dual<X> inv(const fms::dual<X>& x)
{
	return fms::dual(1 / x._0, -x._1 / (x._0 * x._0));
}

template<std::floating_point X>
inline constexpr fms::dual<X> operator/(const fms::dual<X>& x, const fms::dual<X>& y)
{
	return x * inv(y);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator/(const X& x, const fms::dual<X>& y)
{
	return x * inv(y);
}
template<std::floating_point X>
inline constexpr fms::dual<X> operator/(const fms::dual<X>& x, const X& y)
{
	return fms::dual(x._0/y, x._1/y);
}

// Promote to dual functions.
template<class X> inline auto _log =
	fms::_([](X x) { return log(x); }, [](X x) { return 1 / x; });
template<class X> inline auto _exp =
	fms::_([](X x) { return exp(x); }, [](X x) { return exp(x); });

