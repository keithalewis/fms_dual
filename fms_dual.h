// fms_dual.h - dual numbers  x0 + x1 e  where e != 0 but e^2 = 0
#pragma once
#include <concepts>

namespace fms {

	// x = _0 + _1 e
	template<std::floating_point X = double>
	struct dual {
		X _0, _1;

		constexpr dual(const X& _0 = 0, const X& _1 = 0)
			: _0(_0), _1(_1)
		{ }
		dual(const dual&) = default;
		dual& operator=(const dual&) = default;
		~dual()
		{ }

		bool operator==(const dual& x) const
		{
			return _0 == x._0 and _1 == x._1;
		}

#ifdef _DEBUG
		static int test()
		{
			{
				dual d;
				assert(d == 0);
				dual d2{ d };
				assert(d2 == d);
				d = d2;
				assert(!(d != d2));
				assert(d == -d);
			}
			{
				dual d(1);
				assert(d == 1);
			}

			return 0;
		}
#endif // _DEBUG

	};

} // namespace fms

// negative of a dual number
template<std::floating_point X>
constexpr fms::dual<X> operator-(const fms::dual<X>& x)
{
	return fms::dual(-x._0, -x._1);
}

// inverse of a dual number, a*inv(a) == dual(1, 0) == inv(a)*a
template<std::floating_point X>
constexpr fms::dual<X> inv(const fms::dual<X>& a)
{
	return fms::dual(1/a._0, -a._0/(a._1 * a._1));
}
/*
template<std::floating_point X>
constexpr fms::dual<X> operator+(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return fms::dual(a._0 + b._0, a._1 + b._1);
}
template<std::floating_point X>
constexpr fms::dual<X> operator-(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return a + -b;
}
template<std::floating_point X>
constexpr fms::dual<X> operator*(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return fms::dual(a._0 * b._0, a._0 * b._1 + a._1 * b._0);
}
template<std::floating_point X>
constexpr fms::dual<X> operator/(const fms::dual<X>& a, const fms::dual<X>& b)
{
	return a * inv(b);
}

// derivative of f at x, f(x + e) = f(x) + f'(x) e
template<class F, std::floating_point X>
constexpr X dF(const F& f, const fms::dual<X>& x)
{
	return f(x)._1;
}
*/
