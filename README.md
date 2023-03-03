# fms_dual

Compute derivatives to machine precision.

## Build

Requires gcc version 11.1.0 or later.

```
$ git clone https://github.com/keithalewis/fms_dual
$ cd fms_dual
$ make
```

## Run

```
$ make test
```

## Dual Numbers

Suppose $\epsilon \not= 0$ but $\epsilon^2 = 0$.
By Taylor's theorem we have $f(x + \epsilon) = f(x) + f'(x) \epsilon$.

__Exercise__. _If_ $f(x) = x^2$ _then_ $f'(x) = 2x$.

_Hint_: $(x + \epsilon^2) = x^2 + 2x\epsilon + \epsilon^2$.

Of course $\epsilon$ can't be a real number, but the $2\times 2$
matrix $\epsilon = [0, 1; 0, 0]$ is not zero.

__Exercise__. _Show_ $\epsilon^2 = 0$.

_Hint_: Where $0$ is the $2\times 2$ zero matrix.

The `template<class X> struct fms::dual { X _0, _1; }` corresponds
to the dual number $x = x_0 + x_1 \epsilon$.
Unary negation and `inv`erse are definded in the global namespace
as are all the arithmetic operations and their overloads
for scalar operands. This allows dual numbers to be used
just like [`std::complex`](https://en.cppreference.com/w/cpp/numeric/complex).

A function `X f(X)` can be promoted to a function `dual<X> _f(dual<X>)`
if its derivative is known since 
$f(x_0 + x_1 \epsilon) = f(x_0) + f'(x_0)x_1\epsilon$.
Use `auto _f = _(f, df)` where `df` is the derivative.

Dual number definitions are provide for some standard functions, e.g.,
```
template<class X> _log = _([](X x) { return log(x); }, [](X x) { return -1/x; });
```

__Exercise__. _Implement_ `_N(dual<X>)` _for the standard normal distribution_.

_Hint_: $N(x) = (1 + \operatorname{erf}(x / \sqrt{2}) / 2)$ is the standard normal
cumulative distribution in terms of the C standard libary
[`erf`](https://en.cppreference.com/w/c/numeric/math/erf) and
$N'(x) = \exp(-x^2/2)/\sqrt{2\pi}$.

Dual numbers have far less machinery than Automatic Differentiation.

## Unfiled
