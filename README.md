# fms_dual

Compute derivatives to machine precision.

Suppose $\epsilon \not= 0$ but $\epsilon^2 = 0$.
By Taylor's theorem we have $f(x + \epsilon) = f(x) + f'(x) \epsilon$.

__Exercise__. _If_ $f(x) = x^2$ _then_ $f'(x) = 2x$.

_Hint_: $(x + \epsilon^2) = x^2 + 2x\epsilon + \epsilon^2$.

Of course $\epsilon$ can't be a real number, but the $2\times 2$
matrix $\epsilon = \begin{bmatrix}0 & 1 \\ 0 & 0 \end{bmatrix}$
is not zero.

__Exercise__. _Show $\epsilon^2 = \begin{bmatrix}0 & 0 \\ 0 & 0 \end{bmatrix}$_.

```
norm(a + be) = max(sqrt(-b sqrt(4 a^2 + b^2) + 2 a^2 + b^2)/sqrt(2),
                   sqrt(b sqrt(4 a^2 + b^2) + 2 a^2 + b^2)/sqrt(2))
```