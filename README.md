# Scale-Optimized Geometric Acceleration of Polygonal Approximations to Pi

**Author:** C. Wayne Baker  
**Original date:** May 9, 2026  
**Status:** Revised mathematical note / reproducibility archive

## Read the revised paper

- [`paper/scale_optimized_polygonal_pi_acceleration.pdf`](paper/scale_optimized_polygonal_pi_acceleration.pdf)
- [`paper/scale_optimized_polygonal_pi_acceleration.tex`](paper/scale_optimized_polygonal_pi_acceleration.tex)

## Main result

For

\[
P_N=N\sin\!\left(\frac{\pi}{N}\right),
\]

geometric scale nodes \(b^mN\), and Lagrange weights at
\(x_m=b^{-2m}\), the operator

\[
\widehat P^{[r]}_{N,b}
=\sum_{m=0}^{r-1}w_m^{(r,b)}P_{b^mN}
\]

satisfies

\[
\widehat P^{[r]}_{N,b}-\pi=O(N^{-2r}).
\]

The geometric-node weights obey the uniform stability bound

\[
\sum_m|w_m^{(r,b)}|
\le
S_b=
\prod_{\ell=1}^{\infty}
\frac{1+b^{-2\ell}}{1-b^{-2\ell}}.
\]

## Constructibility boundary

The dyadic branch is exactly constructible by nested square roots. Exact
tripling requires a cubic inversion and is not generally
straightedge-and-compass constructible. The revised paper therefore treats
tripling branches as scale-architecture tests unless a separate approximate
construction meets the stated perturbation bound.

The short April note is preserved in `docs/`, but it is not the authoritative
paper because it described `6 -> 18 -> 54 -> ...` too strongly as exact
constructible refinement.

## Reproduce

```bash
python -m pip install mpmath
python scripts/verify_scale_optimized_polygonal_pi.py --precision 1300 --output-dir output
```

Expected final line:

```text
verification_passed: True
```

## Rights

Copyright © 2026 C. Wayne Baker. All rights reserved unless otherwise stated.
