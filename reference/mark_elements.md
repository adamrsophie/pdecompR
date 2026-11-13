# Mark elements for refinement

`mark_elements()` implements Doerfler (bulk) marking, which selects a
smallest set of elements whose squared indicators sum to at least
`theta` times the total, and maximum marking, which selects elements
with `eta_T >= theta * max(eta)`.

## Usage

``` r
mark_elements(eta, theta = 0.5, strategy = c("doerfler", "maximum"))
```

## Arguments

- eta:

  Element error indicators.

- theta:

  Marking parameter in (0, 1\].

- strategy:

  `"doerfler"` or `"maximum"`.

## Value

Integer indices of marked elements.

## References

Doerfler, W. (1996). A convergent adaptive algorithm for Poisson's
equation. *SIAM Journal on Numerical Analysis*, 33, 1106-1124.

## Examples

``` r
mark_elements(c(0.1, 0.5, 0.2, 0.05), theta = 0.5)
#> [1] 2
```
