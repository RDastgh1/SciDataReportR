# Negated "in" Operator

This custom operator returns TRUE if the element on the left-hand side
of the operator is not found in the vector on the right-hand side when
using the %in% operator. Otherwise, it returns FALSE.

## Usage

``` r
x %!in% y
```

## Arguments

- x:

  The element to be tested for absence in the vector.

- y:

  The vector in which to search for the element.

## Value

TRUE if the element is not found in the vector, otherwise FALSE.
