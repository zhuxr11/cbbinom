# cbbinom (development version)

## Bug fixes

* Fix a bug in `cpp_qcbbinom()` that may change input `p`, `tol` and `max_iter`.

## Enhancements

* Add `const` keywords to `UnirootEqn` input of `cpp_uniroot()`, so that it may take a constant pointer to a `UnirootEqn` object, a pointer to a constant `UnirootEqn` object, or both.

* Use `BH` package to compute generalized hypergeometric functions and their numerical derivatives; also use `NULL` for self-determined tolerance.

* Add precision level in generalized hypergeometric functions for flexible trade-off of computational accuracy and time.

* Use `hypergeo2` package to compute generalized hypergeometric function at high precision.


# cbbinom 0.1.0

* Initial CRAN submission.
