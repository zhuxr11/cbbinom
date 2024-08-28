# cbbinom (development version)

## Bug fixes

* Fix a bug in `gen_hypergeo()` that may change input `U` and `L`.

* Fix a bug in `cpp_qcbbinom()` that may change input `p`, `tol` and `max_iter`.

## Enhancements

* Add `const` keywords to `UnirootEqn` input of `cpp_uniroot()`, so that it may take a constant pointer to a `UnirootEqn` object, a pointer to a constant `UnirootEqn` object, or both.


# cbbinom 0.1.0

* Initial CRAN submission.
