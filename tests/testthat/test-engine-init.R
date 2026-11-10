# Julia engine start-up checks that need no running Julia, plus value passing when it runs.

test_that("a SCUBA.jl without ScubaR is rejected so auto falls back to R (review #1)", {
  expect_error(SCUBA:::.scuba_require_scubar(FALSE, "@scuba-r"), "too old.*@scuba-r")
  expect_invisible(SCUBA:::.scuba_require_scubar(TRUE, "@scuba-r"))
  expect_false(exists(".scuba_julia_glue", asNamespace("SCUBA"), inherits = FALSE))
})

test_that("paths reach Julia as values, backslashes and quotes intact (review #7)", {
  skip_if_not(scuba_julia_available(quiet = TRUE), "Julia engine not available")
  windows_path <- "C:\\Users\\a \"b\"\\SCUBA.jl"
  JuliaCall::julia_assign("scuba_r_path_test", windows_path)
  expect_identical(JuliaCall::julia_eval("scuba_r_path_test"), windows_path)
  engine_source <- deparse(SCUBA:::.scuba_julia_init)
  expect_false(any(grepl("Pkg.develop(path = \"%s\"", engine_source, fixed = TRUE)))
  expect_false(any(grepl("Pkg.add(url = \"%s\"", engine_source, fixed = TRUE)))
})
