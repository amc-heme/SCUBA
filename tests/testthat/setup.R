# Engine selection for the test suite.
#   SCUBA_TEST_ENGINE = "r" | "julia" | "auto" (default "auto": Julia if it can be started)
#   SCUBA_JL_PATH     = development checkout of SCUBA.jl (defaults to ../../../SCUBA.jl if present)
if (!nzchar(Sys.getenv("SCUBA_JL_PATH"))) {
  cand <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", "SCUBA.jl"), mustWork = FALSE)
  if (dir.exists(cand)) Sys.setenv(SCUBA_JL_PATH = cand)
}
options(SCUBA.engine = Sys.getenv("SCUBA_TEST_ENGINE", unset = "auto"))
message("SCUBA tests: engine = ", tryCatch(scuba_engine(), error = function(e) "r (julia unavailable)"))
