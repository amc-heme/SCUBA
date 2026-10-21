# The Linux LD_LIBRARY_PATH check runs before Julia starts; no Julia needed for these tests.
chk <- SCUBA:::.scuba_ld_library_path_issue

test_that("system library directories on LD_LIBRARY_PATH are detected on Linux", {
  expect_match(chk("/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu", "Linux"), "x86_64-linux-gnu")
  expect_match(chk("/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-17/lib/server", "Linux"), "R_LD_LIBRARY_PATH")
  expect_match(chk("/usr/lib64", "Linux"), "/usr/lib64")
  expect_match(chk("/lib/aarch64-linux-gnu/", "Linux"), "aarch64")
  expect_match(chk("/usr/local/lib", "Linux"), "/usr/local/lib")
})

test_that("safe paths pass", {
  expect_null(chk("", "Linux"))
  expect_null(chk("/usr/lib/R/lib", "Linux"))
  expect_null(chk("/usr/local/lib/R/lib:/opt/intel/mkl/lib", "Linux"))
  expect_null(chk("/usr/lib/jvm/java-17/lib/server", "Linux"))
  # Julia's own directories first take precedence over later system entries
  expect_null(chk("/home/u/.julia/juliaup/julia-1.13.2/lib/julia:/home/u/.julia/juliaup/julia-1.13.2/lib:/usr/lib/x86_64-linux-gnu", "Linux"))
  # macOS does not use LD_LIBRARY_PATH this way
  expect_null(chk("/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu", "Darwin"))
})
