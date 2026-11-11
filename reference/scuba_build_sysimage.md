# Build a Julia system image with SCUBA.jl compiled in

Compiles SCUBA.jl and its dependencies into a Julia system image, which
cuts the Julia start-up of each R process from about 5 s to about 3 s.
RCall itself cannot be compiled into the image on Julia 1.13
(PackageCompiler fails with an LLVM calling-convention error), so it and
the RCall extension load from their own precompile caches. Building
takes several minutes once and runs in a separate `julia` process
(equivalent to
`julia <SCUBA.jl>/tools/build_sysimage.jl --no-rcall --no-pythoncall`
from a terminal). Afterwards set `options(SCUBA.sysimage = path)` (or
the `SCUBA_SYSIMAGE` environment variable) before the first SCUBA call
in a session.

## Usage

``` r
scuba_build_sysimage(path = NULL, julia_path = NULL)
```

## Arguments

- path:

  Output path for the image (default:
  `~/.julia/scuba/SCUBA.<os>-<arch>.<ext>`).

- julia_path:

  Path to the Julia binary (default: the one JuliaCall finds).

## Value

The image path, invisibly.
