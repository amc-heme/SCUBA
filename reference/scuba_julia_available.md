# Is the Julia engine available?

Initialises Julia through JuliaCall on first use (which may take a while
the first time, as SCUBA.jl is installed and precompiled) and caches the
result for the session.

## Usage

``` r
scuba_julia_available(quiet = FALSE)
```

## Arguments

- quiet:

  Suppress the message shown while Julia starts.
