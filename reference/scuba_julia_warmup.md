# Warm up the Julia engine

Runs each accessor once on tiny in-memory and BPCells-backed Seurat
objects so the first call on a real object does not pay for compilation
(about 0.6 s with JuliaCall, even with precompiled code). Called
automatically at engine start-up unless
`options(SCUBA.julia_warmup = FALSE)`.

## Usage

``` r
scuba_julia_warmup(objects = list())
```

## Arguments

- objects:

  Optional list of Seurat, SingleCellExperiment or `scuba_source`
  objects to warm.

## Value

With `objects`, a named logical vector (one entry per object: `TRUE`
warmed, `FALSE` failed with a warning naming the object, `NA` not a
Seurat/SingleCellExperiment object), invisibly; otherwise `TRUE`
invisibly (`FALSE` if the engine is not available).

## Details

Pass the objects an application will serve as `objects` to also wrap
each of them and run one expression fetch, so the first user request
does not pay for reading the object's cell names and compiling for its
exact layout (about 0.1-1 s per object depending on size). A Shiny app
should call this once, after loading its objects and before serving
sessions.
