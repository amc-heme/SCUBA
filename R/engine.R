# Engine selection: SCUBA 2.x can run its accessors either through the Julia package SCUBA.jl
# (one implementation of the resolution rules shared with Python and Julia) or through the
# legacy pure-R methods. Nothing Julia-related happens at load time.
#
#   options(SCUBA.engine = "auto")    # default: Julia when available, else R
#   options(SCUBA.engine = "julia")   # require Julia (errors if it cannot be set up)
#   options(SCUBA.engine = "r")       # legacy pure-R implementation
#
# Julia-side package location: option `SCUBA.jl_path` / env `SCUBA_JL_PATH` points at a
# development checkout of SCUBA.jl; otherwise the registered package is installed into a
# dedicated shared Julia environment (`~/.julia/environments/scuba-r`).

.scuba_env <- new.env(parent = emptyenv())
.scuba_env$julia_ready <- NULL

#' Which engine SCUBA is using
#'
#' @return `"julia"` or `"r"`.
#' @export
scuba_engine <- function() {
  opt <- getOption("SCUBA.engine", "auto")
  if (!opt %in% c("auto", "julia", "r")) {
    stop("options(SCUBA.engine) must be 'auto', 'julia' or 'r'", call. = FALSE)
  }
  if (opt == "r") return("r")
  if (opt == "julia") {
    if (!scuba_julia_available()) {
      stop("options(SCUBA.engine = 'julia') but the Julia engine could not be initialised. ",
           "Install Julia (https://julialang.org/downloads) and the R package JuliaCall, or use ",
           "options(SCUBA.engine = 'r').", call. = FALSE)
    }
    return("julia")
  }
  if (isTRUE(scuba_julia_available())) "julia" else "r"
}

#' Is the Julia engine available?
#'
#' Initialises Julia through JuliaCall on first use (which may take a while the first time,
#' as SCUBA.jl is installed and precompiled) and caches the result for the session.
#' @param quiet Suppress the message shown while Julia starts.
#' @export
scuba_julia_available <- function(quiet = FALSE) {
  if (!is.null(.scuba_env$julia_ready)) return(.scuba_env$julia_ready)
  ok <- tryCatch({
    .scuba_julia_init(quiet = quiet)
    TRUE
  }, error = function(e) {
    if (!quiet && !identical(getOption("SCUBA.engine", "auto"), "r")) {
      message("SCUBA: Julia engine unavailable (", conditionMessage(e), "); using the R engine.")
    }
    FALSE
  })
  .scuba_env$julia_ready <- ok
  ok
}

#' Run an expression with a given SCUBA engine
#' @param engine `"julia"` or `"r"`.
#' @param expr Expression to evaluate.
#' @export
with_engine <- function(engine, expr) {
  old <- options(SCUBA.engine = engine)
  on.exit(options(old), add = TRUE)
  force(expr)
}

# On Debian/Ubuntu the `R` launcher script (etc/ldpaths) exports
# LD_LIBRARY_PATH=$R_HOME/lib:/usr/lib/x86_64-linux-gnu. With the system library directory on
# that path, Julia embedded in R resolves libstdc++, libgfortran, OpenBLAS and friends to the
# system copies instead of its bundled ones, and the R process segfaults with no R-level error
# (observed on Ubuntu with R 4.5 / Julia 1.13). The loader reads LD_LIBRARY_PATH once at process
# start, so this cannot be repaired from inside R; it must be set before R starts, e.g.
# R_LD_LIBRARY_PATH=/usr/lib/R/lib in the shiny-server service environment.
# Returns NULL when the path looks safe, otherwise an explanation.
.scuba_ld_library_path_issue <- function(path = Sys.getenv("LD_LIBRARY_PATH"),
                                         sysname = Sys.info()[["sysname"]]) {
  if (!identical(sysname, "Linux") || !nzchar(path)) return(NULL)
  entries <- strsplit(path, ":", fixed = TRUE)[[1]]
  entries <- sub("/+$", "", entries[nzchar(entries)])
  system_dir <- grepl("^/(usr/)?lib(32|64|x32)?(/[^/]+-linux-gnu[^/]*)?$", entries) |
    grepl("^/usr/local/lib(32|64)?$", entries)
  julia_dir <- grepl("/lib/julia$", entries)
  bad <- which(system_dir)
  if (!length(bad)) return(NULL)
  # Julia's own library directories listed before the system ones take precedence: fine
  if (any(julia_dir) && min(which(julia_dir)) < min(bad)) return(NULL)
  paste0(
    "LD_LIBRARY_PATH contains the system library directory ", entries[bad[1]],
    " ahead of Julia's libraries, which makes Julia embedded in R load system copies of its ",
    "bundled libraries and crash. Start R with R_LD_LIBRARY_PATH=", file.path(R.home(), "lib"),
    " (set it in the environment that launches R, e.g. the shiny-server service), or prepend ",
    "<julia>/lib/julia:<julia>/lib to LD_LIBRARY_PATH. Set options(SCUBA.ignore_ld_library_path = TRUE) ",
    "to try anyway."
  )
}

.scuba_julia_init <- function(quiet = FALSE) {
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("package 'JuliaCall' is not installed")
  }
  if (!isTRUE(getOption("SCUBA.ignore_ld_library_path", FALSE))) {
    issue <- .scuba_ld_library_path_issue()
    if (!is.null(issue)) stop(issue, call. = FALSE)
  }
  img <- .scuba_sysimage()
  if (!quiet) message("SCUBA: starting the Julia engine", if (is.null(img)) "" else paste0(" (sysimage ", img, ")"), " ...")
  dev_path <- getOption("SCUBA.jl_path", Sys.getenv("SCUBA_JL_PATH", unset = ""))
  project <- getOption("SCUBA.jl_project",
                       Sys.getenv("SCUBA_JL_PROJECT", unset = "scuba-r"))
  # a bare name is a shared environment: ~/.julia/environments/<name>
  shared <- !grepl("[/\\\\]", project)
  if (!is.null(img) && !nzchar(Sys.getenv("JULIA_PROJECT"))) {
    # With a system image SCUBA is loaded before anything runs, so its RCall extension must be
    # locatable from the default environment when JuliaCall loads RCall during start-up.
    Sys.setenv(JULIA_PROJECT = if (shared) paste0("@", project) else normalizePath(project))
  }
  if (is.null(img)) JuliaCall::julia_setup(installJulia = FALSE, verbose = FALSE)
  else JuliaCall::julia_setup(installJulia = FALSE, verbose = FALSE, sysimage_path = img)
  JuliaCall::julia_command("import Pkg")
  active <- JuliaCall::julia_eval("Base.active_project()")
  want <- if (shared) file.path(path.expand("~"), ".julia", "environments", project, "Project.toml") else file.path(project, "Project.toml")
  if (!identical(normalizePath(active, mustWork = FALSE), normalizePath(want, mustWork = FALSE))) {
    if (shared) {
      JuliaCall::julia_command(sprintf('Pkg.activate("%s"; shared = true, io = devnull)', project))
    } else {
      JuliaCall::julia_command(sprintf('Pkg.activate("%s"; io = devnull)', project))
    }
  }
  have <- JuliaCall::julia_eval('Base.find_package("SCUBA") !== nothing')
  if (!isTRUE(have)) {
    if (nzchar(dev_path)) {
      JuliaCall::julia_command(sprintf('Pkg.develop(path = "%s"; io = devnull)', normalizePath(dev_path)))
    } else {
      JuliaCall::julia_command('Pkg.add("SCUBA"; io = devnull)')
    }
  }
  for (pkg in c("RCall", "DataFrames", "CategoricalArrays")) {
    if (!isTRUE(JuliaCall::julia_eval(sprintf('Base.find_package("%s") !== nothing', pkg)))) {
      JuliaCall::julia_command(sprintf('Pkg.add("%s"; io = devnull)', pkg))
    }
  }
  loaded <- tryCatch({ JuliaCall::julia_library("SCUBA"); TRUE }, error = function(e) FALSE)
  if (!loaded) {
    # a development checkout may have gained dependencies since the environment was resolved
    # (resolving takes ~10 s, so it only runs when loading failed)
    JuliaCall::julia_command("Pkg.resolve(io = devnull); Pkg.instantiate(io = devnull)")
    JuliaCall::julia_library("SCUBA")
  }
  JuliaCall::julia_library("RCall")
  if (!isTRUE(JuliaCall::julia_eval("Base.get_extension(SCUBA, :SCUBARCallExt) !== nothing"))) {
    JuliaCall::julia_command("Base.retry_load_extensions()")
  }
  # The glue module ships (precompiled) with SCUBA.jl's RCall extension; older SCUBA.jl versions
  # get the copy embedded below.
  in_ext <- JuliaCall::julia_eval("isdefined(Base.get_extension(SCUBA, :SCUBARCallExt), :ScubaR)")
  if (isTRUE(in_ext)) {
    JuliaCall::julia_command("isdefined(Main, :ScubaR) || (global ScubaR = Base.get_extension(SCUBA, :SCUBARCallExt).ScubaR)")
  } else {
    JuliaCall::julia_command(.scuba_julia_glue)
  }
  if (isTRUE(getOption("SCUBA.julia_warmup", TRUE))) scuba_julia_warmup()
  invisible(TRUE)
}

#' Warm up the Julia engine
#'
#' Runs each accessor once on a tiny in-memory Seurat object so the first call on a real object
#' does not pay for compilation (about 0.6 s with JuliaCall, even with precompiled code). Called
#' automatically at engine start-up unless `options(SCUBA.julia_warmup = FALSE)`; a Shiny app can
#' call it from `global.R` so sessions start warm.
#'
#' @returns `TRUE` invisibly.
#' @export
scuba_julia_warmup <- function() {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) return(invisible(FALSE))
  m <- Matrix::sparseMatrix(i = c(1, 2, 3, 1, 4, 2, 5), j = c(1, 1, 2, 3, 3, 4, 4), x = c(3, 1, 2, 5, 1, 4, 2),
                            dims = c(5, 4), dimnames = list(paste0("g", 1:5), paste0("c", 1:4)))
  md <- data.frame(group = factor(c("a", "b", "a", "b")), score = c(0.5, 1.5, 2.5, 3.5), row.names = colnames(m))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = m, meta.data = md))
  obj <- SeuratObject::SetAssayData(obj, layer = "data", new.data = m)
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(as.numeric(1:8), 4, 2, dimnames = list(colnames(m), c("UMAP_1", "UMAP_2"))), key = "UMAP_", assay = "RNA")
  objs <- list(obj)
  if (requireNamespace("BPCells", quietly = TRUE)) {
    # on-disk BPCells layers in both storage orders (the layout of cellDIVER objects)
    tmp <- tempfile("scuba_warmup_bp")
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    bp <- tryCatch({
      col <- suppressMessages(BPCells::write_matrix_dir(methods::as(m, "IterableMatrix"), file.path(tmp, "col")))
      row <- suppressMessages(BPCells::transpose_storage_order(col, outdir = file.path(tmp, "row"), tmpdir = tmp))
      o <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = row, meta.data = md))
      o[["RNA"]]$data <- col
      o[["umap"]] <- obj[["umap"]]
      o
    }, error = function(e) NULL)
    if (!is.null(bp)) objs <- c(objs, list(bp))
  }
  for (obj in objs) suppressWarnings(with_engine("julia", {
    .jl_call("get_all_cells", obj)
    .jl_call("fetch_data", obj, "rna_g1", NULL, NULL)
    .jl_call("fetch_data", obj, c("rna_g1", "group", "UMAP_1"), NULL, c("c2", "c1"))
    .jl_call("fetch_feature", obj, c("g1", "g2"), NULL, NULL)
    .jl_call("fetch_metadata", obj, NULL, NULL, TRUE)
    .jl_call("fetch_metadata", obj, c("group", "score"), c("c1", "c3"), FALSE)
    .jl_call("fetch_reduction", obj, "umap", NULL, c(1L, 2L))
    .jl_call("fetch_cells", obj, "group", "a")
    .jl_call("unique_values", obj, "group")
    .jl_call("meta_varnames", obj); .jl_call("all_keys", obj); .jl_call("features_in_assay", obj, NULL)
    .jl_call("default_layer", obj); .jl_call("default_reduction_name", obj)
  }))
  JuliaCall::julia_command("ScubaR.clear_cache()")
  invisible(TRUE)
}

# Julia glue evaluated once per session. Every entry point wraps the R object as a SCUBA
# source, runs the accessor collecting warnings, and returns a plain list for R.
.scuba_julia_glue <- '
module ScubaR
using SCUBA, RCall, DataFrames, CategoricalArrays

const S = SCUBA

# Cache of wrapped sources keyed by the R object pointer. Holding the RObject keeps the SEXP
# alive (RCall preserves it), so an address cannot be reused by a different object while cached.
# R copies on modify, so a modified object arrives with a new address and gets a fresh source.
const SRC_CACHE = Dict{Ptr{Cvoid}, Tuple{Any, Any}}()
const SRC_CACHE_MAX = 4
function getsrc(obj)
    key = Ptr{Cvoid}(obj.p)
    hit = get(SRC_CACHE, key, nothing)
    hit === nothing || return hit[2]
    if length(SRC_CACHE) >= SRC_CACHE_MAX
        empty!(SRC_CACHE)
    end
    src = S.rsource(obj)
    SRC_CACHE[key] = (obj, src)
    return src
end
clear_cache() = (empty!(SRC_CACHE); nothing)

function run(f, obj, args...; kwargs...)
    src = getsrc(obj)
    result, ws = S.with_collected_warnings(() -> f(src, args...; kwargs...))
    return result, String[w.message for w in ws]
end

asvec(x::AbstractString) = String[x]
asvec(x) = x

# Hand a table to R column by column so factors, NA and types survive the JuliaCall conversion.
function totable(df::DataFrame)
    cc = S.cell_column(df)
    names_ = String[n for n in names(df) if n != cc]
    columns = Any[]
    cat_idx = Int[]; cat_codes = Any[]; cat_levels = Any[]; cat_ordered = Bool[]
    for (k, n) in enumerate(names_)
        v = df[!, n]
        if v isa CategoricalVector
            push!(cat_idx, k)
            push!(cat_codes, Int[ismissing(x) ? -1 : levelcode(x) for x in v])
            push!(cat_levels, string.(levels(v)))
            push!(cat_ordered, isordered(v))
            push!(columns, Int[0])
        else
            push!(columns, v)
        end
    end
    (names = names_, columns = columns, cat_idx = cat_idx, cat_codes = cat_codes, cat_levels = cat_levels,
     cat_ordered = cat_ordered, cells = String.(S.cells(df)))
end

function fetch_data(obj, vars, layer, cells)
    df, w = run(S.fetch_data, obj, asvec(vars); layer = layer, cells = cells)
    (totable(df)..., warnings = w)
end
function fetch_feature(obj, features, layer, cells)
    df, w = run(S.fetch_feature, obj, asvec(features); layer = layer, cells = cells)
    (totable(df)..., warnings = w)
end
function fetch_metadata(obj, vars, cells, full_table)
    df, w = if full_table
        run(S.fetch_metadata, obj; cells = cells)
    else
        run(S.fetch_metadata, obj, asvec(vars); cells = cells)
    end
    (totable(df)..., warnings = w)
end
function fetch_reduction(obj, reduction, cells, dims)
    df, w = run(S.fetch_reduction, obj, reduction; cells = cells, dims = dims)
    (totable(df)..., warnings = w)
end
fetch_cells(obj, meta_var, meta_levels) = run(S.fetch_cells, obj, meta_var, asvec(meta_levels))
get_all_cells(obj) = run(S.all_cells, obj)
meta_varnames(obj) = run(S.meta_varnames, obj)
function all_keys(obj)
    # R SCUBA returns names = assay / reduction names and values = keys; mirror that for Seurat
    src = getsrc(obj)
    ks = S.all_keys(src)
    names_ = collect(keys(ks)); vals = collect(values(ks))
    if hasproperty(src, :assay_alias)
        inv = Dict(v => k for (k, v) in src.assay_alias)
        merge!(inv, Dict(v => k for (k, v) in src.red_alias))
        names_ = [get(inv, k, k) for k in names_]
    end
    return names_, vals, String[]
end
features_in_assay(obj, assay) = run(S.features_in_modality, obj, assay)
reduction_dimnames(obj, reduction, dims) = run(S.reduction_dimnames, obj, reduction, dims)
default_layer(obj) = run(S.default_layer, obj)
default_reduction(obj) = run(S.default_reduction, obj)
function default_reduction_name(obj)
    # R SCUBA returns the reduction *name* for Seurat ("umap"), the key for others
    src = getsrc(obj)
    key = S.default_reduction(src)
    name = key
    if hasproperty(src, :reductions) && haskey(src.reductions, key) && hasproperty(src.reductions[key], :name)
        name = src.reductions[key].name
    end
    return name, String[]
end
unique_values(obj, var) = run(S.unique_values, obj, var)
n_unique(obj, var) = run(S.n_unique, obj, var)
end
'

# ---------------------------------------------------------------------------
# R-side shims

# Routing. The Julia engine only wins where Julia reads matrices natively, i.e. expression
# fetches. Metadata, reduction and cell queries on live R objects are served by the R methods
# directly: through Julia they would round-trip R -> Julia -> R (1.8 s for a 1M-cell metadata
# table vs 2 ms in R, see bench/report/summary.md). `fetch_data()` requests that contain only
# metadata / reduction vars stay in R as well; mixed requests go to Julia, which fetches the
# features and appends the metadata columns.
.scuba_julia_object <- function(object) {
  (inherits(object, "Seurat") || inherits(object, "SingleCellExperiment")) &&
    identical(scuba_engine(), "julia")
}
.scuba_use_julia <- function(object) .scuba_julia_object(object)

.r_meta_names <- function(object) {
  if (inherits(object, "Seurat")) colnames(object[[]])
  else if (inherits(object, "SingleCellExperiment")) colnames(SummarizedExperiment::colData(object))
  else character()
}
.r_reduction_keys <- function(object) {
  if (inherits(object, "Seurat")) {
    unlist(lapply(SeuratObject::Reductions(object), function(r) SeuratObject::Key(object[[r]])))
  } else if (inherits(object, "SingleCellExperiment")) {
    paste0(SingleCellExperiment::reducedDimNames(object), "_")
  } else character()
}

#' @keywords internal
.scuba_route_fetch_data <- function(object, vars) {
  if (!.scuba_julia_object(object)) return(FALSE)
  vars <- as.character(vars)
  meta <- as.character(.r_meta_names(object))
  reds <- as.character(.r_reduction_keys(object))
  is_meta <- vars %in% meta | sub("^obs_", "", vars) %in% meta
  is_red <- if (length(reds)) {
    vapply(vars, function(v) any(startsWith(v, reds)) && grepl("_[0-9]+$", v), logical(1))
  } else rep(FALSE, length(vars))
  # metadata- and reduction-only requests are faster on the R methods
  !all(is_meta | is_red)
}

# Legacy fetch_feature resolves a *bare* name that is both a metadata column and a feature of the
# default assay through fetch_data, i.e. metadata first (unless the assay is BPCells-backed). The
# Julia engine never searches metadata in fetch_feature, so such requests stay on the R method.
.scuba_feature_collides <- function(object, features) {
  features <- as.character(features)
  meta <- .r_meta_names(object)
  if (!length(meta)) return(FALSE)
  keys <- tryCatch(as.character(all_keys(object)), error = function(e) character())
  keys <- keys[nzchar(keys)]
  keyed <- if (length(keys)) Reduce(`|`, lapply(keys, function(k) startsWith(features, k))) else rep(FALSE, length(features))
  any(!keyed & features %in% meta)
}

.jl_call <- function(fn, ...) {
  res <- tryCatch(
    JuliaCall::julia_call(paste0("ScubaR.", fn), ...),
    error = function(e) {
      msg <- conditionMessage(e)
      msg <- sub("^Error happens in Julia.\\s*", "", msg)
      msg <- sub("(?s)\\s*Stacktrace:.*$", "", msg, perl = TRUE)
      msg <- sub("^[A-Za-z]+Error: ", "", msg)
      stop(msg, call. = FALSE)
    }
  )
  res
}

.jl_emit_warnings <- function(w) {
  for (m in w) warning(m, call. = FALSE)
  invisible(NULL)
}

.jl_table <- function(res) {
  cols <- res$columns
  if (length(res$cat_idx) > 0) {
    for (k in seq_along(res$cat_idx)) {
      codes <- res$cat_codes[[k]]
      lv <- res$cat_levels[[k]]
      codes[codes < 1] <- NA
      cols[[res$cat_idx[k]]] <- factor(lv[codes], levels = lv, ordered = isTRUE(res$cat_ordered[k]))
    }
  }
  if (length(cols) == 0) {
    df <- data.frame(row.names = res$cells)
  } else {
    names(cols) <- res$names
    df <- as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE, optional = TRUE)
    names(df) <- res$names
    rownames(df) <- res$cells
  }
  .jl_emit_warnings(res$warnings)
  df
}

.jl_value <- function(res) {
  .jl_emit_warnings(res[[2]])
  res[[1]]
}

.jl_fetch_data <- function(object, vars, layer = NULL, cells = NULL, slot = NULL) {
  layer <- slot %||% layer
  .jl_table(.jl_call("fetch_data", object, as.character(vars), layer, cells))
}
.jl_fetch_feature <- function(object, features, layer = NULL, cells = NULL) {
  .jl_table(.jl_call("fetch_feature", object, as.character(features), layer, cells))
}
.jl_fetch_metadata <- function(object, vars = NULL, cells = NULL, full_table = FALSE,
                               return_class = "dataframe") {
  if (!return_class %in% c("dataframe", "vector")) {
    stop("return_class must be 'dataframe' or 'vector'", call. = FALSE)
  }
  if (is.null(vars) && !isTRUE(full_table)) {
    stop("`vars` must be provided unless full_table = TRUE", call. = FALSE)
  }
  if (!is.null(vars) && isTRUE(full_table)) {
    warning("`vars` is ignored when full_table = TRUE", call. = FALSE)
  }
  df <- .jl_table(.jl_call("fetch_metadata", object, as.character(vars), cells, isTRUE(full_table)))
  if (return_class == "vector") {
    v <- df[[1]]
    names(v) <- rownames(df)
    return(v)
  }
  df
}
.jl_fetch_reduction <- function(object, reduction, cells = NULL, dims = c(1, 2)) {
  .jl_table(.jl_call("fetch_reduction", object, reduction, cells, as.integer(dims)))
}
.jl_fetch_cells <- function(object, meta_var, meta_levels) {
  if (is.null(meta_levels)) stop("`meta_levels` must be provided", call. = FALSE)
  .jl_value(.jl_call("fetch_cells", object, meta_var, as.character(meta_levels)))
}
.jl_get_all_cells <- function(object) .jl_value(.jl_call("get_all_cells", object))
.jl_meta_varnames <- function(object) .jl_value(.jl_call("meta_varnames", object))
.jl_all_keys <- function(object) {
  res <- .jl_call("all_keys", object)
  .jl_emit_warnings(res[[3]])
  out <- unlist(res[[2]])
  names(out) <- unlist(res[[1]])
  out
}
.jl_features_in_assay <- function(object, assay) .jl_value(.jl_call("features_in_assay", object, assay))
.jl_reduction_dimnames <- function(object, reduction, dims) {
  .jl_value(.jl_call("reduction_dimnames", object, reduction, as.integer(dims)))
}
.jl_default_layer <- function(object) .jl_value(.jl_call("default_layer", object))
.jl_default_reduction <- function(object) .jl_value(.jl_call("default_reduction_name", object))
.jl_unique_values <- function(object, var) .jl_value(.jl_call("unique_values", object, var))
.jl_n_unique <- function(object, var) .jl_value(.jl_call("n_unique", object, var))

# ---------------------------------------------------------------------------
# System image support (cold-start latency)

#' Build a Julia system image with SCUBA.jl compiled in
#'
#' Compiles SCUBA.jl and its dependencies into a Julia system image, which cuts the Julia
#' start-up of each R process from about 5 s to about 3 s. RCall itself cannot be compiled into
#' the image on Julia 1.13 (PackageCompiler fails with an LLVM calling-convention error), so it
#' and the RCall extension load from their own precompile caches. Building takes several
#' minutes once and runs in a separate `julia` process (equivalent to
#' `julia <SCUBA.jl>/tools/build_sysimage.jl --no-rcall --no-pythoncall` from a terminal).
#' Afterwards set `options(SCUBA.sysimage = path)` (or the `SCUBA_SYSIMAGE` environment
#' variable) before the first SCUBA call in a session.
#'
#' @param path Output path for the image (default: `~/.julia/scuba/SCUBA.<os>-<arch>.<ext>`).
#' @param julia_path Path to the Julia binary (default: the one JuliaCall finds).
#' @return The image path, invisibly.
#' @export
scuba_build_sysimage <- function(path = NULL, julia_path = NULL) {
  dev_path <- getOption("SCUBA.jl_path", Sys.getenv("SCUBA_JL_PATH", unset = ""))
  if (!nzchar(dev_path)) {
    stop("scuba_build_sysimage() needs a SCUBA.jl checkout: set options(SCUBA.jl_path = ...) or SCUBA_JL_PATH", call. = FALSE)
  }
  julia <- julia_path %||% Sys.which("julia")
  if (!nzchar(julia)) stop("julia binary not found", call. = FALSE)
  script <- file.path(normalizePath(dev_path), "tools", "build_sysimage.jl")
  args <- c(paste0("--project=", file.path(normalizePath(dev_path), "tools", "sysimage")), script, "--no-rcall", "--no-pythoncall")
  if (!is.null(path)) args <- c(args, "--out", path.expand(path))
  message("Building the SCUBA system image (this takes several minutes) ...")
  out <- system2(julia, args, stdout = TRUE, stderr = "")
  img <- tail(out[nzchar(out)], 1)
  message("System image written to ", img, "\nUse it with options(SCUBA.sysimage = \"", img, "\") before the first SCUBA call.")
  invisible(img)
}

.scuba_sysimage <- function() {
  p <- getOption("SCUBA.sysimage", Sys.getenv("SCUBA_SYSIMAGE", unset = ""))
  if (!nzchar(p)) return(NULL)
  if (!file.exists(path.expand(p))) {
    warning("SCUBA.sysimage / SCUBA_SYSIMAGE points to a file that does not exist (", p,
            "); starting Julia without a system image", call. = FALSE)
    return(NULL)
  }
  normalizePath(path.expand(p))
}
