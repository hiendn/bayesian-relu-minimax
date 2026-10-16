# Native sampler bindings. Compile once before workers fork.
# The caller supplies native_dir, eta=0.9 and the clipping bound Bclip.
.compiled <- getOption('elementary_native_functions')
if (is.null(.compiled)) {
  stopifnot(requireNamespace('Rcpp',quietly=TRUE),
            requireNamespace('RcppArmadillo',quietly=TRUE))
  .compiled <- new.env(parent=baseenv())
  .old_env <- Sys.getenv(c('PKG_CPPFLAGS','R_MAKEVARS_USER'),unset=NA_character_)
  .restore <- function() for (key in names(.old_env)) {
    if (is.na(.old_env[[key]])) Sys.unsetenv(key) else
      do.call(Sys.setenv,setNames(list(.old_env[[key]]),key))
  }
  tryCatch({
    if (Sys.info()[['sysname']]=='Darwin') {
      .sdk <- '/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk'
      .headers <- file.path(.sdk,'usr/include/c++/v1')
      if (file.exists(file.path(.headers,'cmath'))) Sys.setenv(
        PKG_CPPFLAGS=paste(Sys.getenv('PKG_CPPFLAGS'),paste0('-I',.headers)))
      # Some binary R installations retain an obsolete external Fortran path.
      # The bundled BLAS/LAPACK already link R's bundled Fortran libraries.
      if (!dir.exists('/opt/gfortran/lib') &&
          file.exists(file.path(R.home('lib'),'libgfortran.5.dylib')))
        Sys.setenv(R_MAKEVARS_USER=file.path(native_dir,'local_Makevars'))
    }
    for (.file in c('gaussian_kernel.cpp','direct_kernel.cpp'))
      Rcpp::sourceCpp(file.path(native_dir,.file),env=.compiled,
        cacheDir=file.path(native_dir,'cache'),showOutput=FALSE,rebuild=FALSE)
  },finally=.restore())
  options(elementary_native_functions=.compiled)
}
run_tempered_chain <- function(x,y,grid,seed,Kmax=64L,
    iterations=12000L,burn=6000L,thin=20L,
    lambdas=study_ladder,
    swap_every=5L,init_k=1L) {
  stopifnot(tail(lambdas,1)==eta,all(diff(lambdas)>0),iterations>burn,
            Kmax>=2L,thin>=1L,swap_every>=1L)
  started <- proc.time()[['elapsed']]
  result <- .compiled$gaussian_chain_native(x,y,grid,as.integer(seed),
    as.integer(Kmax),as.integer(iterations),as.integer(burn),as.integer(thin),
    lambdas,as.integer(swap_every),as.integer(init_k),Bclip)
  result$predictive_mean <- drop(result$predictive_mean)
  result$elapsed_seconds <- proc.time()[['elapsed']]-started
  result
}
run_chain <- function(x,y,grid,f0,family,seed,
    iterations=8000L,burn=4000L,thin=10L,Kmax=64L,
    lambdas=study_ladder,
    swap_every=5L) {
  stopifnot(tail(lambdas,1)==eta,all(diff(lambdas)>0),iterations>burn,
            Kmax>=2L,thin>=1L,swap_every>=1L)
  set.seed(seed)
  .compiled$native_direct_chain(x,y,grid,f0,family,as.integer(iterations),
    as.integer(burn),as.integer(thin),as.integer(Kmax),lambdas,
    as.integer(swap_every),Bclip)
}
