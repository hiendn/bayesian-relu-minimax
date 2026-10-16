#!/usr/bin/env Rscript
# Fast reconstruction from final numerical records and retained band summaries.
args <- commandArgs(trailingOnly=TRUE)
if(!length(args)||!all(args[-1L]%in%'--check-only'))
  stop('Usage: Rscript Replot.R OUTPUT [--check-only]')
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1L])
root <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(root,'source/study.R'))
source(file.path(root,'source/plots.R'))
source(file.path(root,'source/credible.R'))
source(file.path(root,'source/uniform_plot.R'))
read <- function(name)read.csv(file.path(root,'records',name),stringsAsFactors=FALSE)
native <- read('selected_native_fits.csv');risk <- read('selected_risks.csv')
study_validate(native,risk)
cells <- study_cells(risk);rates <- study_rates(risk)
same <- function(actual,retained,key) {
  k <- function(x)do.call(paste,c(x[key],sep=':'))
  retained <- retained[match(k(actual),k(retained)),]
  stopifnot(nrow(actual)==nrow(retained),!anyNA(k(retained)))
  for(column in names(actual)) {
    if(is.numeric(actual[[column]])) {
      a <- actual[[column]];b <- retained[[column]]
      stopifnot(identical(is.na(a),is.na(b)),
        max(abs(a-b),na.rm=TRUE)<1e-10)
    } else stopifnot(identical(actual[[column]],retained[[column]]))
  }
}
same(cells,read('risk_cells.csv'),c('target','n','variant'))
same(rates,read('rate_intervals.csv'),c('target','variant','n_from','n_to','type'))
uniform <- read('uniform_credible_intervals_data.csv')
gaussian <- read('gaussian_credible_intervals_data.csv')
check_band <- function(z,targets,ns,count) {
  stopifnot(nrow(z)==count*1201L,setequal(z$target,targets),setequal(z$n,ns),
    all(z$replicate==1L),all(is.finite(unlist(z[c('x','truth','posterior_mean',
      'credible_lower','credible_upper')]))),all(z$credible_lower<=z$credible_upper))
  for(t in targets)for(n in ns) {
    p <- z[z$target==t&z$n==n,]
    stopifnot(nrow(p)==1201L,max(abs(p$x-seq.int(0L,1200L)/1200))<1e-14,
      max(abs(p$truth-study_plot_truth(p$x,t)))<1e-12)
  }
}
check_band(uniform,study_targets,c(1024L,4096L),4L)
check_band(gaussian,study_targets,c(1024L,4096L),4L)
if('--check-only'%in%args) {
  cat('Verified2100 native records,3500 risks,70 complete50-fit cases,70 paired rate intervals and8 credible panels. No plots written.\n')
  quit(save='no')
}
out <- args[1L];if(file.exists(out))stop('Output already exists; choose a new directory.')
dir.create(out,recursive=TRUE);out <- normalizePath(out,mustWork=TRUE)
study_plot(cells,file.path(out,'adaptive_squared_error.pdf'))
study_write_tables(cells,out)
write.csv(cells,file.path(out,'risk_cells.csv'),row.names=FALSE)
write.csv(rates,file.path(out,'rate_intervals.csv'),row.names=FALSE)
credible_plot(gaussian,file.path(out,'two_lobes_credible_intervals.pdf'),study_targets)
uniform_pdf <- uniform_plot(uniform,file.path(out,'uniform_lipschitz_credible_intervals'))
python <- Sys.getenv('PYTHON','python3')
for(figure_pdf in c(uniform_pdf,file.path(out,'two_lobes_credible_intervals.pdf'))) {
  status <- system2(python,c(shQuote(file.path(root,'source/normalize_figure1.py')),shQuote(figure_pdf)))
  if(status!=0L)stop('Credible-figure page-box normalization failed; Python with pypdf is required.')
}
poppler <- Sys.getenv('PDFTOPPM',Sys.which('pdftoppm'))
if(nzchar(poppler)) {
  status <- system2(poppler,c('-singlefile','-r','210','-png',shQuote(uniform_pdf),
    shQuote(sub('[.]pdf$','',uniform_pdf))),stdout=FALSE,stderr=file.path(out,'figure1_render.log'))
  if(status!=0L)stop('Figure1 PNG rendering failed.')
}
cat('Reconstructed the final study figures and tables.\n')
