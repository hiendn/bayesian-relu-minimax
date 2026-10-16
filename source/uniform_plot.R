# Physical Figure 1 panel size and font settings.
uniform_plot <- function(records,stem) {
config <- list(width_inches=2*7.4/3,height_inches=3*2/2.17,pointsize=10,font_family='Helvetica',baseline_width_inches=7.4,baseline_height_inches=3,baseline_layout_heights=c(1,1,.17))
yrange <- range(records$truth,records$posterior_mean,records$credible_lower,records$credible_upper)
yrange <- yrange+c(-1,1)*diff(yrange)*.045
settings <- function()par(oma=c(0,0,0,0),mar=c(2.6,2.7,1.75,.55),mgp=c(1.45,.45,0),
  tcl=-.22,las=1,cex.axis=.90,cex.lab=.95,xaxs='i',yaxs='i')
geometry <- function()par(c('fin','pin','mai','cex','mex','ps','family'))
pdf(NULL,width=config$baseline_width_inches,height=config$baseline_height_inches,
  family=config$font_family,useDingbats=FALSE,pointsize=config$pointsize)
layout(matrix(c(1,2,3,4,5,6,7,7,7),nrow=3,byrow=TRUE),heights=config$baseline_layout_heights)
settings();plot.new();baseline_geometry <- geometry();invisible(dev.off())
pdf(paste0(stem,'.pdf'),width=config$width_inches,height=config$height_inches,
  family=config$font_family,useDingbats=FALSE,pointsize=config$pointsize)
layout(matrix(1:4,nrow=2,byrow=TRUE),heights=c(1,1));settings()
par(cex=baseline_geometry$cex,mex=baseline_geometry$mex)
actual_geometry <- list()
for(n in c(1024L,4096L))for(target in c('half','one')) {
  z <- records[records$n==n&records$target==target,]
  plot(z$x,z$truth,type='n',xlim=c(0,1),ylim=yrange,xlab='x',ylab=expression(f(x)),axes=FALSE)
  axis(1,at=c(0,.25,.5,.75,1),labels=c('0','.25','.5','.75','1'))
  axis(2,at=c(0,.5,1),labels=c('0','.5','1'));box(col='grey45')
  polygon(c(z$x,rev(z$x)),c(z$credible_lower,rev(z$credible_upper)),col='#DCDCDC',border=NA)
  lines(z$x,z$truth,col='#171717',lwd=1.25);lines(z$x,z$posterior_mean,col='#2166AC',lwd=1.45)
  title(main=substitute(paste(beta==B,',  ',n==N),list(B=if(target=='half')quote(1/2)else 1,N=n)),cex.main=.98,font.main=1,line=.55)
  actual_geometry[[paste(target,n)]] <- geometry()
}
invisible(dev.off())
for(g in actual_geometry)stopifnot(max(abs(g$fin-baseline_geometry$fin))<1e-12,
  max(abs(g$pin-baseline_geometry$pin))<1e-12,max(abs(g$mai-baseline_geometry$mai))<1e-12,
  identical(g$cex,baseline_geometry$cex),identical(g$mex,baseline_geometry$mex),
  identical(g$ps,baseline_geometry$ps),identical(g$family,baseline_geometry$family))

  invisible(paste0(stem,'.pdf'))
}
