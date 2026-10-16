# Plots and tables of all50 original replications.
study_reps <- 50L
study_write_tables <- function(cells,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  select<-function(b,n,v)cells[cells$beta==b&cells$n==n&cells$variant==v,]
  val<-function(x,col,digits=4L)if(x$complete)sprintf(paste0('%.',digits,'f'),x[[col]])else'--'
  withse<-function(x,normalized=FALSE)if(!x$complete)sprintf('-- [%d/50]',x$accepted_count)else
    if(normalized)sprintf('%.3f (%.3f)',x$mean_risk_over_q,x$se_risk_over_q)else
      sprintf('%.4f (%.4f)',x$mean_risk,x$se_risk)
  common<-r"(Each estimate uses 50 independent replications. Parentheses give Monte Carlo standard errors. Supplementary Tables~\ref{tab:supp-half} and~\ref{tab:supp-lipschitz} give pointwise 95\% normal Monte Carlo confidence intervals.)"
  main_caption<-sub('Parentheses give Monte Carlo standard errors.',
    'The SE column gives Monte Carlo standard errors.',common,fixed=TRUE)
  main<-c(r"(\begin{table}[t])",paste0(r"(\caption{Squared prediction error with Gaussian parameters and geometric width. )",main_caption,'}'),
    r"(\label{tab:simulation-adaptive-risk})",r"(\centering\small)",
    r"(\begin{tabular}{@{}crcccc@{}})",r"(\toprule)",
    r"($\beta$ & $n$ & Mean $\widetilde B_n$ & Error & SE & Error/$q_{n,\beta}$\\)",r"(\midrule)")
  width<-c(r"(\begin{table}[htbp])",paste0(r"(\caption{Squared prediction error divided by $q_{n,\beta}=(\log n/n)^{2\beta/(2\beta+1)}$ under three width laws with Gaussian parameters. )",common,'}'),
    r"(\label{tab:width-robustness})",r"(\centering\small)",r"(\begin{tabular}{@{}crrrr@{}})",
    r"(\toprule)",r"($\beta$ & $n$ & Geometric & Shifted Poisson & Power law\\)",r"(\midrule)")
  direct<-c(r"(\begin{table}[htbp])",paste0(r"(\caption{Squared prediction error for uniform and Student $t_{10}$ parameters with geometric width. )",common,'}'),
    r"(\label{tab:direct-family-adaptive-risk})",r"(\centering\small)",r"(\begin{tabular}{@{}crrrrr@{}})",
    r"(\toprule)",r"(& & \multicolumn{2}{c}{Uniform parameters} & \multicolumn{2}{c}{Student $t_{10}$ parameters}\\)",
    r"(\cmidrule(lr){3-4}\cmidrule(l){5-6})",r"($\beta$ & $n$ & Error & Error$/q_{n,\beta}$ & Error & Error$/q_{n,\beta}$\\)",r"(\midrule)")
  for(b in c(.5,1))for(n in study_ns) {
    if(b==1&&n==min(study_ns)){main<-c(main,r"(\midrule)");width<-c(width,r"(\midrule)");direct<-c(direct,r"(\midrule)")}
    label<-if(n==min(study_ns))if(b==.5)'$1/2$'else'$1$'else''
    x<-select(b,n,'N/G')
    main<-c(main,paste0(paste(c(label,n,val(x,'mean_Bhat',2),
      val(x,'mean_risk'),val(x,'se_risk'),val(x,'mean_risk_over_q',3)),collapse=' & '),r"(\\)"))
    width<-c(width,paste0(paste(c(label,n,vapply(study_variants[1:3],function(v)withse(select(b,n,v),TRUE),'')),collapse=' & '),r"(\\)"))
    direct<-c(direct,paste0(paste(c(label,n,unlist(lapply(study_variants[4:5],function(v){x<-select(b,n,v);c(withse(x),withse(x,TRUE))}))),collapse=' & '),r"(\\)"))
  }
  finish<-function(z,star=FALSE)c(z,r"(\bottomrule)",r"(\end{tabular})",if(star)r"(\end{table*})"else r"(\end{table})")
  writeLines(finish(main),file.path(out,'adaptive_rate_table.tex'))
  writeLines(finish(width),file.path(out,'width_robustness_table.tex'))
  writeLines(finish(direct),file.path(out,'direct_family_adaptive_table.tex'))
  study_write_supplementary_tables(cells,out)
  write.csv(cells,file.path(out,'final_adaptive_risk_cells.csv'),row.names=FALSE)
}
study_write_supplementary_tables <- function(cells,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  target_names<-c(half=r"(H\"older $1/2$)",one='Lipschitz')
  target_labels<-c(half='half',one='lipschitz')
  for(target in study_targets) {
    suffix<-target_labels[[target]]
    caption<-paste0('Supplementary squared-error estimates for the ',target_names[[target]],
      r"( target. Each mean uses 50 replications; SE is the Monte Carlo standard error and CI is the pointwise 95\% normal Monte Carlo confidence interval. N/U/$t_{10}$ denote Gaussian/uniform/Student parameters; G/P/H denote geometric/shifted Poisson/power-law width.)")
    table<-c(r"(\begin{table}[H])",paste0(r"(\caption{)",caption,'}'),
      paste0(r"(\label{tab:supp-)",suffix,'}'),r"(\centering\small)",
      r"(\begin{tabular}{@{}rlrrrr@{}})",r"(\toprule)",
      r"($n$ & Prior & Mean error & SE & \multicolumn{2}{c}{95\% CI}\\)",
      r"(\cmidrule(l){5-6})",r"( & & & & Lower & Upper\\)",r"(\midrule)")
    for(n in study_ns) {
      if(n!=min(study_ns))table<-c(table,r"(\addlinespace[3pt])")
      for(variant in study_variants) {
        x<-cells[cells$target==target&cells$n==n&cells$variant==variant,]
        stopifnot(nrow(x)==1L,isTRUE(x$complete),x$accepted_count==50L)
        label<-if(variant=='t10/G')'$t_{10}$/G'else variant
        values<-sprintf('%.6f',unlist(x[c('mean_risk','se_risk','ci_lower','ci_upper')],use.names=FALSE))
        table<-c(table,paste0(paste(c(if(variant==study_variants[1L])n else '',label,values),
          collapse=' & '),r"(\\)"))
      }
    }
    table<-c(table,r"(\bottomrule)",r"(\end{tabular})",r"(\end{table})")
    writeLines(table,file.path(out,paste0('supplementary_',suffix,'_table.tex')))
  }
}
study_draw_plot <- function(cells,normalized=FALSE) {
  colors<-c('#0072B2','#E69F00','#009E73','#CC79A7','#D55E00')
  symbols<-c(16,17,15,18,8);types<-c(1,5,4,3,6);dx<-seq(-.065,.065,length.out=5L)
  layout(matrix(c(1,2,3,3),2,byrow=TRUE),heights=c(1,.17))
  par(oma=c(0,0,2.2,0),mar=c(2.5,4.6,2.3,1),mgp=c(2.7,.65,0),tcl=-.25)
  for(b in c(.5,1)) {
    c<-cells[cells$beta==b,];good<-c[c$complete,]
    anchor<-c[c$variant=='N/G'&c$n==min(study_ns),]
    upper_guide<-lower_guide<-numeric()
    if(anchor$complete) {
      upper_guide<-anchor$mean_risk*study_q(study_ns,b)/study_q(min(study_ns),b)
      lower_guide<-anchor$mean_risk*(study_ns/min(study_ns))^(-2*b/(2*b+1))
      if(normalized) {
        upper_guide<-upper_guide/study_q(study_ns,b)
        lower_guide<-lower_guide/study_q(study_ns,b)
      }
    }
    if(normalized) {
      yy<-c(good$normalized_ci_lower,good$normalized_ci_upper,upper_guide,lower_guide)
      yy<-yy[is.finite(yy)&yy>0]
      ylab<-'Mean squared error / q';ylim<-if(length(yy))range(yy)*c(.85,1.18)else c(.1,10)
    }else{
      yy<-c(good$ci_lower,good$ci_upper,upper_guide,lower_guide)
      yy<-yy[is.finite(yy)&yy>0]
      ylab<-'Mean integrated squared error';ylim<-if(length(yy))range(yy)*c(.8,1.24)else c(.001,1)
    }
    plot(NA,xlim=range(log2(study_ns))+c(-.19,.19),ylim=ylim,log='y',xaxt='n',xlab='',ylab=ylab,
         main=if(b==.5)expression(beta==1/2)else expression(beta==1),cex.main=1.05)
    axis(1,at=log2(study_ns),labels=study_ns,cex.axis=.77)
    abline(h=10^seq(-6,3),col='#EEEEEE',lwd=.6)
    if(anchor$complete) {
      lines(log2(study_ns),lower_guide,col='#333333',lty=2,lwd=1.8)
      lines(log2(study_ns),upper_guide,col='#777777',lty=3,lwd=1.8)
    }
    for(j in seq_along(study_variants)) {
      z<-c[c$variant==study_variants[j],];z<-z[match(study_ns,z$n),]
      y<-if(normalized)z$mean_risk_over_q else z$mean_risk
      lo<-if(normalized)z$normalized_ci_lower else z$ci_lower
      hi<-if(normalized)z$normalized_ci_upper else z$ci_upper
      xx<-log2(z$n)+dx[j]
      lines(xx,y,col=colors[j],lty=types[j],lwd=1.5)
      ok<-z$complete&is.finite(lo)&lo>0
      segments(xx[ok],lo[ok],xx[ok],hi[ok],col=colors[j],lwd=1)
      segments(xx[ok]-.025,lo[ok],xx[ok]+.025,lo[ok],col=colors[j])
      segments(xx[ok]-.025,hi[ok],xx[ok]+.025,hi[ok],col=colors[j])
      points(xx,y,pch=symbols[j],col=colors[j],cex=.76)
    }
    if(any(!c$complete)) {
      counts<-vapply(study_variants[4:5],function(v)sum(c$complete[c$variant==v]),0L)
      mtext(sprintf('U/G: %d/7; t10/G: %d/7 cells complete; unavailable means withheld',
        counts[1L],counts[2L]),side=3,line=.05,cex=.65)
    }
  }
  par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
  text(.5,.87,'Sample size n',cex=.91)
  entry<-legend(x=.5,y=.50,xjust=.5,legend=c(study_variants,'Lower-rate guide','Upper-rate guide'),
    col=c(colors,'#333333','#777777'),pch=c(symbols,NA,NA),lty=c(types,2,3),
    lwd=c(rep(1.5,5L),1.8,1.8),bty='n',horiz=TRUE,cex=.74,x.intersp=1.0,seg.len=1.5,y.intersp=.9)
  stopifnot(entry$rect$left>=0,entry$rect$left+entry$rect$w<=1)
}
study_plot <- function(cells,pdf_path,allow_pending=FALSE) {
  stopifnot(is.logical(allow_pending),length(allow_pending)==1L,!is.na(allow_pending))
  if(!allow_pending)stopifnot(nrow(cells)==70L,all(cells$complete),
    all(cells$accepted_count==50L),all(cells$failed_count==0L),all(cells$pending_count==0L))
  stopifnot(all(is.na(cells$mean_risk[!cells$complete])),
    all(cells$accepted_count[cells$complete]==50L))
  dir.create(dirname(pdf_path),recursive=TRUE,showWarnings=FALSE)
  pdf(pdf_path,width=9.2,height=4.35,useDingbats=FALSE)
  study_draw_plot(cells);dev.off()
  png(sub('\\.pdf$','.png',pdf_path,ignore.case=TRUE),width=2208,height=1044,res=240)
  study_draw_plot(cells);dev.off()
  normalized<-file.path(dirname(pdf_path),'adaptive_normalized_error.png')
  png(normalized,width=2208,height=1044,res=240)
  study_draw_plot(cells,TRUE);dev.off()
}
