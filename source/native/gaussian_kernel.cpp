// Native implementation of the reference collapsed Gaussian transitions.
// Distribution draws and singleton indices use R's RNG and R_unif_index.
// No prior or tuning constants are changed.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Random.h>
using namespace Rcpp;

int g_sample(int n) { return int(R_unif_index(double(n))); }

struct GState {
  arma::vec u,b,m;
  arma::mat R;
  double lm;
  int k() const { return u.n_elem; }
};

arma::mat g_design(const arma::vec &x,const arma::vec &u,const arma::vec &b) {
  arma::mat p(x.n_elem,u.n_elem);
  for(unsigned j=0;j<u.n_elem;j++) for(unsigned i=0;i<x.n_elem;i++)
    p(i,j)=std::max(0.0,x[i]*u[j]+b[j]);
  return p;
}

GState g_eval(const arma::vec &x,const arma::vec &y,
              const arma::vec &u,const arma::vec &b,double lambda) {
  GState s; s.u=u;s.b=b;
  int k=u.n_elem;
  arma::mat p=g_design(x,u,b);
  double tau2=std::pow(double(k),3);
  arma::mat Q=lambda*(p.t()*p);
  Q.diag()+=1.0/tau2;
  if(!arma::chol(s.R,Q)) stop("Gaussian precision Cholesky failed");
  arma::vec rhs=lambda*(p.t()*y);
  s.m=arma::solve(arma::trimatu(s.R),
        arma::solve(arma::trimatl(s.R.t()),rhs));
  double ld=k*std::log(tau2)+2.0*arma::accu(arma::log(s.R.diag()));
  s.lm=-0.5*ld+0.5*arma::dot(rhs,s.m);
  return s;
}

double g_bp(int k,int cap) { return k>=cap?0.0:(k<=1?1.0:0.5); }
double g_dp(int k,int cap) { return 1.0-g_bp(k,cap); }
arma::vec g_insert(const arma::vec &a,double z,int j) {
  arma::vec out(a.n_elem+1);
  for(int i=0;i<int(out.n_elem);i++) out[i]=i<j?a[i]:(i==j?z:a[i-1]);
  return out;
}
arma::vec g_remove(const arma::vec &a,int j) {
  arma::vec out(a.n_elem-1);
  for(int i=0;i<int(out.n_elem);i++) out[i]=a[i<j?i:i+1];
  return out;
}

arma::vec g_predict(const GState &s,const arma::vec &grid,double B) {
  arma::mat P=g_design(grid,s.u,s.b);
  arma::vec mu=P*s.m;
  arma::mat tmp=arma::solve(arma::trimatl(s.R.t()),P.t());
  arma::vec v=arma::sum(arma::square(tmp),0).t();
  arma::vec out(grid.n_elem);
  for(unsigned i=0;i<grid.n_elem;i++) {
    double sd=std::sqrt(std::max(v[i],0.0));
    if(sd<=1e-10) out[i]=std::max(-B,std::min(B,mu[i]));
    else {
      double a=(-B-mu[i])/sd,b=(B-mu[i])/sd;
      double pa=R::pnorm(a,0,1,1,0),pb=R::pnorm(b,0,1,1,0);
      out[i]=-B*pa+B*(1-pb)+mu[i]*(pb-pa)+
        sd*(R::dnorm(a,0,1,0)-R::dnorm(b,0,1,0));
    }
  }
  return out;
}

// [[Rcpp::export]]
List gaussian_chain_native(arma::vec x,arma::vec y,arma::vec grid,int seed,
  int Kmax,int iterations,int burn,int thin,NumericVector lambdas,
  int swap_every,int init_k,double Bclip) {
  Function setseed("set.seed");setseed(seed);
  RNGScope scope;GetRNGstate();
  int nt=lambdas.size(),nk=(iterations-burn)/thin,ng=grid.n_elem;
  std::vector<GState> ss;ss.reserve(nt);
  for(int t=0;t<nt;t++) {
    int k=t==nt-1?init_k:std::min(Kmax,1+int(R::rgeom(0.5)));
    arma::vec u(k),b(k);
    for(int j=0;j<k;j++) u[j]=R::rnorm(0,1);
    for(int j=0;j<k;j++) b[j]=R::rnorm(0,1);
    ss.push_back(g_eval(x,y,u,b,lambdas[t]));
  }
  IntegerVector label=seq(0,nt-1),stage(nt),roundtrip(nt),coldvisit(nt),
    poststage(nt),posttrip(nt),postvisit(nt),propSwap(nt-1),accSwap(nt-1),
    propMove(3),accMove(3),widthTrace(iterations);
  stage[label[0]]=1;
  NumericVector widths(nk),predTrace(nk);
  NumericMatrix probes(nk,5),predictions(ng,nk);
  List states(nk);
  arma::vec sum(ng,arma::fill::zeros);
  int idx[5];double frac[]={.1,.3,.5,.7,.9};
  for(int j=0;j<5;j++) idx[j]=std::max(1,std::min(ng,int(std::nearbyint(frac[j]*ng))))-1;
  int capcount=0,kept=0;
  for(int iter=1;iter<=iterations;iter++) {
    for(int t=0;t<nt;t++) {
      GState cand;GState &s=ss[t];double lr;int mt;
      if(R::runif(0,1)<.25) {
        int k=s.k();bool up=R::runif(0,1)<g_bp(k,Kmax);
        int j=g_sample(up?k+1:k);
        if(up) {
          // R evaluates the u argument before b.
          double nu=R::rnorm(0,1),nb=R::rnorm(0,1);
          cand=g_eval(x,y,g_insert(s.u,nu,j),g_insert(s.b,nb,j),lambdas[t]);
          lr=cand.lm-s.lm-std::log(2.0)+std::log(g_dp(k+1,Kmax))-std::log(g_bp(k,Kmax));
        } else {
          cand=g_eval(x,y,g_remove(s.u,j),g_remove(s.b,j),lambdas[t]);
          lr=cand.lm-s.lm+std::log(2.0)+std::log(g_bp(k-1,Kmax))-std::log(g_dp(k,Kmax));
        }
        mt=0;
      } else {
        int j=g_sample(s.k());
        arma::vec u=s.u,b=s.b;
        if(R::runif(0,1)<.2) {
          u[j]=R::rnorm(0,1);b[j]=R::rnorm(0,1);
          cand=g_eval(x,y,u,b,lambdas[t]);lr=cand.lm-s.lm;mt=2;
        } else {
          bool cu=R::runif(0,1)<.5;
          double old=cu?u[j]:b[j],nv=old+R::rnorm(0,.65);
          if(cu)u[j]=nv;else b[j]=nv;
          cand=g_eval(x,y,u,b,lambdas[t]);
          lr=cand.lm-s.lm-.5*(nv*nv-old*old);mt=1;
        }
      }
      bool accepted=std::log(R::runif(0,1))<lr;
      propMove[mt]++;if(accepted){ss[t]=std::move(cand);accMove[mt]++;}
    }
    if(iter%swap_every==0) {
      int first=(iter/swap_every)%2==1?0:1;
      for(int t=first;t<nt-1;t+=2) {
        GState left=g_eval(x,y,ss[t+1].u,ss[t+1].b,lambdas[t]);
        GState right=g_eval(x,y,ss[t].u,ss[t].b,lambdas[t+1]);
        double lr=left.lm+right.lm-ss[t].lm-ss[t+1].lm;
        propSwap[t]++;
        if(std::log(R::runif(0,1))<lr) {
          ss[t]=std::move(left);ss[t+1]=std::move(right);accSwap[t]++;
          int q=label[t];label[t]=label[t+1];label[t+1]=q;
        }
      }
    }
    int hi=label[0],ci=label[nt-1];
    if(stage[hi]==2){roundtrip[hi]++;stage[hi]=1;}else if(stage[hi]==0)stage[hi]=1;
    if(stage[ci]==1){stage[ci]=2;coldvisit[ci]++;}
    if(iter==burn+1)poststage[hi]=1;
    if(iter>burn) {
      if(poststage[hi]==2){posttrip[hi]++;poststage[hi]=1;}else if(poststage[hi]==0)poststage[hi]=1;
      if(poststage[ci]==1){poststage[ci]=2;postvisit[ci]++;}
    }
    GState &cold=ss[nt-1];widthTrace[iter-1]=cold.k();
    if(iter>burn&&cold.k()==Kmax)capcount++;
    if(iter>burn&&(iter-burn)%thin==0) {
      arma::vec pp=g_predict(cold,grid,Bclip);sum+=pp;
      widths[kept]=cold.k();states[kept]=List::create(_["u"]=cold.u,_["b"]=cold.b);
      predTrace[kept]=pp[int(std::ceil(ng*.37))-1];
      for(int j=0;j<5;j++)probes(kept,j)=pp[idx[j]];
      for(int i=0;i<ng;i++)predictions(i,kept)=pp[i];
      kept++;
    }
    if(iter%2000==0)checkUserInterrupt();
  }
  NumericVector swapacc(nt-1),moveacc(3);
  for(int t=0;t<nt-1;t++)swapacc[t]=double(accSwap[t])/std::max(propSwap[t],1);
  for(int j=0;j<3;j++)moveacc[j]=double(accMove[j])/propMove[j];
  moveacc.attr("names")=CharacterVector::create("dimension","coordinate","refresh");
  return List::create(_["predictive_mean"]=sum/kept,_["kept_width"]=widths,
    _["kept_pred"]=predTrace,_["kept_probe"]=probes,_["kept_states"]=states,
    _["kept_predictions"]=predictions,_["width_trace"]=widthTrace,
    _["at_cap"]=double(capcount)/(iterations-burn),_["move_accept"]=moveacc,
    _["swap_accept"]=swapacc,_["round_trips"]=int(Rcpp::sum(roundtrip)),
    _["cold_visits_from_hot"]=int(Rcpp::sum(coldvisit)),
    _["post_burn_round_trips"]=int(Rcpp::sum(posttrip)),
    _["post_burn_cold_visits"]=int(Rcpp::sum(postvisit)),_["seed"]=seed);
}
