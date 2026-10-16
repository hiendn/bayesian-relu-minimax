// Native implementation of the reference direct R sampler.
// All proposal distributions, likelihood powers and Gibbs updates are unchanged.
// Index sampling uses R_unif_index, matching R 4.4 sample.int without replacement.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <R_ext/Random.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#ifdef __clang__
#pragma clang fp contract(off)
#endif
using namespace Rcpp;

namespace {
struct State {
  std::vector<double> a, u, b, pred;
  double rss, lp;
};
inline double relu(double x, double u, double b) {
  return std::max(u * x + b, 0.0);
}
inline int index_draw(int n) { return static_cast<int>(R_unif_index(n)); }
std::vector<int> permutation(int n) {
  std::vector<int> pool(n), out(n);
  for (int i=0; i<n; ++i) pool[i]=i;
  int left=n;
  for (int i=0; i<n; ++i) {
    int j=index_draw(left);
    out[i]=pool[j]; pool[j]=pool[--left];
  }
  return out;
}
double coord_lp(double z, bool uniform, int k, bool output) {
  if (!uniform) return R::dt(z, 10.0, true);
  double bound=output ? double(k)*k : 2.0;
  return std::abs(z)<=bound ? -std::log(2.0*bound) : R_NegInf;
}
double draw_coord(bool uniform, int k, bool output) {
  if (!uniform) return R::rt(10.0);
  double bound=output ? double(k)*k : 2.0;
  return R::runif(-bound,bound);
}
double log_prior(const State &s, bool uniform) {
  int k=s.a.size();
  long double aa=0, uu=0, bb=0;
  for (int j=0;j<k;++j) {
    aa += coord_lp(s.a[j],uniform,k,true);
    uu += coord_lp(s.u[j],uniform,k,false);
    bb += coord_lp(s.b[j],uniform,k,false);
  }
  return -k*std::log(2.0) + double(aa) + double(uu) + double(bb);
}
void evaluate(State &s, const NumericVector &x, const NumericVector &y,
              bool uniform) {
  int n=x.size(), k=s.a.size(); s.pred.assign(n,0.0);
  for (int j=0;j<k;++j)
    for (int i=0;i<n;++i) s.pred[i] += relu(x[i],s.u[j],s.b[j])*s.a[j];
  long double ss=0.0;
  for (int i=0;i<n;++i) { double r=y[i]-s.pred[i]; ss += r*r; }
  s.rss=double(ss); s.lp=log_prior(s,uniform);
}
double birth_prob(int k,int cap) { return k==cap ? 0.0 : k==1 ? 1.0 : .5; }
double death_prob(int k,int cap) { return 1.0-birth_prob(k,cap); }
bool dimension_move(State &s,const NumericVector &x,const NumericVector &y,
                    bool uniform,double lambda,int cap) {
  int k=s.a.size(); bool up=R::runif(0,1)<birth_prob(k,cap);
  State z=s; double log_q_ratio;
  if(up) {
    int j=index_draw(k+1);
    double a=draw_coord(uniform,k+1,true);
    double u=draw_coord(uniform,k+1,false);
    double b=draw_coord(uniform,k+1,false);
    z.a.insert(z.a.begin()+j,a); z.u.insert(z.u.begin()+j,u);
    z.b.insert(z.b.begin()+j,b);
    double lp=coord_lp(a,uniform,k+1,true) +
      coord_lp(u,uniform,k+1,false) + coord_lp(b,uniform,k+1,false);
    log_q_ratio=std::log(death_prob(k+1,cap))-
      std::log(birth_prob(k,cap))-lp;
  } else {
    int j=index_draw(k);
    double lp=coord_lp(s.a[j],uniform,k,true) +
      coord_lp(s.u[j],uniform,k,false) + coord_lp(s.b[j],uniform,k,false);
    z.a.erase(z.a.begin()+j); z.u.erase(z.u.begin()+j); z.b.erase(z.b.begin()+j);
    log_q_ratio=std::log(birth_prob(k-1,cap))-
      std::log(death_prob(k,cap))+lp;
  }
  evaluate(z,x,y,uniform);
  double ratio=-lambda*(z.rss-s.rss)/2.0+z.lp-s.lp+log_q_ratio;
  bool accepted=std::isfinite(ratio) && std::log(R::runif(0,1))<ratio;
  if(accepted) s=std::move(z);
  return accepted;
}
bool hidden_move(State &s,const NumericVector &x,const NumericVector &y,
                 bool uniform,double lambda) {
  int k=s.a.size(), j=index_draw(k); bool is_u=index_draw(2)==0;
  double old=is_u?s.u[j]:s.b[j]; bool refresh=R::runif(0,1)<.25;
  double proposal=refresh ? draw_coord(uniform,k,false) : old+R::rnorm(0,.4);
  State z=s; if(is_u) z.u[j]=proposal; else z.b[j]=proposal;
  evaluate(z,x,y,uniform);
  double ratio=-lambda*(z.rss-s.rss)/2.0+z.lp-s.lp;
  if(refresh) ratio += coord_lp(old,uniform,k,false)-coord_lp(proposal,uniform,k,false);
  bool accepted=std::isfinite(ratio) && std::log(R::runif(0,1))<ratio;
  if(accepted) s=std::move(z);
  return accepted;
}
double truncnorm_draw(double mu,double sd,double lower,double upper) {
  if(sd<=0) return std::max(lower,std::min(upper,mu));
  double lo=R::pnorm(lower,mu,sd,true,false);
  double hi=R::pnorm(upper,mu,sd,true,false);
  if(hi>lo+1e-12) return R::qnorm(R::runif(lo,hi),mu,sd,true,false);
  if(upper<mu) return -truncnorm_draw(-mu,sd,-upper,-lower);
  double alpha=(lower-mu)/sd;
  if(alpha>0) {
    double rate=(alpha+std::sqrt(alpha*alpha+4))/2;
    while(true) {
      double zz=alpha+R::rexp(1/rate);
      if(zz<=(upper-mu)/sd && R::runif(0,1)<std::exp(-.5*(zz-rate)*(zz-rate)))
        return mu+sd*zz;
    }
  }
  stop("Truncated-normal numeric failure"); return NA_REAL;
}
double uniform_output_draw(double lambda,double norm2,double cross,double bound) {
  double q=lambda*norm2;
  if(q==0) return R::runif(-bound,bound);
  double mu=cross/norm2;
  if(q<1e-6) {
    double mode=std::max(-bound,std::min(bound,mu));
    while(true) {
      double z=R::runif(-bound,bound);
      double log_accept=lambda*cross*(z-mode)-q*(z*z-mode*mode)/2;
      if(std::log(R::runif(0,1))<std::min(0.0,log_accept)) return z;
    }
  }
  return truncnorm_draw(mu,1/std::sqrt(q),-bound,bound);
}
void output_gibbs(State &s,const NumericVector &x,const NumericVector &y,
                  bool uniform,double lambda) {
  int n=x.size(), k=s.a.size(); std::vector<int> order=permutation(k);
  std::vector<double> phi(n);
  for(int j:order) {
    long double nn=0.0, cc=0.0;
    for(int i=0;i<n;++i) {
      phi[i]=relu(x[i],s.u[j],s.b[j]);
      nn += phi[i]*phi[i];
      double resid=y[i]-s.pred[i]+s.a[j]*phi[i];
      cc += phi[i]*resid;
    }
    double norm2=double(nn), cross=double(cc), new_a;
    if(uniform) new_a=uniform_output_draw(lambda,norm2,cross,double(k)*k);
    else {
      double prior_var=(10+s.a[j]*s.a[j])/(2*R::rgamma(5.5,1));
      double precision=lambda*norm2+1/prior_var;
      new_a=R::rnorm(lambda*cross/precision,1/std::sqrt(precision));
    }
    for(int i=0;i<n;++i) s.pred[i] += (new_a-s.a[j])*phi[i];
    s.a[j]=new_a;
  }
  long double ss=0.0;
  for(int i=0;i<n;++i) { double r=y[i]-s.pred[i]; ss += r*r; }
  s.rss=double(ss); s.lp=log_prior(s,uniform);
}
}

// [[Rcpp::export]]
List native_direct_chain(NumericVector x,NumericVector y,NumericVector grid,
                         NumericVector f0,std::string family,int iterations,
                         int burn,int thin,int Kmax,NumericVector lambdas,
                         int swap_every,double Bclip) {
  if(family!="uniform" && family!="student10") stop("Unknown family");
  if(x.size()!=y.size() || grid.size()!=f0.size() || iterations<=burn || thin<1)
    stop("Invalid chain inputs");
  if(Kmax<2 || lambdas.size()<2 || !(lambdas[lambdas.size()-1]>0 && lambdas[lambdas.size()-1]<1) || swap_every<1)
    stop("Invalid sampler parameters");
  bool uniform=family=="uniform"; int nt=lambdas.size(), ng=grid.size();
  std::vector<State> states(nt);
  for(int t=0;t<nt;++t) {
    int k=std::min(Kmax,1+int(R::rgeom(.5))); State &s=states[t];
    s.a.resize(k); s.u.resize(k); s.b.resize(k);
    for(int j=0;j<k;++j) s.a[j]=draw_coord(uniform,k,true);
    for(int j=0;j<k;++j) s.u[j]=draw_coord(uniform,k,false);
    for(int j=0;j<k;++j) s.b[j]=draw_coord(uniform,k,false);
    evaluate(s,x,y,uniform);
  }
  int kept=(iterations-burn)/thin, ix=0, at_cap=0;
  IntegerVector width(kept); NumericMatrix probes(kept,5);
  List kept_states(kept); // Recording only: no random draws or transition changes.
  NumericVector mean_f(ng), accepted(2), proposed(2), acc_swap(nt-1), prop_swap(nt-1);
  int probe_idx[5]; double fracs[5]={.1,.3,.5,.7,.9};
  for(int j=0;j<5;++j) probe_idx[j]=std::max(1,std::min(ng,int(std::nearbyint(fracs[j]*ng))))-1;
  auto t0=std::chrono::steady_clock::now();
  for(int iter=1;iter<=iterations;++iter) {
    if(iter%1000==0) checkUserInterrupt();
    for(int t=0;t<nt;++t) {
      bool dim=R::runif(0,1)<.25; int jj=dim?0:1;
      bool accept=dim ? dimension_move(states[t],x,y,uniform,lambdas[t],Kmax) :
        hidden_move(states[t],x,y,uniform,lambdas[t]);
      proposed[jj]++; accepted[jj] += accept;
      output_gibbs(states[t],x,y,uniform,lambdas[t]);
    }
    if(iter%swap_every==0) {
      int first=((iter/swap_every)%2==1)?0:1;
      for(int t=first;t<nt-1;t+=2) {
        double ratio=(lambdas[t]-lambdas[t+1])*(states[t].rss-states[t+1].rss)/2;
        prop_swap[t]++;
        if(std::log(R::runif(0,1))<ratio) {
          std::swap(states[t],states[t+1]); acc_swap[t]++;
        }
      }
    }
    const State &s=states[nt-1];
    if(iter>burn && int(s.a.size())==Kmax) at_cap++;
    if(iter>burn && (iter-burn)%thin==0) {
      width[ix]=s.a.size();
      kept_states[ix]=List::create(
        _["a"]=NumericVector(s.a.begin(),s.a.end()),
        _["u"]=NumericVector(s.u.begin(),s.u.end()),
        _["b"]=NumericVector(s.b.begin(),s.b.end()));
      std::vector<double> f(ng,0.0);
      for(size_t j=0;j<s.a.size();++j)
        for(int i=0;i<ng;++i) f[i] += s.a[j]*relu(grid[i],s.u[j],s.b[j]);
      for(int j=0;j<5;++j) probes(ix,j)=f[probe_idx[j]];
      for(int i=0;i<ng;++i) mean_f[i] += std::max(-Bclip,std::min(Bclip,f[i]));
      ix++;
    }
  }
  long double se=0;
  for(int i=0;i<ng;++i) { mean_f[i]/=kept; double r=mean_f[i]-f0[i]; se+=r*r; }
  for(int j=0;j<2;++j) accepted[j]/=proposed[j];
  for(int t=0;t<nt-1;++t) acc_swap[t]/=std::max(1.0,prop_swap[t]);
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  return List::create(_["width"]=width,_["pred_probes"]=probes,
    _["kept_states"]=kept_states,
    _["posterior_mean"]=mean_f,_["squared_error"]=double(se)/ng,
    _["dim_accept"]=accepted[0],_["hidden_accept"]=accepted[1],
    _["swap_accept"]=acc_swap,_["cap_fraction"]=double(at_cap)/(iterations-burn),
    _["seconds"]=seconds);
}
