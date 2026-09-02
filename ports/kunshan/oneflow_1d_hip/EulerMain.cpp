#include "OneDEulerBackend.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {
constexpr double G=1.4, AT=1e-15, RT=1e-15;
constexpr double SmoothJumpTolerance=0.02;
struct M{double a=0,r=0;std::uint64_t u=0;bool ok=true;};
std::uint64_t bits(double x){std::uint64_t b;std::memcpy(&b,&x,8);return b>>63?~b:(b|(std::uint64_t(1)<<63));}
void add(M&m,double x,double y){if(!std::isfinite(x)||!std::isfinite(y)){m.ok=false;return;}double a=std::abs(y-x);m.a=std::max(m.a,a);m.r=std::max(m.r,a/std::max(1.0,std::abs(x)));auto p=bits(x),q=bits(y);m.u=std::max(m.u,p>q?p-q:q-p);if(a>AT+RT*std::max(1.0,std::abs(x)))m.ok=false;}
void merge(M&a,const M&b){a.a=std::max(a.a,b.a);a.r=std::max(a.r,b.r);a.u=std::max(a.u,b.u);a.ok=a.ok&&b.ok;return;}
std::vector<double> smooth(int n){std::vector<double>q(3*n);for(int i=0;i<n;++i){double x=(i+.5)/n,r=1+.1*sin(2*M_PI*x),u=.2+.05*cos(2*M_PI*x),p=1+.08*sin(4*M_PI*x);q[i]=r;q[n+i]=r*u;q[2*n+i]=p/(G-1)+.5*r*u*u;}return q;}
std::vector<double> sod(int n){std::vector<double>q(3*n);for(int i=0;i<n;++i){double x=(i+.5)/n,r=x<=.5?1:.125,p=x<=.5?1:.1;q[i]=r;q[n+i]=0;q[2*n+i]=p/(G-1);}return q;}
double primitive(const std::vector<double>&q,int off,int c,int n,int k){double r=q[off+c],m=q[off+n+c],e=q[off+2*n+c];return k==0?r:k==1?m/r:(G-1)*(e-.5*m*m/r);}
bool physical(const std::vector<double>&q,int n){for(int s=0;s<=3;++s)for(int i=0;i<n;++i){double r=primitive(q,s*3*n,i,n,0),p=primitive(q,s*3*n,i,n,2);if(!std::isfinite(r)||!std::isfinite(p)||r<=0||p<=0)return false;}return true;}
bool smooth_cell(const std::vector<double>&q,int off,int cell,int n){if(cell<=0||cell>=n-1)return false;for(int k:{0,2}){double v=primitive(q,off,cell,n,k);double vl=primitive(q,off,cell-1,n,k),vr=primitive(q,off,cell+1,n,k);if(std::abs(v-vl)>SmoothJumpTolerance*std::max(1.0,std::abs(v))||std::abs(v-vr)>SmoothJumpTolerance*std::max(1.0,std::abs(v)))return false;}return true;}
#ifdef ONEFLOW_1D_USE_HIP
bool run(const char*name,const oneflow_1d::EulerBackend&cpuBackend,const oneflow_1d::EulerBackend&hipBackend,std::vector<double>cpu,std::vector<double>hip,int n,int steps,oneflow_1d::EulerBoundary b){int fv=3*(n+1),cv=3*n;std::array<M,3> ml{},mr{},mf{},mres{},ms{},mp{};int excludedResidual=0;bool pass=true;for(int s=0;s<steps;++s){oneflow_1d::EulerTrace a,z;cpuBackend.Step(cpu.data(),n,G,.05/n,1.0/n,b,a);hipBackend.Step(hip.data(),n,G,.05/n,1.0/n,b,z);for(int c=0;c<3;++c){for(int stage=0;stage<3;++stage)for(int i=0;i<n+1;++i){add(ml[c],a.faceLeft[stage*fv+c*(n+1)+i],z.faceLeft[stage*fv+c*(n+1)+i]);add(mr[c],a.faceRight[stage*fv+c*(n+1)+i],z.faceRight[stage*fv+c*(n+1)+i]);add(mf[c],a.numericalFlux[stage*fv+c*(n+1)+i],z.numericalFlux[stage*fv+c*(n+1)+i]);}for(int stage=0;stage<3;++stage)for(int i=0;i<n;++i){if(b==oneflow_1d::EulerBoundary::Transmissive&&!smooth_cell(a.state,stage*cv,i,n)){if(c==0)++excludedResidual;continue;}add(mres[c],a.residual[stage*cv+c*n+i],z.residual[stage*cv+c*n+i]);}for(int stage=0;stage<=3;++stage)for(int i=0;i<n;++i){add(ms[c],a.state[stage*cv+c*n+i],z.state[stage*cv+c*n+i]);for(int k=0;k<3;++k)add(mp[k],primitive(a.state,stage*cv,i,n,k),primitive(z.state,stage*cv,i,n,k));}}for(int c=0;c<3;++c)pass=pass&&ml[c].ok&&mr[c].ok&&mf[c].ok&&mres[c].ok&&ms[c].ok;for(int k=0;k<3;++k)pass=pass&&mp[k].ok;pass=pass&&physical(a.state,n)&&physical(z.state,n);std::copy(a.state.begin()+3*cv,a.state.end(),cpu.begin());std::copy(z.state.begin()+3*cv,z.state.end(),hip.begin());}
std::cout<<name<<" steps="<<steps<<" criterion abs <= 1e-15 + 1e-15*max(1,abs(reference))\n";if(b==oneflow_1d::EulerBoundary::Transmissive)std::cout<<"  Sod residual smooth mask: jump <= "<<SmoothJumpTolerance<<", excluded cell-stage comparisons="<<excludedResidual<<"\n";for(int c=0;c<3;++c){std::cout<<"  U"<<c<<" left abs="<<ml[c].a<<" rel="<<ml[c].r<<" ulp="<<ml[c].u<<"\n";std::cout<<"  U"<<c<<" right abs="<<mr[c].a<<" rel="<<mr[c].r<<" ulp="<<mr[c].u<<"\n";std::cout<<"  U"<<c<<" flux abs="<<mf[c].a<<" rel="<<mf[c].r<<" ulp="<<mf[c].u<<"\n";std::cout<<"  U"<<c<<" residual(smooth) abs="<<mres[c].a<<" rel="<<mres[c].r<<" ulp="<<mres[c].u<<"\n";std::cout<<"  U"<<c<<" state abs="<<ms[c].a<<" rel="<<ms[c].r<<" ulp="<<ms[c].u<<"\n";}for(int k=0;k<3;++k)std::cout<<"  primitive "<<k<<" abs="<<mp[k].a<<" rel="<<mp[k].r<<" ulp="<<mp[k].u<<"\n";std::cout<<name<<" result: "<<(pass?"PASS":"FAIL")<<"\n";return pass;}
#else
bool cpu_self_test(const oneflow_1d::EulerBackend&backend){constexpr int n=31;std::vector<double>q(3*n);for(int i=0;i<n;++i){q[i]=1.25;q[n+i]=.375;q[2*n+i]=.9/(G-1.0)+.5*1.25*.09;}oneflow_1d::EulerTrace t;backend.Step(q.data(),n,G,.05/n,1.0/n,oneflow_1d::EulerBoundary::Periodic,t);return physical(t.state,n);}
#endif
}
int main(){std::cout<<std::setprecision(17)<<std::scientific;
#ifdef ONEFLOW_1D_USE_HIP
oneflow_1d::CpuEulerBackend cpuBackend;oneflow_1d::HipEulerBackend hipBackend;int n=257;bool a=run("smooth-periodic",cpuBackend,hipBackend,smooth(n),smooth(n),n,10,oneflow_1d::EulerBoundary::Periodic);bool b=run("sod-transmissive",cpuBackend,hipBackend,sod(n),sod(n),n,10,oneflow_1d::EulerBoundary::Transmissive);std::cout<<"OneFLOW 1D compressible Euler CPU/HIP validation: "<<((a&&b)?"PASS":"FAIL")<<"\n";return a&&b?0:1;
#else
oneflow_1d::CpuEulerBackend cpuBackend;bool ok=cpu_self_test(cpuBackend);std::cout<<"OneFLOW 1D Euler CPU self-test: "<<(ok?"PASS":"FAIL")<<"\n";return ok?0:1;
#endif
}
