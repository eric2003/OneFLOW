#include "OneDEulerBackend.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double Gamma = 1.4;
constexpr double DtScale = 0.05;
using Clock = std::chrono::steady_clock;
int ParsePositive(const char *text,const char *name){char *end=nullptr;long v=std::strtol(text,&end,10);if(end==text||*end!='\0'||v<=0||v>100000000)throw std::invalid_argument(std::string("invalid ")+name);return (int)v;}
std::vector<double> InitialState(int nx){std::vector<double>q(3*nx);for(int i=0;i<nx;++i){double x=(i+.5)/nx,r=1+.1*std::sin(2*M_PI*x),u=.2+.05*std::cos(2*M_PI*x),p=1+.08*std::sin(4*M_PI*x);q[i]=r;q[nx+i]=r*u;q[2*nx+i]=p/(Gamma-1)+.5*r*u*u;}return q;}
double MaxAbs(const std::vector<double>&a,const std::vector<double>&b){double e=0;for(size_t i=0;i<a.size();++i)e=std::max(e,std::abs(a[i]-b[i]));return e;}
double Checksum(const std::vector<double>&v){double s=0;for(size_t i=0;i<v.size();++i)s+=v[i]*(1+(int)(i%17));return s;}
struct Result{std::vector<double> state;double createMs=0,uploadMs=0,advanceMs=0,downloadMs=0,kernelMs=0;int launches=0,syncs=0;};
Result Run(const oneflow_1d::EulerBackend&backend,const std::vector<double>&initial,int nx,int steps){oneflow_1d::EulerProblem p{nx,Gamma,DtScale/nx,1.0/nx,oneflow_1d::EulerBoundary::Periodic};auto a=Clock::now();auto state=backend.CreateState(p);auto b=Clock::now();backend.Upload(*state,initial.data());auto c=Clock::now();oneflow_1d::EulerRunStats stats;backend.Advance(*state,steps,{oneflow_1d::EulerRunMode::NoTrace,nullptr,&stats});auto d=Clock::now();std::vector<double>final(initial.size());backend.Download(*state,final.data());auto e=Clock::now();return{std::move(final),std::chrono::duration<double,std::milli>(b-a).count(),std::chrono::duration<double,std::milli>(c-b).count(),std::chrono::duration<double,std::milli>(d-c).count(),std::chrono::duration<double,std::milli>(e-d).count(),stats.kernelMilliseconds,stats.kernelLaunches,stats.synchronizationCount};}
}
int main(int argc,char**argv){try{int nx=argc>1?ParsePositive(argv[1],"nx"):65536,steps=argc>2?ParsePositive(argv[2],"steps"):100,repeats=argc>3?ParsePositive(argv[3],"repeats"):3,warmup=argc>4?ParsePositive(argv[4],"warmup"):1;auto initial=InitialState(nx);oneflow_1d::CpuEulerBackend cpu;oneflow_1d::HipEulerBackend hip;for(int i=0;i<warmup;++i){Run(cpu,initial,nx,steps);Run(hip,initial,nx,steps);}double cpuAdvance=0,hipAdvance=0,hipCreate=0,hipUpload=0,hipDownload=0,hipKernel=0;int hipLaunches=0,hipSyncs=0;std::vector<double>cpuFinal,hipFinal;for(int i=0;i<repeats;++i){auto c=Run(cpu,initial,nx,steps);auto h=Run(hip,initial,nx,steps);cpuAdvance+=c.advanceMs;hipAdvance+=h.advanceMs;hipKernel+=h.kernelMs;hipLaunches+=h.launches;hipSyncs+=h.syncs;hipCreate+=h.createMs;hipUpload+=h.uploadMs;hipDownload+=h.downloadMs;cpuFinal=std::move(c.state);hipFinal=std::move(h.state);}std::cout<<std::fixed<<std::setprecision(6);std::cout<<"OneFLOW 1D Euler stateful backend benchmark\n"<<"nx="<<nx<<" steps="<<steps<<" repeats="<<repeats<<" warmup="<<warmup<<"\n";std::cout<<"cpu_advance_ms="<<cpuAdvance<<" cpu_per_step_ms="<<cpuAdvance/repeats/steps<<"\n";std::cout<<"hip_advance_ms="<<hipAdvance<<" hip_per_step_ms="<<hipAdvance/repeats/steps<<"\n";std::cout<<"hip_kernel_ms="<<hipKernel<<" hip_kernel_ratio="<<hipKernel/hipAdvance<<" hip_kernel_launches="<<hipLaunches<<" hip_syncs="<<hipSyncs<<"\n";std::cout<<"hip_create_ms="<<hipCreate<<" hip_upload_ms="<<hipUpload<<" hip_download_ms="<<hipDownload<<"\n";std::cout<<"hip_steady_speedup="<<cpuAdvance/hipAdvance<<"\n";std::cout<<"final_max_abs_error="<<MaxAbs(cpuFinal,hipFinal)<<"\n";std::cout<<"cpu_checksum="<<Checksum(cpuFinal)<<" hip_checksum="<<Checksum(hipFinal)<<"\n";return 0;}catch(const std::exception&e){std::cerr<<"Euler stateful benchmark failed: "<<e.what()<<"\n";return 1;}}
