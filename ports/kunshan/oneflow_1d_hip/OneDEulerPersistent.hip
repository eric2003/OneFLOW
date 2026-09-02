#include "OneDEulerBackend.h"
#include <hip/hip_runtime.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace oneflow_1d {
namespace {
void Check(hipError_t e,const char *s){if(e!=hipSuccess)throw std::runtime_error(std::string(s)+": "+hipGetErrorString(e));}
void CheckProblem(const EulerProblem&p){if(p.nx<=0||p.gamma<=1.0||p.dt<=0.0||p.dx<=0.0)throw std::invalid_argument("invalid Euler problem");}
void CheckOptions(int steps,const EulerRunOptions&o){if(steps<=0)throw std::invalid_argument("Euler steps must be positive");if(o.mode==EulerRunMode::FullTrace&&(steps!=1||o.trace==nullptr))throw std::invalid_argument("FullTrace requires one step and a trace");if(o.mode==EulerRunMode::NoTrace&&o.trace!=nullptr)throw std::invalid_argument("NoTrace cannot provide a trace");}
__device__ int CI(int c,int i,int nx){return c*nx+i;}
__device__ int FI(int c,int i,int nx){return c*(nx+1)+i;}
__device__ void PF(const double*q,int i,int nx,double g,double f[3],double&a){double r=q[CI(0,i,nx)],m=q[CI(1,i,nx)],e=q[CI(2,i,nx)],u=m/r,p=(g-1.0)*(e-0.5*m*m/r);f[0]=m;f[1]=m*u+p;f[2]=(e+p)*u;a=fabs(u)+sqrt(g*p/r); return;}
__global__ void Face(const double*q,int nx,double g,int b,double*l,double*r,double*f){int face=blockIdx.x*blockDim.x+threadIdx.x;if(face>nx)return;int il=face-1,ir=face;if(b==0){il=(face+nx-1)%nx;ir=face%nx;}else{il=il<0?0:il;ir=ir>=nx?nx-1:ir;}double fl[3],fr[3],al,ar;PF(q,il,nx,g,fl,al);PF(q,ir,nx,g,fr,ar);double a=al>ar?al:ar;for(int c=0;c<3;++c){double ql=q[CI(c,il,nx)],qr=q[CI(c,ir,nx)];int k=FI(c,face,nx);if(l!=nullptr)l[k]=ql;if(r!=nullptr)r[k]=qr;f[k]=0.5*(fl[c]+fr[c])-0.5*a*(qr-ql);}}
__global__ void ResidualUpdate(const double*f,int nx,double id,const double*base,const double*cur,int n,double dt,double bw,double uw,double*res,double*next){int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=n)return;int c=k/nx,i=k-c*nx;double value=-(f[FI(c,i+1,nx)]-f[FI(c,i,nx)])*id;if(res!=nullptr)res[k]=value;next[k]=bw*base[k]+uw*(cur[k]+dt*value);}

struct HipState final: EulerState {
  explicit HipState(const EulerProblem&p):problem(p){
    CheckProblem(problem); int cv=3*p.nx,fv=3*(p.nx+1); size_t cb=(size_t)cv*sizeof(double),fb=(size_t)fv*sizeof(double);
    try{Check(hipStreamCreate(&stream),"stream");Check(hipEventCreate(&kernelBegin),"kernel begin event");Check(hipEventCreate(&kernelEnd),"kernel end event");Check(hipMalloc(&base,cb),"base");Check(hipMalloc(&work,cb),"work");Check(hipMalloc(&scratch,cb),"scratch");Check(hipMalloc(&flux,fb),"flux");}catch(...){Release();throw;}
  }
  ~HipState() override{Release();}
  void Release() noexcept{if(stream)static_cast<void>(hipStreamSynchronize(stream));for(double*p:{residual,flux,right,left,scratch,work,base})if(p)static_cast<void>(hipFree(p));if(kernelEnd)static_cast<void>(hipEventDestroy(kernelEnd));if(kernelBegin)static_cast<void>(hipEventDestroy(kernelBegin));if(stream)static_cast<void>(hipStreamDestroy(stream));stream=nullptr;kernelBegin=kernelEnd=nullptr;base=work=scratch=left=right=flux=residual=nullptr; return;}
  void EnsureTraceBuffers(){int fv=3*(problem.nx+1),cv=3*problem.nx;size_t fb=(size_t)fv*sizeof(double),cb=(size_t)cv*sizeof(double);try{if(left==nullptr)Check(hipMalloc(&left,fb),"left");if(right==nullptr)Check(hipMalloc(&right,fb),"right");if(residual==nullptr)Check(hipMalloc(&residual,cb),"residual");}catch(...){if(residual){static_cast<void>(hipFree(residual));residual=nullptr;}if(right){static_cast<void>(hipFree(right));right=nullptr;}if(left){static_cast<void>(hipFree(left));left=nullptr;}throw;}}
  EulerProblem problem; hipStream_t stream=nullptr; hipEvent_t kernelBegin=nullptr; hipEvent_t kernelEnd=nullptr; double*base=nullptr;double*work=nullptr;double*scratch=nullptr;double*left=nullptr;double*right=nullptr;double*flux=nullptr;double*residual=nullptr;bool uploaded=false;
};
HipState& HS(EulerState&s){auto*p=dynamic_cast<HipState*>(&s);if(!p)throw std::invalid_argument("not HIP Euler state");return *p;}
const HipState& HS(const EulerState&s){auto*p=dynamic_cast<const HipState*>(&s);if(!p)throw std::invalid_argument("not HIP Euler state");return *p;}
}

const char*HipEulerBackend::Name()const{return "HIP";}
bool HipEulerBackend::IsAccelerator()const{return true;}
std::unique_ptr<EulerState>HipEulerBackend::CreateState(const EulerProblem&p)const{return std::make_unique<HipState>(p);}
void HipEulerBackend::Upload(EulerState&s,const double*h)const{auto&x=HS(s);if(!h)throw std::invalid_argument("null Euler host state");size_t b=(size_t)(3*x.problem.nx)*sizeof(double);Check(hipMemcpyAsync(x.base,h,b,hipMemcpyHostToDevice,x.stream),"H2D");Check(hipStreamSynchronize(x.stream),"upload sync");x.uploaded=true;}

void HipEulerBackend::Advance(EulerState&s,int steps,const EulerRunOptions&o)const{
  CheckOptions(steps,o);auto&x=HS(s);if(!x.uploaded)throw std::logic_error("HIP Euler state was not uploaded");int nx=x.problem.nx,cv=3*nx,fv=3*(nx+1),th=256,fb=(fv+th-1)/th,cb=(cv+th-1)/th;size_t cbytes=(size_t)cv*sizeof(double),fbytes=(size_t)fv*sizeof(double);EulerTrace*t=o.trace;if(t){ResizeEulerTrace(*t,nx);x.EnsureTraceBuffers();}
  if(t){Check(hipMemcpyAsync(t->state.data(),x.base,cbytes,hipMemcpyDeviceToHost,x.stream),"initial D2H");Check(hipStreamSynchronize(x.stream),"initial trace sync");}
  if(o.stats && !t) Check(hipEventRecord(x.kernelBegin,x.stream),"kernel begin");
  for(int step=0;step<steps;++step){
    for(int stage=0;stage<3;++stage){double bw=stage==1?.75:(stage==2?1.0/3.0:0.0),uw=stage==0?1.0:(stage==1?.25:2.0/3.0);double*cur=stage==0?x.base:(stage==1?x.work:x.scratch);double*next=stage==0?x.work:(stage==1?x.scratch:x.work);
      hipLaunchKernelGGL(Face,dim3(fb),dim3(th),0,x.stream,cur,nx,x.problem.gamma,(int)x.problem.boundary,t?x.left:nullptr,t?x.right:nullptr,x.flux);Check(hipGetLastError(),"Face launch");
      hipLaunchKernelGGL(ResidualUpdate,dim3(cb),dim3(th),0,x.stream,x.flux,nx,1.0/x.problem.dx,x.base,cur,cv,x.problem.dt,bw,uw,t?x.residual:nullptr,next);Check(hipGetLastError(),"ResidualUpdate launch");
      if(t){Check(hipStreamSynchronize(x.stream),"trace sync");int fo=stage*fv,co=stage*cv;Check(hipMemcpy(t->faceLeft.data()+fo,x.left,fbytes,hipMemcpyDeviceToHost),"left D2H");Check(hipMemcpy(t->faceRight.data()+fo,x.right,fbytes,hipMemcpyDeviceToHost),"right D2H");Check(hipMemcpy(t->numericalFlux.data()+fo,x.flux,fbytes,hipMemcpyDeviceToHost),"flux D2H");Check(hipMemcpy(t->residual.data()+co,x.residual,cbytes,hipMemcpyDeviceToHost),"residual D2H");Check(hipMemcpy(t->state.data()+(stage+1)*cv,next,cbytes,hipMemcpyDeviceToHost),"state D2H");}
    }
    std::swap(x.base,x.work);
  }
  if(!t){
    if(o.stats){Check(hipEventRecord(x.kernelEnd,x.stream),"kernel end");Check(hipEventSynchronize(x.kernelEnd),"kernel sync");float ms=0.0F;Check(hipEventElapsedTime(&ms,x.kernelBegin,x.kernelEnd),"kernel elapsed");o.stats->kernelMilliseconds=ms;o.stats->kernelLaunches=steps*EulerRkStages*2;o.stats->synchronizationCount=1;}
    else Check(hipStreamSynchronize(x.stream),"advance sync");
  }
}
void HipEulerBackend::Download(const EulerState&s,double*h)const{const auto&x=HS(s);if(!x.uploaded)throw std::logic_error("HIP Euler state was not uploaded");if(!h)throw std::invalid_argument("null Euler host state");size_t b=(size_t)(3*x.problem.nx)*sizeof(double);Check(hipMemcpy(h,x.base,b,hipMemcpyDeviceToHost),"final D2H");}
void HipEulerBackend::Step(const double*state,int nx,double gamma,double dt,double dx,EulerBoundary b,EulerTrace&trace)const{EulerBackend::Step(state,nx,gamma,dt,dx,b,trace);}
}
