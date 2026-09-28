#include "OneDEuler.h"
#include <hip/hip_runtime.h>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace oneflow_1d {
namespace {
void Check(hipError_t e, const char *s) { if (e != hipSuccess) throw std::runtime_error(std::string(s)+": "+hipGetErrorString(e)); return; }
__device__ int CI(int c,int i,int nx){return c*nx+i;}
__device__ int FI(int c,int i,int nx){return c*(nx+1)+i;}
__device__ void PF(const double *q,int i,int nx,double g,double f[3],double &a){
  double r=q[CI(0,i,nx)],m=q[CI(1,i,nx)],e=q[CI(2,i,nx)],u=m/r,p=(g-1.0)*(e-0.5*m*m/r);
  f[0]=m; f[1]=m*u+p; f[2]=(e+p)*u; a=fabs(u)+sqrt(g*p/r);
  return;
}
__global__ void Face(const double *q,int nx,double g,int b,double *l,double *r,double *f){
  int face=blockIdx.x*blockDim.x+threadIdx.x; if(face>nx)return; int il=face-1,ir=face;
  if(b==0){il=(face+nx-1)%nx;ir=face%nx;} else {il=il<0?0:il;ir=ir>=nx?nx-1:ir;}
  double fl[3],fr[3],al,ar; PF(q,il,nx,g,fl,al); PF(q,ir,nx,g,fr,ar); double a=al>ar?al:ar;
  for(int c=0;c<3;++c){double ql=q[CI(c,il,nx)],qr=q[CI(c,ir,nx)];int k=FI(c,face,nx);l[k]=ql;r[k]=qr;f[k]=0.5*(fl[c]+fr[c])-0.5*a*(qr-ql);}
  return;
}
__global__ void Residual(const double *f,int nx,double id,double *res){int k=blockIdx.x*blockDim.x+threadIdx.x,n=3*nx;if(k>=n)return;int c=k/nx,i=k-c*nx;res[k]=-(f[FI(c,i+1,nx)]-f[FI(c,i,nx)])*id;return;}
__global__ void Update(const double *base,const double *cur,const double *res,int n,double dt,double bw,double uw,double *next){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)next[i]=bw*base[i]+uw*(cur[i]+dt*res[i]);return;}
}
void OneDHipEulerStep(const double *in,int nx,double g,double dt,double dx,EulerBoundary b,EulerTrace &t){
 if(!in||nx<=0||g<=1.0||dt<=0||dx<=0)throw std::invalid_argument("invalid HIP Euler input"); ResizeEulerTrace(t,nx);
 int fv=3*(nx+1),cv=3*nx,th=256,fb=(nx+1+th-1)/th,cb=(cv+th-1)/th; size_t fbytes=fv*sizeof(double),cbytes=cv*sizeof(double);
 double *base=nullptr,*cur=nullptr,*next=nullptr,*l=nullptr,*r=nullptr,*f=nullptr,*res=nullptr;
 try { Check(hipMalloc(&base,cbytes),"malloc base");Check(hipMalloc(&cur,cbytes),"malloc cur");Check(hipMalloc(&next,cbytes),"malloc next");Check(hipMalloc(&l,fbytes),"malloc left");Check(hipMalloc(&r,fbytes),"malloc right");Check(hipMalloc(&f,fbytes),"malloc flux");Check(hipMalloc(&res,cbytes),"malloc residual");
  Check(hipMemcpy(base,in,cbytes,hipMemcpyHostToDevice),"copy base");Check(hipMemcpy(cur,in,cbytes,hipMemcpyHostToDevice),"copy current");
  std::copy(in,in+cv,t.state.begin());
  for(int s=0;s<3;++s){hipLaunchKernelGGL(Face,dim3(fb),dim3(th),0,0,cur,nx,g,(int)b,l,r,f);Check(hipGetLastError(),"face launch");hipLaunchKernelGGL(Residual,dim3(cb),dim3(th),0,0,f,nx,1.0/dx,res);Check(hipGetLastError(),"residual launch");double bw=s==1?0.75:(s==2?1.0/3.0:0.0),uw=s==0?1.0:(s==1?0.25:2.0/3.0);hipLaunchKernelGGL(Update,dim3(cb),dim3(th),0,0,base,cur,res,cv,dt,bw,uw,next);Check(hipGetLastError(),"update launch");Check(hipDeviceSynchronize(),"sync");Check(hipMemcpy(t.faceLeft.data()+s*fv,l,fbytes,hipMemcpyDeviceToHost),"copy left");Check(hipMemcpy(t.faceRight.data()+s*fv,r,fbytes,hipMemcpyDeviceToHost),"copy right");Check(hipMemcpy(t.numericalFlux.data()+s*fv,f,fbytes,hipMemcpyDeviceToHost),"copy flux");Check(hipMemcpy(t.residual.data()+s*cv,res,cbytes,hipMemcpyDeviceToHost),"copy residual");Check(hipMemcpy(t.state.data()+(s+1)*cv,next,cbytes,hipMemcpyDeviceToHost),"copy state");std::swap(cur,next);}
 } catch(...) {(void)hipFree(res);(void)hipFree(f);(void)hipFree(r);(void)hipFree(l);(void)hipFree(next);(void)hipFree(cur);(void)hipFree(base);throw;} (void)hipFree(res);(void)hipFree(f);(void)hipFree(r);(void)hipFree(l);(void)hipFree(next);(void)hipFree(cur);(void)hipFree(base); return;
}
}
