#include "OneDEulerBackend.h"
#include "OneDWeno5.h"

#include "AccelRuntime.h"
#include "DeviceBuffer.h"

#include <hip/hip_runtime.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace oneflow_1d {
namespace {

// ---------------------------------------------------------------------------
// HIP execution-model helpers (stream/event/kernel — NOT abstracted by AccelBackend)
// ---------------------------------------------------------------------------
void Check(hipError_t e,const char *s){
    if(e!=hipSuccess) throw std::runtime_error(std::string(s)+": "+hipGetErrorString(e));
}

// ---------------------------------------------------------------------------
// AccelBackend memory helpers (unified with codes/accel)
// ---------------------------------------------------------------------------
ONEFLOW::AccelBackend & GetAccelBackend() {
    return ONEFLOW::AccelRuntime::Instance().Backend();
}

void CheckProblem(const EulerProblem&p){
    if(p.nx<=0||p.gamma<=1.0||p.dt<=0.0||p.dx<=0.0)
        throw std::invalid_argument("invalid Euler problem");
}

void CheckOptions(int steps,const EulerRunOptions&o){
    if(steps<=0) throw std::invalid_argument("Euler steps must be positive");
    if(o.mode==EulerRunMode::FullTrace&&(steps!=1||o.trace==nullptr))
        throw std::invalid_argument("FullTrace requires one step and a trace");
    if(o.mode==EulerRunMode::NoTrace&&o.trace!=nullptr)
        throw std::invalid_argument("NoTrace cannot provide a trace");
}

// ---------------------------------------------------------------------------
// Device kernels (unchanged — HIP-specific execution)
// ---------------------------------------------------------------------------
__device__ int CI(int c,int i,int nx){return c*nx+i;}
__device__ int FI(int c,int i,int nx){return c*(nx+1)+i;}
__device__ void PF(const double*q,int i,int nx,double g,double f[3],double&a){
    double r=q[CI(0,i,nx)],m=q[CI(1,i,nx)],e=q[CI(2,i,nx)];
    double u=m/r,p=(g-1.0)*(e-0.5*m*m/r);
    f[0]=m; f[1]=m*u+p; f[2]=(e+p)*u;
    a=fabs(u)+sqrt(g*p/r);
}
__global__ void Face(const double*q,int nx,double g,int b,double*l,double*r,double*f){
    int face=blockIdx.x*blockDim.x+threadIdx.x;
    if(face>nx)return;
    int il=face-1,ir=face;
    if(b==0){il=(face+nx-1)%nx;ir=face%nx;}
    else{il=il<0?0:il;ir=ir>=nx?nx-1:ir;}
    double fl[3],fr[3],al,ar;
    PF(q,il,nx,g,fl,al);PF(q,ir,nx,g,fr,ar);
    double a=al>ar?al:ar;
    for(int c=0;c<3;++c){
        double ql=q[CI(c,il,nx)],qr=q[CI(c,ir,nx)];
        int k=FI(c,face,nx);
        if(l!=nullptr)l[k]=ql;if(r!=nullptr)r[k]=qr;
        f[k]=0.5*(fl[c]+fr[c])-0.5*a*(qr-ql);
    }
}
__global__ void ResidualUpdate(const double*f,int nx,double id,const double*base,
    const double*cur,int n,double dt,double bw,double uw,double*res,double*next){
    int k=blockIdx.x*blockDim.x+threadIdx.x;
    if(k>=n)return;
    int c=k/nx,i=k-c*nx;
    double value=-(f[FI(c,i+1,nx)]-f[FI(c,i,nx)])*id;
    if(res!=nullptr)res[k]=value;
    next[k]=bw*base[k]+uw*(cur[k]+dt*value);
}

// ---------------------------------------------------------------------------
// HipState — persistent device state using AccelBackend + DeviceBuffer
// ---------------------------------------------------------------------------
struct HipState final: EulerState {
    explicit HipState(const EulerProblem&p):problem(p){
        CheckProblem(problem);
        int cv=3*p.nx, fv=3*(p.nx+1);
        bool weno=(problem.method==EulerMethod::Weno5);
        int pv=3*(p.nx+2*WenoGhostCells);

        try{
            Check(hipStreamCreate(&stream),"stream");
            Check(hipEventCreate(&kernelBegin),"kernel begin event");
            Check(hipEventCreate(&kernelEnd),"kernel end event");

            // Memory via AccelBackend (unified with codes/accel)
            base.Resize(cv);
            work.Resize(cv);
            scratch.Resize(cv);
            flux.Resize(fv);

            if(weno){
                padded.Resize(pv);
                positive.Resize(pv);
                negative.Resize(pv);
            }
        }catch(...){Release();throw;}
    }

    ~HipState() override{Release();}

    void Release() noexcept{
        // Synchronize before freeing
        if(stream) static_cast<void>(hipStreamSynchronize(stream));

        // DeviceBuffer auto-frees on destruction — just reset explicitly
        negative.Reset(); positive.Reset(); padded.Reset();
        residual.Reset(); flux.Reset(); right.Reset(); left.Reset();
        scratch.Reset(); work.Reset(); base.Reset();

        // Stream/event still manual (not abstracted by AccelBackend)
        if(kernelEnd)   static_cast<void>(hipEventDestroy(kernelEnd));
        if(kernelBegin) static_cast<void>(hipEventDestroy(kernelBegin));
        if(stream)      static_cast<void>(hipStreamDestroy(stream));
        stream=nullptr; kernelBegin=kernelEnd=nullptr;
    }

    void EnsureTraceBuffers(){
        int fv=3*(problem.nx+1), cv=3*problem.nx;
        try{
            if(left.Size()==0)   left.Resize(fv);
            if(right.Size()==0)  right.Resize(fv);
            if(residual.Size()==0) residual.Resize(cv);
        }catch(...){
            residual.Reset(); right.Reset(); left.Reset();
            throw;
        }
    }

    EulerProblem problem;
    hipStream_t stream=nullptr;
    hipEvent_t kernelBegin=nullptr;
    hipEvent_t kernelEnd=nullptr;

    // Device memory — now managed by DeviceBuffer (RAII, unified with codes/accel)
    ONEFLOW::DeviceBuffer<double> base, work, scratch;
    ONEFLOW::DeviceBuffer<double> left, right, flux, residual;
    ONEFLOW::DeviceBuffer<double> padded, positive, negative;  // WENO5 only

    bool uploaded=false;
};

HipState& HS(EulerState&s){
    auto*p=dynamic_cast<HipState*>(&s);
    if(!p) throw std::invalid_argument("not HIP Euler state");
    return *p;
}
const HipState& HS(const EulerState&s){
    auto*p=dynamic_cast<const HipState*>(&s);
    if(!p) throw std::invalid_argument("not HIP Euler state");
    return *p;
}

// ---------------------------------------------------------------------------
// Copy helper — unified path through AccelBackend
// ---------------------------------------------------------------------------
void CopyToDevice(const void* host, void* device, size_t bytes){
    GetAccelBackend().Copy(device, host, bytes, ONEFLOW::AccelCopyKind::HostToDevice);
}
void CopyToHost(void* host, const void* device, size_t bytes){
    GetAccelBackend().Copy(host, device, bytes, ONEFLOW::AccelCopyKind::DeviceToHost);
}

} // namespace

// ---------------------------------------------------------------------------
// HipEulerBackend — public API
// ---------------------------------------------------------------------------
const char*HipEulerBackend::Name()const{return "HIP";}
bool HipEulerBackend::IsAccelerator()const{return true;}

std::unique_ptr<EulerState>HipEulerBackend::CreateState(const EulerProblem&p)const{
    return std::make_unique<HipState>(p);
}

void HipEulerBackend::Upload(EulerState&s,const double*h)const{
    auto&x=HS(s);
    if(!h) throw std::invalid_argument("null Euler host state");
    size_t b=(size_t)(3*x.problem.nx)*sizeof(double);
    CopyToDevice(h, x.base.Data(), b);
    Check(hipStreamSynchronize(x.stream),"upload sync");
    x.uploaded=true;
}

void HipEulerBackend::Advance(EulerState&s,int steps,const EulerRunOptions&o)const{
    CheckOptions(steps,o);
    auto&x=HS(s);
    if(!x.uploaded) throw std::logic_error("HIP Euler state was not uploaded");

    const bool weno=x.problem.method==EulerMethod::Weno5;
    int nx=x.problem.nx, cv=3*nx, fv=3*(nx+1), th=256;
    int fb=(fv+th-1)/th, cb=(cv+th-1)/th;
    size_t cbytes=(size_t)cv*sizeof(double), fbytes=(size_t)fv*sizeof(double);

    if(weno){
        if(o.trace)
            throw std::invalid_argument("WENO5 FullTrace not yet supported via EulerBackend");
        if(o.stats) Check(hipEventRecord(x.kernelBegin,x.stream),"weno kernel begin");

        for(int step=0;step<steps;++step){
            HipWeno5StepRaw(
                x.base.Data(), x.base.Data(), x.work.Data(),
                x.padded.Data(), x.positive.Data(), x.negative.Data(),
                x.left.Data(), x.right.Data(), x.flux.Data(), x.residual.Data(),
                nx, x.problem.gamma, x.problem.dx,
                (int)x.problem.boundary, x.problem.dt,
                x.stream);
            std::swap(x.base, x.work);
        }

        if(o.stats){
            Check(hipEventRecord(x.kernelEnd,x.stream),"weno kernel end");
            Check(hipEventSynchronize(x.kernelEnd),"weno kernel sync");
            float ms=0.0F;
            Check(hipEventElapsedTime(&ms,x.kernelBegin,x.kernelEnd),"weno kernel elapsed");
            o.stats->kernelMilliseconds=ms;
            o.stats->kernelLaunches=steps*EulerRkStages*5;
            o.stats->synchronizationCount=1;
        }else{
            Check(hipStreamSynchronize(x.stream),"weno advance sync");
        }
        return;
    }

    /* Rusanov path */
    EulerTrace*t=o.trace;
    if(t){ResizeEulerTrace(*t,nx); x.EnsureTraceBuffers();}

    if(t){
        CopyToHost(t->state.data(), x.base.Data(), cbytes);
        Check(hipStreamSynchronize(x.stream),"initial trace sync");
    }

    if(o.stats && !t) Check(hipEventRecord(x.kernelBegin,x.stream),"kernel begin");

    for(int step=0;step<steps;++step){
        for(int stage=0;stage<3;++stage){
            double bw=stage==1?.75:(stage==2?1.0/3.0:0.0);
            double uw=stage==0?1.0:(stage==1?.25:2.0/3.0);
            double*cur=stage==0?x.base.Data():(stage==1?x.work.Data():x.scratch.Data());
            double*next=stage==0?x.work.Data():(stage==1?x.scratch.Data():x.work.Data());

            hipLaunchKernelGGL(Face,dim3(fb),dim3(th),0,x.stream,
                cur,nx,x.problem.gamma,(int)x.problem.boundary,
                t?x.left.Data():nullptr,
                t?x.right.Data():nullptr,
                x.flux.Data());
            Check(hipGetLastError(),"Face launch");

            hipLaunchKernelGGL(ResidualUpdate,dim3(cb),dim3(th),0,x.stream,
                x.flux.Data(),nx,1.0/x.problem.dx,
                x.base.Data(),cur,cv,x.problem.dt,bw,uw,
                t?x.residual.Data():nullptr,
                next);
            Check(hipGetLastError(),"ResidualUpdate launch");

            if(t){
                Check(hipStreamSynchronize(x.stream),"trace sync");
                int fo=stage*fv, co=stage*cv;
                CopyToHost(t->faceLeft.data()+fo,       x.left.Data(),     fbytes);
                CopyToHost(t->faceRight.data()+fo,      x.right.Data(),    fbytes);
                CopyToHost(t->numericalFlux.data()+fo,  x.flux.Data(),     fbytes);
                CopyToHost(t->residual.data()+co,       x.residual.Data(), cbytes);
                CopyToHost(t->state.data()+(stage+1)*cv, next,              cbytes);
            }
        }
        std::swap(x.base, x.work);
    }

    if(!t){
        if(o.stats){
            Check(hipEventRecord(x.kernelEnd,x.stream),"kernel end");
            Check(hipEventSynchronize(x.kernelEnd),"kernel sync");
            float ms=0.0F;
            Check(hipEventElapsedTime(&ms,x.kernelBegin,x.kernelEnd),"kernel elapsed");
            o.stats->kernelMilliseconds=ms;
            o.stats->kernelLaunches=steps*EulerRkStages*2;
            o.stats->synchronizationCount=1;
        }else{
            Check(hipStreamSynchronize(x.stream),"advance sync");
        }
    }
}

void HipEulerBackend::Download(const EulerState&s,double*h)const{
    const auto&x=HS(s);
    if(!x.uploaded) throw std::logic_error("HIP Euler state was not uploaded");
    if(!h) throw std::invalid_argument("null Euler host state");
    size_t b=(size_t)(3*x.problem.nx)*sizeof(double);
    CopyToHost(h, x.base.Data(), b);
}

void HipEulerBackend::Step(const double*state,int nx,double gamma,double dt,double dx,
    EulerBoundary b,EulerTrace&trace)const{
    EulerBackend::Step(state,nx,gamma,dt,dx,b,trace);
}

} // namespace oneflow_1d
