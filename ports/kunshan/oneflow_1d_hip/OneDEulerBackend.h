#pragma once

#include "OneDEuler.h"
#include <memory>

namespace oneflow_1d
{

struct EulerProblem
{
    int nx = 0;
    double gamma = 1.4;
    double dt = 0.0;
    double dx = 0.0;
    EulerBoundary boundary = EulerBoundary::Periodic;
};

enum class EulerRunMode { NoTrace, FullTrace };

struct EulerRunStats
{
    double kernelMilliseconds = 0.0;
    int kernelLaunches = 0;
    int synchronizationCount = 0;
};

struct EulerRunOptions
{
    EulerRunMode mode = EulerRunMode::NoTrace;
    EulerTrace * trace = nullptr;
    EulerRunStats * stats = nullptr;
};

class EulerState
{
public:
    virtual ~EulerState() = default;
};

class EulerBackend
{
public:
    virtual ~EulerBackend() = default;
    virtual const char * Name() const = 0;
    virtual bool IsAccelerator() const = 0;
    virtual std::unique_ptr< EulerState > CreateState(
        const EulerProblem & problem ) const = 0;
    virtual void Upload(
        EulerState & state, const double * hostState ) const = 0;
    virtual void Advance(
        EulerState & state, int steps, const EulerRunOptions & options ) const = 0;
    virtual void Download(
        const EulerState & state, double * hostState ) const = 0;
    virtual void Step(
        const double * state, int nx, double gamma, double dt, double dx,
        EulerBoundary boundary, EulerTrace & trace ) const;
};

class CpuEulerBackend final : public EulerBackend
{
public:
    const char * Name() const override;
    bool IsAccelerator() const override;
    std::unique_ptr< EulerState > CreateState(
        const EulerProblem & problem ) const override;
    void Upload( EulerState & state, const double * hostState ) const override;
    void Advance(
        EulerState & state, int steps,
        const EulerRunOptions & options ) const override;
    void Download(
        const EulerState & state, double * hostState ) const override;
    void Step(
        const double * state, int nx, double gamma, double dt, double dx,
        EulerBoundary boundary, EulerTrace & trace ) const override;
};

#ifdef ONEFLOW_1D_USE_HIP
class HipEulerBackend final : public EulerBackend
{
public:
    const char * Name() const override;
    bool IsAccelerator() const override;
    std::unique_ptr< EulerState > CreateState(
        const EulerProblem & problem ) const override;
    void Upload( EulerState & state, const double * hostState ) const override;
    void Advance(
        EulerState & state, int steps,
        const EulerRunOptions & options ) const override;
    void Download(
        const EulerState & state, double * hostState ) const override;
    void Step(
        const double * state, int nx, double gamma, double dt, double dx,
        EulerBoundary boundary, EulerTrace & trace ) const override;
};
#endif

} // namespace oneflow_1d
