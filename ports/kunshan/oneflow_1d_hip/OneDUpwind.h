#pragma once

void OneDHipUpwindResidual(
    const double * u, int uSize, int uStart, int nx,
    double c, double dx, double * uLeft, double * uRight, double * residual );

void OneDCpuUpwindResidual(
    const double * u, int uSize, int uStart, int nx,
    double c, double dx, double * uLeft, double * uRight, double * residual );
