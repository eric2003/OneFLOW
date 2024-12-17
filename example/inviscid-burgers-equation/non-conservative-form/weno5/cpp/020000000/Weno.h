#pragma once
#include <vector>

template<typename... Args>
auto sum(Args... args) {
    return (... + args);
}

template<typename... Args>  
auto SQR(Args... args) {  
    return (... + (args * args)); // 使用折叠表达式计算平方和  
}

double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );

void wenoL( int N, std::vector<double> & u, std::vector<double> & f );
void wenoR( int N, std::vector<double> & u, std::vector<double> & f );
void WenoRhs( int N, double dx, std::vector<double> & u, std::vector<double> & r );

