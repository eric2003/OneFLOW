#pragma once
#include <vector>

template<typename... Args>  
auto SQR(Args... args) {  
    return (... + (args * args));
}

void thomas_algorithm( const std::vector<double> & a,
    const std::vector<double> & b,
    const std::vector<double> & c,
    const std::vector<double> & d,
    std::vector<double> & x );

using RhsPtr = void(*)(int, double, std::vector<double> &, std::vector<double> &);
