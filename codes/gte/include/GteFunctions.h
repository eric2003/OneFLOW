// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.2 (2017/06/18)

#pragma once

//#include <Mathematics/GteConstants.h>
#include <GteConstants.h>

#include <cmath>
#include <limits>

namespace gte
{
    // The variadic template definition is required in order to support
    // template types such as BSNumber<UInteger> and BSRational<UInteger>.
    template <typename... Types>
    class Function
    {
        // Each specialization of the template must implement the following
        // functions.  The input value can be passed by value or by reference.

        // static T ACos(T x);  // acos(x)
        // static T ACosh(T x);  // acosh(x)
        // static T ASin(T x);  // asin(x)
        // static T ASinh(T x);  // asinh(x)
        // static T ATan(T x);  // atan(x)
        // static T ATanh(T x);  // atanh(x)
        // static T ATan2(T y, T x);  // atan2(y,x)
        // static T ATanpi(T x);  // atan(x)/pi
        // static T ATan2pi(T y, T x);  // atan2(y,x)/pi
        // static T Ceil(T x);  // ceil(x)
        // static T Cos(T x);  // cos(x)
        // static T Cosh(T x);  // cosh(x)
        // static T Cospi(T x);  // cos(pi*x)
        // static T Exp(T x);  // e^x
        // static T Exp2(T x);  // 2^x
        // static T Exp10(T x);  // 10^x
        // static T FAbs(T x);  // |x|
        // static T Floor(T x);  // floor(x)
        // static T FMod(T x, T y);  // fmod(x,y)
        // static T FRExp(T x, int* exponent);  // frexp(x, exponent)
        // static T InvSqrt(T x);  // 1/sqrt(x)
        // static T LDExp(T x, int exponent);  // ldexp(x, exponent)
        // static T Log(T x);  // log_e(x)
        // static T Log2(T x);  // log_2(x)
        // static T Log10(T x);  // log_10(x)
        // static T Pow(T x, T y);  // x^y
        // static T Sin(T x);  // sin(x)
        // static T Sinh(T x);  // sinh(x)
        // static T Sinpi(T x);  // sin(pi*x)
        // static T Sqr(T x);  // x^2
        // static T Sqrt(T x);  // x^{1/2}
        // static T Tan(T x);  // tan(x)
        // static T Tanh(T x);  // tanh(x)
        // static T Sign(T x);  // sign of x as a T number
        // static int ISign(T x);  // sign of x as an integer
        // static T Clamp(T x, T min, T max);  // Clamp x to the interval [min,max].
        // static T Saturate(T x);  // Clamp x to the interval [0,1].

        // The maximum number of bisections required for root finding on an
        // interval of finite floating-point numbers.  After this number of
        // iterations, the root-bounding interval consists of two consecutive
        // floating-point numbers.  Thus, you have reached the limit of the
        // precision of the numbers.  WARNING:  The return value for BSNumber
        // or BSRational types is std::numeric_limits<unsigned int>::max(),
        // because there is no limit on precision (other than the practical
        // limitations of memory usage and computational costs).  We recommend
        // you use a reasonable number of iterations for the BS* types.
        //
        // static unsigned int GetMaxBisections();
    };

    // A specialization for 32-bit IEEE floating-point numbers.
    template <>
    class Function<float>
    {
    public:
        inline static float ACos(float x)
        {
            return acos(x);
        }

        inline static float ACosh(float x)
        {
            return log(x + sqrt(x * x - 1.0f));
        }

        inline static float ASin(float x)
        {
            return asin(x);
        }

        inline static float ASinh(float x)
        {
            return log(x + sqrt(x * x + 1.0f));
        }

        inline static float ATan(float x)
        {
            return atan(x);
        }

        inline static float ATanh(float x)
        {
            return log((1.0f + x) / (1.0f - x)) * 0.5f;
        }

        inline static float ATan2(float y, float x)
        {
            return atan2(y, x);
        }

        inline static float ATanpi(float x)
        {
            return atan(x) * (float)GTE_C_INV_PI;
        }

        inline static float ATan2pi(float y, float x)
        {
            return atan2(y, x) * (float)GTE_C_INV_PI;
        }

        inline static float Ceil(float x)
        {
            return ceil(x);
        }

        inline static float Cos(float x)
        {
            return cos(x);
        }

        inline static float Cosh(float x)
        {
            return cosh(x);
        }

        inline static float Cospi(float x)
        {
            return cos(x * (float)GTE_C_PI);
        }

        inline static float Exp(float x)
        {
            return exp(x);
        }

        inline static float Exp2(float x)
        {
            return exp(x * (float)GTE_C_LN_2);
        }

        inline static float Exp10(float x)
        {
            return exp(x * (float)GTE_C_LN_10);
        }

        inline static float FAbs(float x)
        {
            return fabs(x);
        }

        inline static float Floor(float x)
        {
            return floor(x);
        }

        inline static float FMod(float x, float y)
        {
            return fmod(x, y);
        }

        inline static float FRExp(float x, int* exponent)
        {
            return frexp(x, exponent);
        }

        inline static float InvSqrt(float x)
        {
            return 1.0f / sqrt(x);
        }

        inline static float LDExp(float x, int exponent)
        {
            return ldexp(x, exponent);
        }

        inline static float Log(float x)
        {
            return log(x);
        }

        inline static float Log2(float x)
        {
            return log(x) * (float)GTE_C_INV_LN_2;
        }

        inline static float Log10(float x)
        {
            return log10(x);
        }

        inline static float Pow(float x, float y)
        {
            return pow(x, y);
        }

        inline static float Sin(float x)
        {
            return sin(x);
        }

        inline static float Sinh(float x)
        {
            return sinh(x);
        }

        inline static float Sinpi(float x)
        {
            return sin(x * (float)GTE_C_PI);
        }

        inline static float Sqr(float x)
        {
            return x * x;
        }

        inline static float Sqrt(float x)
        {
            return sqrt(x);
        }

        inline static float Tan(float x)
        {
            return tan(x);
        }

        inline static float Tanh(float x)
        {
            return tanh(x);
        }

        inline static float Sign(float x)
        {
            return (x > 0.0f ? 1.0f : (x < 0.0f ? -1.0f : 0.0f));
        }

        inline static int ISign(float x)
        {
            return (x > 0.0f ? 1 : (x < 0.0f ? -1 : 0));
        }

        inline static float Clamp(float x, float min, float max)
        {
            return (x <= min ? min : (x >= max ? max : x));
        }

        inline static float Saturate(float x)
        {
            return (x <= 0.0f ? 0.0f : (x >= 1.0f ? 1.0f : x));
        }

        inline static void ChebyshevRatios(float t, float cosA, float& f0, float& f1)
        {
            if (cosA < 1.0f)
            {
                // The angle A is in (0,pi/2].
                float A = acos(cosA);
                float invSinA = 1.0f / sin(A);
                f0 = sin((1.0f - t) * A) * invSinA;
                f1 = sin(t * A) * invSinA;
            }
            else
            {
                // The angle theta is 0.
                f0 = 1.0f - t;
                f1 = t;
            }
        }

        inline static unsigned int GetMaxBisections()
        {
            return (unsigned int)(3 + std::numeric_limits<float>::digits - std::numeric_limits<float>::min_exponent);
        }
    };

    // A specialization for 64-bit IEEE floating-point numbers.
    template <>
    class Function<double>
    {
    public:
        inline static double ACos(double x)
        {
            return acos(x);
        }

        inline static double ACosh(double x)
        {
            return log(x + sqrt(x * x - 1.0));
        }

        inline static double ASin(double x)
        {
            return asin(x);
        }

        inline static double ASinh(double x)
        {
            return log(x + sqrt(x * x + 1.0));
        }

        inline static double ATan(double x)
        {
            return atan(x);
        }

        inline static double ATanh(double x)
        {
            return log((1.0 + x) / (1.0 - x)) * 0.5;
        }

        inline static double ATan2(double y, double x)
        {
            return atan2(y, x);
        }

        inline static double ATanpi(double x)
        {
            return atan(x) * GTE_C_INV_PI;
        }

        inline static double ATan2pi(double y, double x)
        {
            return atan2(y, x) * GTE_C_INV_PI;
        }

        inline static double Ceil(double x)
        {
            return ceil(x);
        }

        inline static double Cos(double x)
        {
            return cos(x);
        }

        inline static double Cosh(double x)
        {
            return cosh(x);
        }

        inline static double Cospi(double x)
        {
            return cos(x * GTE_C_PI);
        }

        inline static double Exp(double x)
        {
            return exp(x);
        }

        inline static double Exp2(double x)
        {
            return exp(x * GTE_C_LN_2);
        }

        inline static double Exp10(double x)
        {
            return exp(x * GTE_C_LN_10);
        }

        inline static double FAbs(double x)
        {
            return fabs(x);
        }

        inline static double Floor(double x)
        {
            return floor(x);
        }

        inline static double FMod(double x, double y)
        {
            return fmod(x, y);
        }

        inline static double FRExp(double x, int* exponent)
        {
            return frexp(x, exponent);
        }

        inline static double InvSqrt(double x)
        {
            return 1.0 / sqrt(x);
        }

        inline static double LDExp(double x, int exponent)
        {
            return ldexp(x, exponent);
        }

        inline static double Log(double x)
        {
            return log(x);
        }

        inline static double Log2(double x)
        {
            return log(x) * GTE_C_INV_LN_2;
        }

        inline static double Log10(double x)
        {
            return log10(x);
        }

        inline static double Pow(double x, double y)
        {
            return pow(x, y);
        }

        inline static double Sin(double x)
        {
            return sin(x);
        }

        inline static double Sinh(double x)
        {
            return sinh(x);
        }

        inline static double Sinpi(double x)
        {
            return sin(x * GTE_C_PI);
        }

        inline static double Sqr(double x)
        {
            return x * x;
        }

        inline static double Sqrt(double x)
        {
            return sqrt(x);
        }

        inline static double Tan(double x)
        {
            return tan(x);
        }

        inline static double Tanh(double x)
        {
            return tanh(x);
        }

        inline static double Sign(double x)
        {
            return (x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0));
        }

        inline static int ISign(double x)
        {
            return (x > 0.0 ? 1 : (x < 0.0 ? -1 : 0));
        }

        inline static double Clamp(double x, double min, double max)
        {
            return (x <= min ? min : (x >= max ? max : x));
        }

        inline static double Saturate(double x)
        {
            return (x <= 0.0 ? 0.0 : (x >= 1.0 ? 1.0 : x));
        }

        inline static unsigned int GetMaxBisections()
        {
            return (unsigned int)(3 + std::numeric_limits<double>::digits - std::numeric_limits<double>::min_exponent);
        }
    };
}
