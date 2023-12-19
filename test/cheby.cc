#include <stdio.h>
#include <complex>
#include "chebyshev.hh"

using Comx = std::complex<double>;
using std::real;
using std::imag;
using std::conj;
using std::arg;

constexpr Comx im { 0,1 };

inline Comx e_i(const double s) { return { cos(s), sin(s) }; }

constexpr double kPi = 3.14159265358979323846264338327950288;

inline Comx _fm(double s)
{
    const Comx es = e_i(s);
    const double s2 = s*s;
    return (es*6. - 6. - im*s*6. + s2*3 + im*s2*s) / (s2*s2*3.);
}

inline Comx _f1(double s)
{
    const Comx es = e_i(s);
    const double s2 = s*s, s3 = s2*s, s4 = s2*s2, s5 = s3*s2;
    return (3. + im*s + s4/24. - im*s5/60. + es*(-3. + im*s*2. + s2*.5)) / s4;
}

inline Comx _f2(double s)
{
    const Comx es = e_i(s);
    const double s2 = s*s, s3 = s2*s, s4 = s2*s2, s5 = s3*s2;
    return (4. + im*s - im*s3/6. - s4/12. + im*s5/30. + es*(-4. + im*s*3. + s2)) / s4;
}

static double _fm_r(double s) { s = (s+1.)/2.*kPi; return real(_fm(s)); }
static double _fm_i(double s) { s = (s+1.)/2.*kPi; return imag(_fm(s)); }
static double _f1_r(double s) { s = (s+1.)/2.*kPi; return real(_f1(s)); }
static double _f1_i(double s) { s = (s+1.)/2.*kPi; return imag(_f1(s)); }
static double _f2_r(double s) { s = (s+1.)/2.*kPi; return real(_f2(s)); }
static double _f2_i(double s) { s = (s+1.)/2.*kPi; return imag(_f2(s)); }

static void eval_chebyshev_series(double (*f)(double))
{
    const int nt = 22;
    double cs[nt] {};
    const int nx = 10;

    csterm(f, nt, cs);

    printf("========\n");
    for (int i=0; i<nt; ++i)
        printf("%d: %.30e\n", i, cs[i]);

    for (int i = 1; i < nx; ++i)
    {
        double x = (2./(nx))*i - 1;
        double y0 = f(x);
        double y1 = cseval(x, nt, cs);
        printf("%e\n", abs(y0-y1));
    }
}

static void eval_chebyshev_series()
{
    eval_chebyshev_series(_fm_r);
    eval_chebyshev_series(_fm_i);
    eval_chebyshev_series(_f1_r);
    eval_chebyshev_series(_f1_i);
    eval_chebyshev_series(_f2_r);
    eval_chebyshev_series(_f2_i);
}
