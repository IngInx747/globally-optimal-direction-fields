#ifndef CHEBYSHEV_SERIES_HH
#define CHEBYSHEV_SERIES_HH

#include <math.h>
#include <assert.h>

/// Generate Chebyshev series in [-1, 1] of degree N
/// Adapted from fortran algo. TOMS446
/// func: the function to be fit
/// nt  : number of terms
/// cs  : output series
template <class Func>
inline void csterm(const Func &func, const int nt, double *cs)
{
    constexpr double kPi = 3.14159265358979323846264338327950288;

    for (int i = 0; i < nt; ++i)
        cs[i] = 0;

    //for (int i = 0; i < nt; ++i)
    //    for (int j = 0; j < nt; ++j)
    //        cs[i] += func(cos(kPi * (j*2 + 1) / (nt*2)))
    //                * cos(kPi * i * (j*2 + 1) / (nt*2));

    for (int i = 0; i < nt; ++i)
    {
        const double y = func(cos(kPi * (i*2 + 1) / (nt*2)));
        for (int j = 0; j < nt; ++j)
            cs[j] += y * cos (kPi * j * (i*2 + 1) / (nt*2));
    }

    for (int i = 0; i < nt; ++i)
        cs[i] *= 2. / nt;
}

/// Evaluate Chebyshev series at a point
/// Adapted from fortran algo. TOMS446
/// x  : xcoord in [-1, 1]
/// nt : number of terms in the series
/// cs : Chebyshev series generated by csterm()
inline double cseval(const double x, const int nt, const double *cs)
{
    double b0 {}, b1 {}, b2 {};
    const double x2 = x * 2.;

    for (int i = 0; i < nt; ++i)
    {
        int j = nt - i - 1;
        b2 = b1;
        b1 = b0;
        b0 = x2*b1 - b2 + cs[j];
    }

    return (b0 - b2) * 0.5;
}

/// From the original fortran inits
/// April 1977 version.  W. Fullerton, c3, Los Alamos Scientific Lab.
///
/// Initialize the orthogonal series so that inits is the number of terms
/// needed to insure the error is no larger than eta. Ordinarily, eta
/// will be chosen to be one-tenth machine precision.
inline int csinit(int nt, const double *cs, const double eta)
{
    double err { };
    int i = nt - 1;

    for ( ; err <= eta && i >= 0; --i)
        err += fabs(cs[i]);

    return i + 1;
}

#endif