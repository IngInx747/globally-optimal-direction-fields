#ifndef GEOMETRY_HH
#define GEOMETRY_HH

template <typename T> constexpr T pi()
{
    return (T)(3.14159265358979323846264338327950288);
}

template <typename T> inline int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> inline T round(T val, T mtp)
{
    return mtp * std::round(val / mtp);
}

template <typename T> inline T radian(T deg)
{
    return (T)(deg * pi<T>() / (T)(180.));
}

template <typename T> inline T degree(T rad)
{
    return (T)(rad * (T)(180.) / pi<T>());
}

inline double cosine(double a, double b, double c)
{
    const double cs = (a*a + b*b - c*c) / (a*b*2);
    return cs < -1 ? -1 : cs > 1 ? 1 : cs;
}

#endif