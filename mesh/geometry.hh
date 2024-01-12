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

#include <OpenMesh/Core/Geometry/VectorT.hh>

template <typename T, size_t N>
using VecN = OpenMesh::VectorT<T, N>;

using Vec2 = VecN<double, 2>;
using Vec3 = VecN<double, 3>;
using Vec4 = VecN<double, 4>;
using Vec5 = VecN<double, 5>;
using Vec6 = VecN<double, 6>;

using Int2 = VecN<int, 2>;
using Int3 = VecN<int, 3>;
using Int4 = VecN<int, 4>;
using Int5 = VecN<int, 5>;
using Int6 = VecN<int, 6>;

inline double cross(const Vec2 &a, const Vec2 &b)
{
    return a[0]*b[1] - a[1]*b[0];
}

inline double dihedral_angle(const Vec3 &e, const Vec3 &n0, const Vec3 &n1)
{
    return std::atan2(dot(e.normalized(), cross(n0, n1)), dot(n0, n1)); // (-pi, pi]
}

#endif