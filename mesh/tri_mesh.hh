#ifndef TRI_MESH_HH
#define TRI_MESH_HH

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "property.hh"

////////////////////////////////////////////////////////////////
/// Basic types
////////////////////////////////////////////////////////////////

template <typename T, size_t N>
using VecN = OpenMesh::VectorT<T, N>;

using Vec2 = VecN<double, 2>;
using Vec3 = VecN<double, 3>;
using Vec4 = VecN<double, 4>;
using Vec5 = VecN<double, 5>;
using Vec6 = VecN<double, 6>;

using Int2 = VecN<int, 2>;
using Int3 = VecN<int, 3>;

using Vh = OpenMesh::VertexHandle;
using Fh = OpenMesh::FaceHandle;
using Eh = OpenMesh::EdgeHandle;
using Hh = OpenMesh::HalfedgeHandle;

////////////////////////////////////////////////////////////////
/// Mesh
////////////////////////////////////////////////////////////////

struct TriMeshTraits : public OpenMesh::DefaultTraitsDouble
{
    // Default types
    typedef Vec2 TexCoord2D;
    typedef Vec3 TexCoord3D;

    // Default attributes
    VertexAttributes   (OpenMesh::Attributes::Status);
    FaceAttributes     (OpenMesh::Attributes::Status);
    EdgeAttributes     (OpenMesh::Attributes::Status);
    HalfedgeAttributes (OpenMesh::Attributes::Status);

    // Customized attributes
    VertexTraits   {};
    FaceTraits     {};
    EdgeTraits     {};
    HalfedgeTraits {};
};

struct TriMesh : public OpenMesh::TriMesh_ArrayKernelT<TriMeshTraits>
{};

////////////////////////////////////////////////////////////////
/// Topology
////////////////////////////////////////////////////////////////

// Tell if a halfedge is of the same orientation as its associated edge
inline bool is_sync(const TriMesh &mesh, const Hh &hh)
{
    auto hdge = make_smart(hh, mesh);
    return hdge.edge().v1() == hdge.to();
}

// Represent an edge-based value with respect of halfedge
template <typename T>
inline T sync(const TriMesh &mesh, const Hh &hh, const T &val)
{
    return is_sync(mesh, hh) ? val : -val;
}

////////////////////////////////////////////////////////////////
/// Debug
////////////////////////////////////////////////////////////////

inline void print_vert(const TriMesh &mesh, const Vh &vh)
{
    printf("(%d)", vh.idx());
}

inline void print_hdge(const TriMesh &mesh, const Hh &hh)
{
    auto hdge = make_smart(hh, mesh);
    printf("(%d, %d)", hdge.from().idx(), hdge.to().idx());
}

inline void print_edge(const TriMesh &mesh, const Eh &eh)
{
    auto edge = make_smart(eh, mesh);
    printf("(%d, %d)", edge.v0().idx(), edge.v1().idx());
}

inline void print_face(const TriMesh &mesh, const Fh &fh)
{
    auto hdge = make_smart(fh, mesh).halfedge();
    printf("(%d, %d, %d)", hdge.from().idx(), hdge.to().idx(), hdge.next().to().idx());
}

#endif