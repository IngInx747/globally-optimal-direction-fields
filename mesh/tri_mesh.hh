#ifndef TRI_MESH_HH
#define TRI_MESH_HH

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "handle.hh"
#include "property.hh"
#include "geometry.hh"

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

inline void print_handle(const TriMesh &mesh, const Vh &vh, const int offset = 0)
{
    printf("(%d)", vh.idx() + offset);
}

inline void print_handle(const TriMesh &mesh, const Hh &hh, const int offset = 0)
{
    printf("(%d, %d)",
        mesh.from_vertex_handle(hh).idx() + offset,
        mesh.to_vertex_handle(hh).idx()   + offset);
}

inline void print_handle(const TriMesh &mesh, const Eh &eh, const int offset = 0)
{
    const auto hh = mesh.halfedge_handle(eh, 0);
    printf("(%d, %d)",
        mesh.from_vertex_handle(hh).idx() + offset,
        mesh.to_vertex_handle(hh).idx()   + offset);
}

inline void print_handle(const TriMesh &mesh, const Fh &fh, const int offset = 0)
{
    const auto hh = mesh.halfedge_handle(fh);
    printf("(%d, %d, %d)",
        mesh.from_vertex_handle(hh).idx() + offset,
        mesh.to_vertex_handle(hh).idx()   + offset,
        mesh.to_vertex_handle(mesh.next_halfedge_handle(hh)).idx() + offset);
}

#endif