#ifndef MESH_TRAITS_HH
#define MESH_TRAITS_HH

#include <OpenMesh/Core/Mesh/Traits.hh>
#include "geometry.hh"

struct MeshTraits : public OpenMesh::DefaultTraitsDouble
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

#endif