#ifndef MESH_HH
#define MESH_HH

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Casts.hh>
#include "property.hh"
#include "geometry.hh"

using Hh = OpenMesh::HalfedgeHandle;
using Vh = OpenMesh::VertexHandle;
using Fh = OpenMesh::FaceHandle;
using Eh = OpenMesh::EdgeHandle;

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

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Fh &fh) { return mesh.status(fh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Fh &fh, const bool val) { mesh.status(fh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Hh &hh) { return mesh.status(hh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Hh &hh, const bool val) { mesh.status(hh).set_selected(val); }

template <class MeshT>
inline bool is_sharp(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).feature(); }
template <class MeshT>
inline void set_sharp(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_feature(val); }

template <class MeshT>
inline bool is_sharp(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).feature(); }
template <class MeshT>
inline void set_sharp(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_feature(val); }

template <class MeshT>
inline bool is_fixed(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).locked(); }
template <class MeshT>
inline void set_fixed(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_locked(val); }

template <class MeshT>
inline bool is_fixed(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).locked(); }
template <class MeshT>
inline void set_fixed(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_locked(val); }

//struct PolyMesh : public OpenMesh::PolyMesh_ArrayKernelT<MeshTraits> {};

using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<MeshTraits>;

//struct TriMesh : public OpenMesh::TriMesh_ArrayKernelT<MeshTraits> {};

using TriMesh = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;

inline const char *var_v_index() { return "vert:index"; }
inline const char *var_f_index() { return "face:index"; }
inline const char *var_e_index() { return "edge:index"; }
inline const char *var_h_index() { return "hdge:index"; }

#endif