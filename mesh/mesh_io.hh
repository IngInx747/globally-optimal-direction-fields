#ifndef MESH_IO_HH
#define MESH_IO_HH

#include "mesh.hh"

template <class MeshT>
int read_mesh(MeshT&, const char*);

template <class MeshT>
int save_mesh(const MeshT&, const char*);

template <class MeshT>
int save_marked_face_centroids(const MeshT&, const char*, const double offset);

template <class MeshT>
int save_marked_vertices(const MeshT&, const char*, const double offset);

int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char   *filename,
    const std::streamsize prec = 17i64);

#endif