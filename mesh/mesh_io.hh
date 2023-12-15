#ifndef MESH_IO_HH
#define MESH_IO_HH

#include "tri_mesh.hh"

int read_mesh(TriMesh&, const char*);

int save_mesh(const TriMesh&, const char*);

#endif