#ifndef FRAME_FIELD_IO_HH
#define FRAME_FIELD_IO_HH

#include "mesh.hh"

int read_face_vector(TriMesh&, const char*, const char*);

int save_face_vector(const TriMesh&, const char*, const char*);

int dump_face_vector  (const TriMesh&, const char*, const char*, const double scale=1, const double offset=0);

int dump_vertex_vector(const TriMesh&, const char*, const char*, const double scale=1, const double offset=0);

int dump_face_n_rosy  (const TriMesh&, const char*, const int n_rosy, const char *filename, const double scale=1, const double offset=0);

int dump_vertex_n_rosy(const TriMesh&, const char*, const int n_rosy, const char *filename, const double scale=1, const double offset=0);

#endif