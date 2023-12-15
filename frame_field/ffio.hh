#ifndef FRAME_FIELD_IO_HH
#define FRAME_FIELD_IO_HH

#include "tri_mesh.hh"

int save_face_vector  (const TriMesh&, const char*, const char*, const double scale=1, const double offset=0);

int save_vertex_vector(const TriMesh&, const char*, const char*, const double scale=1, const double offset=0);

int save_face_vector  (const TriMesh&, const char*, const int n_rosy, const char *filename, const double scale=1, const double offset=0);

int save_vertex_n_rosy(const TriMesh&, const char*, const int n_rosy, const char *filename, const double scale=1, const double offset=0);

int save_selected_face_centroids(const TriMesh&, const char *filename, const double offset);

int save_selected_vertices(const TriMesh&, const char *filename, const double offset);

#endif