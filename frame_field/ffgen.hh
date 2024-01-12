#ifndef FRAME_FIELD_GENERATOR_HH
#define FRAME_FIELD_GENERATOR_HH

#include "tri_mesh.hh"

int generate_n_rosy_free(TriMesh&, const int n, const double s, const double lambda);

int generate_n_rosy_curvature_aligned(TriMesh&, const int n, const double s, const double lambda);

int calculate_n_rosy_singularities(TriMesh&, const int n); // singularity in face

void pull_back_vertex_space(TriMesh&, const char*, const int n);

void pull_back_face_space(TriMesh&, const char*, const int n);

int calculate_n_rosy_singularities(TriMesh&, const char*, const int n); // singularity on vertex

#endif