#ifndef FRAME_FIELD_GENERATOR_HH
#define FRAME_FIELD_GENERATOR_HH

#include "tri_mesh.hh"

int generate_n_rosy_free(TriMesh&, const int n, const double s, const double lambda);

int generate_n_rosy_curvature_aligned(TriMesh&, const int n, const double s, const double lambda);

void pull_back_vertex_space(TriMesh&, const char*, const int n);

#endif