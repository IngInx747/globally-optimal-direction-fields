#ifndef MESH_IO_HH
#define MESH_IO_HH

template <class Mesh>
int read_mesh(Mesh&, const char*);

template <class Mesh>
int save_mesh(const Mesh&, const char*);

#endif