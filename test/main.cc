#include "tri_mesh.hh"
#include "mesh_io.hh"
#include "ffgen.hh"
#include "ffio.hh"

using namespace OpenMesh;

static std::string mesh_filename {};
static std::string mesh_prefix {};
static std::string mesh_path {};

static const char *var_fvecs[] { "face:vec0", "face:vec1", "face:vec2", "face:vec3", "face:vec4", "face:vec5" };
static const char *var_vvecs[] { "vert:vec0", "vert:vec1", "vert:vec2", "vert:vec3", "vert:vec4", "vert:vec5" };

static double get_average_edge_length(const TriMesh &mesh)
{
    double al {};
    for (auto face : mesh.faces()) al += mesh.calc_face_area(face);
    al = sqrt(al / (double)mesh.n_faces() * 4./sqrt(3.)) * .25;
    return al;
}

static void generate_random_cross_field(TriMesh &mesh, const char *var_vec0, const char *var_vec1)
{
    auto f_v0 = getOrMakeProperty<Fh, Vec3>(mesh, var_vec0);
    auto f_v1 = getOrMakeProperty<Fh, Vec3>(mesh, var_vec1);

    for (auto face : mesh.faces())
    {
        const auto d = mesh.calc_edge_vector(face.halfedge());
        const auto n = mesh.calc_normal(face);
        f_v0[face] = d.normalized();
        f_v1[face] = normalize(cross(n, d));
    }
}

static void test_random_cross_field(TriMesh &mesh)
{
    generate_random_cross_field(mesh, var_fvecs[0], var_fvecs[1]);
    double al = get_average_edge_length(mesh);
    save_face_vector(mesh, var_fvecs[0], (mesh_prefix + ".vec0.obj").c_str(), al, al*1e-1);
    save_face_vector(mesh, var_fvecs[1], (mesh_prefix + ".vec1.obj").c_str(), al, al*1e-1);
}

static void test_global_optimal_cross_field(TriMesh &mesh)
{
    const int n = 4;
    const double s = 0.;
    const double lambda = 0.;
    // generate n-rosy complex
    //int err = generate_n_rosy_free(mesh, n, s, lambda);
    int err = generate_n_rosy_curvature_aligned(mesh, n, s, lambda);
    if (err) { printf("err = %d\n", err); if (err != 2) return; }
    // calculate one direction of the complex
    pull_back_vertex_space(mesh, var_vvecs[0], n);
    // save n-rosy cross field
    double al = get_average_edge_length(mesh);
    std::stringstream ss; ss << mesh_prefix << ".nrosy" << std::to_string(n);
    save_vertex_n_rosy(mesh, var_vvecs[0], n, (ss.str() + ".obj").c_str(), al, al*1e-1);
    save_selected_face_centroids(mesh, (ss.str() + ".singularity.obj").c_str(), al*1e-1);
    //save_vertex_vector(mesh, var_vvecs[0], (mesh_prefix + ".vec0.obj").c_str(), al, al*1e-1);
}

int main(const int argc, const char **argv)
{
    if (argc < 2) return 1;

    mesh_filename.append(argv[1]);
    mesh_prefix = mesh_filename.substr(0, mesh_filename.find_last_of("."));
    mesh_path = mesh_filename.substr(0, mesh_filename.find_last_of("/\\"));

    TriMesh mesh;

    if (read_mesh(mesh, mesh_filename.c_str()))
    { printf("Cannot open mesh: %s\n", mesh_filename.c_str()); return 1; }

    //test_random_cross_field(mesh);
    test_global_optimal_cross_field(mesh);

    return 0;
}