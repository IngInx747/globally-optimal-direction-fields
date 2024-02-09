#include "mesh.hh"
#include "mesh_io.hh"
#include "ffgen.hh"
#include "ffio.hh"

using namespace OpenMesh;

static std::string mesh_filename {};
static std::string mesh_prefix {};
static std::string mesh_path {};

static const char *_var_fvec { "face:vec" };
static const char *_var_vvec { "vert:vec" };

static const char* getopt(const char **begin, const char **end, const char *option)
{
    const char **iter = std::find(begin, end, std::string { option });
    if (iter != end && ++iter != end) return *iter;
    return nullptr;
}

static bool is_opt(const char **begin, const char **end, const char *option)
{
    return std::find(begin, end, std::string { option }) != end;
}

static double get_average_edge_length(const TriMesh &mesh)
{
    double al {};
    for (auto face : mesh.faces()) al += mesh.calc_face_area(face);
    al = sqrt(al / (double)mesh.n_faces() * 4./sqrt(3.)) * .25;
    return al;
}

inline bool is_good_to_fix(const TriMesh &mesh, const Vh &vh, const double da)
{
    int nl {};
    for (auto edge : mesh.ve_range(vh))
        if (is_marked(mesh, edge)) ++nl;

    if (nl == 1) return true;

    if (nl == 2)
    {
        Vec3 de[2]; int ne {};
        for (auto edge : mesh.ve_range(vh)) if (is_marked(mesh, edge))
            de[ne++] = mesh.calc_edge_vector(edge).normalized();
        // check if the feature line is too much bended
        return abs(dot(de[0], de[1])) >= cos(da);
    }

    return false;
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

int main(const int argc, const char **argv)
{
    if (argc < 2) return 1;

    mesh_filename.append(argv[1]);
    mesh_prefix = mesh_filename.substr(0, mesh_filename.find_last_of("."));
    mesh_path = mesh_filename.substr(0, mesh_filename.find_last_of("/\\"));

    TriMesh mesh;

    if (read_mesh(mesh, mesh_filename.c_str()))
    { printf("Cannot open mesh: %s\n", mesh_filename.c_str()); return 1; }

    int n_rosy = 4;
    double rosy_s = 0.;
    double lambda = 0.;
    int alignment = 0;
    double sharp_dihedral_angle = 0.;
    double max_bending_angle = 60.0;
    const char *arg {};
    int err {};

    if ((arg = getopt(argv, argv + argc, "-n")) != nullptr) { n_rosy = atoi(arg); }
    if ((arg = getopt(argv, argv + argc, "-s")) != nullptr) { rosy_s = atof(arg); }
    if ((arg = getopt(argv, argv + argc, "-l")) != nullptr) { lambda = atof(arg); }

    if (is_opt(argv, argv + argc, "-alignment")) { alignment = 1; }
    if ((arg = getopt(argv, argv + argc, "-dihedral-angle")) != nullptr) { sharp_dihedral_angle = atof(arg); }
    if ((arg = getopt(argv, argv + argc, "-max-bending-angle")) != nullptr) { max_bending_angle = atof(arg); }

    //if (alignment && !(n_rosy == 4 || n_rosy == 2))
    //{ printf("Curvature alignment only applies for N = 2 or 4\n"); alignment = false; }

    printf("ROSY: %d\n", n_rosy);
    printf("s = %.1f\n", rosy_s);
    printf("lambda = %.1f\n", lambda);
    printf("Align? %s\n", alignment ? "YES" : "NO");

    // mark feature edges
    if (alignment)
    {
        for (auto edge : mesh.edges()) if (edge.is_boundary())
            set_marked(mesh, edge, true);

        if (sharp_dihedral_angle > 0)
        {
            printf("Sharp dihedral angle: %.1f\n", sharp_dihedral_angle);
            for (auto edge : mesh.edges())
                if (abs(mesh.calc_dihedral_angle(edge)) > radian(sharp_dihedral_angle))
                    set_marked(mesh, edge, true);
        }
    }

    // mark fixed vertices
    if (alignment)
    {
        for (auto hdge : mesh.halfedges()) if (is_marked(mesh, hdge.edge()))
            set_marked(mesh, hdge.to(), true);

        // remove bended feature vertices
        printf("Maximum bending angle: %.1f\n", max_bending_angle);
        for (auto vert : mesh.vertices()) if (is_marked(mesh, vert))
            if (!is_good_to_fix(mesh, vert, radian(max_bending_angle)))
                set_marked(mesh, vert, false);
    }

    // generate n-rosy complex
    if (alignment) err = generate_n_rosy_aligned(mesh, n_rosy, rosy_s, lambda);
    else           err = generate_n_rosy_free(mesh, n_rosy, rosy_s, lambda);
    if (err) { printf("err = %d\n", err); if (err != 2) return err; }

    int si = calculate_n_rosy_singularities(mesh, n_rosy);
    std::cout << "x_eular = " << si/(n_rosy*2) << std::endl;

//#define FF_ON_VERTEX

    // calculate one direction of the complex
#ifdef FF_ON_VERTEX
    pull_back_vertex_space(mesh, _var_vvec, n_rosy);
#else
    pull_back_face_space(mesh, _var_fvec, n_rosy);
    calculate_n_rosy_singularities(mesh, _var_fvec, n_rosy);
#endif

    double al = get_average_edge_length(mesh);
    std::stringstream ss; ss << mesh_prefix << ".nrosy" << std::to_string(n_rosy);

    // save n-rosy cross field
#ifdef FF_ON_VERTEX
    dump_vertex_n_rosy(mesh, _var_vvec, n_rosy, (ss.str() + ".obj").c_str(), al, al*1e-1);
    save_marked_face_centroids(mesh, (ss.str() + ".singularity.obj").c_str(), al*1e-1);
#else
    dump_face_n_rosy(mesh, _var_fvec, n_rosy, (ss.str() + ".obj").c_str(), al, al*1e-1);
    save_marked_vertices(mesh, (ss.str() + ".singularity.obj").c_str(), al*1e-2);
#endif

    save_face_vector(mesh, _var_fvec, ss.str().c_str());

    release_n_rosy(mesh);

    return 0;
}