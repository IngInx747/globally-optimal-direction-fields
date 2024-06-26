#include <fstream>
#include <iomanip> // std::setprecision
#include "mesh.hh"
#include "geometry.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// Save data
////////////////////////////////////////////////////////////////

int read_face_vector(TriMesh &mesh, const char *var_vec, const char *filename)
{
    std::ifstream in(filename, std::ios::in);
    if (!in) return in.bad();

    auto f_v = getOrMakeProperty<Fh, Vec3>(mesh, var_vec);

    for (std::string line; std::getline(in, line); )
    {
        std::stringstream ss(line);
        int id {};
        Vec3 v {};

        ss >> id;
        auto fh = mesh.face_handle(id);

        for (int i = 0; i < 3; ++i)
            ss >> v[i];

        if (fh.is_valid()) f_v[fh] = v;
    }

    return 0;
}

int save_face_vector(const TriMesh &mesh, const char *var_vec, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::fixed << std::setprecision(17);

    auto f_v = getProperty<Fh, Vec3>(mesh, var_vec);

    for (auto face : mesh.faces())
    {
        const int id = face.idx();
        const auto &v = f_v[face];
        out << id   << " "
            << v[0] << " "
            << v[1] << " "
            << v[2] << " "
            << "\n";
    }

    return 0;
}

////////////////////////////////////////////////////////////////
/// Dump for visualization
////////////////////////////////////////////////////////////////

inline Vec3 rotate(const Vec3 &d, const Vec3 &n, const double t)
{
    const auto nz = n.normalized();
    const double pr = dot(d, nz);
    const auto nx = (d - nz*pr).normalized();
    const auto ny = cross(nz, nx).normalized();
    const double cs = cos(t), ss = sin(t);
    return nx*cs + ny*ss + nz*pr;
}

static int save_obj(
    const double *vs, const int nv,
    const int    *fs, const int nf,
    const int    *es, const int ne,
    const char   *filename,
    const std::streamsize prec = 17i64)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return out.bad();

    out << std::fixed << std::setprecision(prec);

    for (int i = 0; i < nv; ++i)
        out << "v "
            << vs[i*3 + 0] << " "
            << vs[i*3 + 1] << " "
            << vs[i*3 + 2] << " "
            << "\n";

    for (int i = 0; i < nf; ++i)
        out << "f "
            << fs[i*3 + 0] << " "
            << fs[i*3 + 1] << " "
            << fs[i*3 + 2] << " "
            << "\n";

    if (es)
    for (int i = 0; i < ne; ++i)
        out << "l "
            << es[i*2 + 0] << " "
            << es[i*2 + 1] << " "
            << "\n";

    return 0;
}

int dump_face_vector(const TriMesh &mesh, const char *var_vec, const char *filename, const double scale, const double offset)
{
    if (!hasProperty<Fh, Vec3>(mesh, var_vec)) return 1;
    auto f_v = getProperty<Fh, Vec3>(mesh, var_vec);

    std::vector<Vec3> ps {};
    std::vector<Int2> es {};
    int np {};

    for (const auto face : mesh.faces())
    {
        const auto p = mesh.calc_face_centroid(face);
        const auto n = mesh.calc_normal(face) * offset;
        const auto d = f_v[face] * scale;
        ps.push_back(p + n);
        ps.push_back(p + n + d);
        es.push_back({ np+1, np+2 });
        np += 2;
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        (const int*)es.data(),    (const int)es.size(),
        filename);
}

int dump_vertex_vector(const TriMesh &mesh, const char *var_vec, const char *filename, const double scale, const double offset)
{
    if (!hasProperty<Vh, Vec3>(mesh, var_vec)) return 1;
    auto v_v = getProperty<Vh, Vec3>(mesh, var_vec);

    std::vector<Vec3> ps {};
    std::vector<Int2> es {};
    int np {};

    for (const auto vert : mesh.vertices())
    {
        const auto p = mesh.point(vert);
        const auto n = mesh.calc_normal(vert) * offset;
        const auto d = v_v[vert] * scale;
        ps.push_back(p + n);
        ps.push_back(p + n + d);
        es.push_back({ np+1, np+2 });
        np += 2;
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        (const int*)es.data(),    (const int)es.size(),
        filename);
}

static int dump_face_n_rosy_odd(const TriMesh &mesh, const char *var_vec, const int n_rosy, const char *filename, const double scale, const double offset)
{
    if (!hasProperty<Fh, Vec3>(mesh, var_vec)) return 1;
    auto f_v = getProperty<Fh, Vec3>(mesh, var_vec);

    const double dt = pi<double>()*2. / n_rosy;

    std::vector<Vec3> ps {};
    std::vector<Int2> es {};
    int np {};

    for (int r = 0; r < n_rosy; ++r)
    {
        const double t = dt * r;

        for (const auto face : mesh.faces())
        {
            const auto p = mesh.calc_face_centroid(face);
            const auto n = mesh.calc_normal(face) * offset;
            const auto d = rotate(f_v[face], n, t) * scale;
            ps.push_back(p + n);
            ps.push_back(p + n + d);
            es.push_back({ np+1, np+2 });
            np += 2;
        }
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        (const int*)es.data(),    (const int)es.size(),
        filename);
}

static int dump_face_n_rosy_even(const TriMesh &mesh, const char *var_vec, const int n_rosy, const char *filename, const double scale, const double offset)
{
    if (!hasProperty<Fh, Vec3>(mesh, var_vec)) return 1;
    auto f_v = getProperty<Fh, Vec3>(mesh, var_vec);

    const double dt = pi<double>()*2. / n_rosy;

    std::vector<Vec3> ps {};
    std::vector<Int2> es {};
    int np {};

    for (int r = 0; r < n_rosy / 2; ++r)
    {
        const double t = dt * r;

        for (const auto face : mesh.faces())
        {
            const auto p = mesh.calc_face_centroid(face);
            const auto n = mesh.calc_normal(face) * offset;
            const auto d = rotate(f_v[face], n, t) * scale;
            ps.push_back(p + n - d);
            ps.push_back(p + n + d);
            es.push_back({ np+1, np+2 });
            np += 2;
        }
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        (const int*)es.data(),    (const int)es.size(),
        filename);
}

int dump_face_n_rosy(const TriMesh &mesh, const char *var_vec, const int n_rosy, const char *filename, const double scale, const double offset)
{
    if (n_rosy % 2) // odd
        return dump_face_n_rosy_odd (mesh, var_vec, n_rosy, filename, scale, offset);
    else // even
        return dump_face_n_rosy_even(mesh, var_vec, n_rosy, filename, scale, offset);
}

static int dump_vertex_n_rosy_odd(const TriMesh &mesh, const char *var_vec, const int n_rosy, const char *filename, const double scale, const double offset)
{
    if (!hasProperty<Vh, Vec3>(mesh, var_vec)) return 1;
    auto v_v = getProperty<Vh, Vec3>(mesh, var_vec);

    const double dt = pi<double>()*2. / n_rosy;

    std::vector<Vec3> ps {};
    std::vector<Int2> es {};
    int np {};

    for (int r = 0; r < n_rosy; ++r)
    {
        const double t = dt * r;

        for (const auto vert : mesh.vertices())
        {
            const auto p = mesh.point(vert);
            const auto n = mesh.calc_normal(vert) * offset;
            const auto d = rotate(v_v[vert], n, t) * scale;
            ps.push_back(p + n);
            ps.push_back(p + n + d);
            es.push_back({ np+1, np+2 });
            np += 2;
        }
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        (const int*)es.data(),    (const int)es.size(),
        filename);
}

static int dump_vertex_n_rosy_even(const TriMesh &mesh, const char *var_vec, const int n_rosy, const char *filename, const double scale, const double offset)
{
    if (!hasProperty<Vh, Vec3>(mesh, var_vec)) return 1;
    auto v_v = getProperty<Vh, Vec3>(mesh, var_vec);

    const double dt = pi<double>()*2. / n_rosy;

    std::vector<Vec3> ps {};
    std::vector<Int2> es {};
    int np {};

    for (int r = 0; r < n_rosy / 2; ++r)
    {
        const double t = dt * r;

        for (const auto vert : mesh.vertices())
        {
            const auto p = mesh.point(vert);
            const auto n = mesh.calc_normal(vert) * offset;
            const auto d = rotate(v_v[vert], n, t) * scale;
            ps.push_back(p + n - d);
            ps.push_back(p + n + d);
            es.push_back({ np+1, np+2 });
            np += 2;
        }
    }

    return save_obj(
        (const double*)ps.data(), (const int)ps.size(),
        nullptr,                  0,
        (const int*)es.data(),    (const int)es.size(),
        filename);
}

int dump_vertex_n_rosy(const TriMesh &mesh, const char *var_vec, const int n_rosy, const char *filename, const double scale, const double offset)
{
    if (n_rosy % 2) // odd
        return dump_vertex_n_rosy_odd (mesh, var_vec, n_rosy, filename, scale, offset);
    else // even
        return dump_vertex_n_rosy_even(mesh, var_vec, n_rosy, filename, scale, offset);
}
