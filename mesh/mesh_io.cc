#include <string>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "tri_mesh.hh"

using namespace OpenMesh;

enum IO_ERR_TYPE
{
    NO_ERROR,
    INTERNAL,
    CANNOT_OPEN,
    UNSUPPORTED_FORMAT,
};

static const char *__io_err_msg[] = {
    "",
    "internal error",
    "cannot open file",
    "unsupported format"
};

inline std::string get_file_extension(const char *filename)
{
    std::string s { filename };
    auto found = s.find_last_of(".");
    return s.substr(found + 1);
}

inline bool is_file_extension(const char *filename, const char *extension)
{
    return get_file_extension(filename).compare(extension) == 0;
}

static int read_mesh_detri2(TriMesh &mesh, const char *filename)
{
    constexpr size_t kInf = std::numeric_limits<std::streamsize>::max();

    std::ifstream in(filename, std::ios::in);
    if (!in) return IO_ERR_TYPE::CANNOT_OPEN;

    std::vector<std::tuple<int, int>> edge_list;
    int ver {}, dim {}, nv {}, nf {};

    while (in)
    {
        std::string buf {};
        in >> buf;

        if (buf.compare("MeshVersionFormatted") == 0) // version
        {
            in >> ver;
        }
        else if (buf.compare("Dimension") == 0) // dimension
        {
            in >> dim;
        }
        else if (buf.compare("Vertices") == 0) // vertices
        {
            in >> nv;
            for (int i = 0; i < nv; ++i)
            {
                TriMesh::Point p { 0,0,0 };
                for (int j = 0; j < dim; ++j)
                    in >> p[j];
                mesh.new_vertex(p);
                // optional: process vertex tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("Triangles") == 0) // faces
        {
            in >> nf;
            in.ignore(kInf, '\n');
            for (int i = 0; i < nf; ++i)
            {
                int ids[3];
                for (int j = 0; j < 3; ++j)
                {
                    in >> ids[j];
                    ids[j] -= 1;
                }
                std::vector<VertexHandle> vhs {
                    mesh.vertex_handle(ids[0]),
                    mesh.vertex_handle(ids[1]),
                    mesh.vertex_handle(ids[2])
                };
                mesh.add_face(vhs);
                // optional: process facet tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("Edges") == 0) // edge
        {
            int n {};
            in >> n;
            for (int i = 0; i < n; ++i)
            {
                int ids[2];
                in >> ids[0] >> ids[1];
                edge_list.emplace_back(ids[0] - 1, ids[1] - 1);
                // optional: process edge tag
                in.ignore(kInf, '\n');
            }
        }
        else if (buf.compare("End") == 0)
        {}
    }

    for (const auto &vid : edge_list)
    {
        int i0 = std::get<0>(vid);
        int i1 = std::get<1>(vid);
        auto hdge = mesh.find_halfedge(mesh.vertex_handle(i0), mesh.vertex_handle(i1));
        if (hdge.is_valid()) mesh.status(hdge.edge()).set_selected(true);
    }

    return IO_ERR_TYPE::NO_ERROR;
}

static int save_mesh_detri2(const TriMesh &mesh, const char *filename)
{
    std::ofstream out(filename, std::ios::out);
    if (!out) return IO_ERR_TYPE::CANNOT_OPEN;

    out << std::fixed << std::setprecision(17);

    out << "MeshVersionFormatted 1\n\nDimension\n3\n\n";

    out << "Vertices\n" << mesh.n_vertices() << "\n";
    for (auto vert : mesh.vertices())
    {
        const auto &p = mesh.point(vert);
        out << p[0] << "  " << p[1] << "  " << p[2] << "\n";
    }

    out << "Triangles\n" << mesh.n_faces() << "\n";
    for (auto face : mesh.faces())
    {
        for (auto vert : face.vertices())
            out << vert.idx() + 1 << " ";
        out << "1\n"; // arbitrary positive tag
    }

    int ne {};
    for (auto edge : mesh.edges()) if (edge.selected()) ++ne;

    out << "Edges\n" << ne << "\n";
    for (auto edge : mesh.edges())
    {
        if (edge.selected())
        {
            auto vert0 = edge.v0(), vert1 = edge.v1();
            out << vert0.idx() + 1 << " " << vert1.idx() + 1 << " -1\n";
        }
    }

    out << "End\n";

    return IO_ERR_TYPE::NO_ERROR;
}

int read_mesh(TriMesh &mesh, const char *filename)
{
    int err { IO_ERR_TYPE::NO_ERROR };
    auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0 || ext.compare("off") == 0)
    {
        IO::Options opt;

        opt += IO::Options::VertexTexCoord;
        mesh.request_vertex_texcoords2D();

        opt += IO::Options::FaceColor;
        mesh.request_face_colors();

        if (!IO::read_mesh(mesh, filename, opt)) err = IO_ERR_TYPE::INTERNAL;

        if (!opt.vertex_has_texcoord()) mesh.release_vertex_texcoords2D();
        if (!opt.face_has_color()) mesh.release_face_colors();
    }
    else if (ext.compare("mesh") == 0)
    {
        err = read_mesh_detri2(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err != 0) fprintf(stderr, "read_mesh error %d: %s\n", err, __io_err_msg[err]);
    return err;
}

int save_mesh(const TriMesh &mesh, const char *filename)
{
    int err {};
    auto ext = get_file_extension(filename);

    if (ext.compare("obj") == 0 || ext.compare("off") == 0)
    {
        IO::Options opt;

        if (mesh.has_face_colors()) opt += IO::Options::FaceColor;
        //if (mesh.has_vertex_colors()) opt += IO::Options::VertexColor;
        if (mesh.has_vertex_texcoords2D()) opt += IO::Options::VertexTexCoord;

        if (!IO::write_mesh(mesh, filename, opt, 17i64)) err = IO_ERR_TYPE::INTERNAL;
    }
    else if (ext.compare("mesh") == 0)
    {
        err = save_mesh_detri2(mesh, filename);
    }
    else
    {
        err = IO_ERR_TYPE::UNSUPPORTED_FORMAT;
    }

    if (err != 0) fprintf(stderr, "save_mesh error %d: %s\n", err, __io_err_msg[err]);
    return err;
}
