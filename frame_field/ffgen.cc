#include <complex>
#include "mesh.hh"
#include "geometry.hh"
#include "linear_system.hh"

using namespace OpenMesh;

using Comx = std::complex<double>;
using std::real;
using std::imag;
using std::conj;
using std::arg;

constexpr double kPi = pi<double>();
constexpr Comx im { 0,1 };

///
/// C version of r1mach and d1mach from core (netlib)
/// specialized for IEEE arithmetic
/// 
/// MACHINE CONSTANTS (s for single, d for double)
/// {S|D}1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
/// {S|D}1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
/// {S|D}1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
/// {S|D}1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
/// {S|D}1MACH(5) = LOG10(B)
///
static const double _mach[]
{
    2.2250738585072014e-308,
    1.7976931348623157e+308,
    1.1102230246251565e-16,
    2.2204460492503131e-16,
    3.0102999566398120e-01
};

static const double _s11r[]
{
    0.0448875760891932036595562553276,
    0.0278480909574822965157922173757,
    0.00394490790249120295818107628687,
    -0.00157697939158619172562804651751,
    -0.0000886578217796691901712579357311,
    0.0000301708056772263120428135787035,
    9.521839632337438230089618156e-7,
    -3.00028307455805582080773625835e-7,
    -6.14917009583473496433650831019e-9,
    1.85133588988085286010092653662e-9,
    2.67848449041765751590373973224e-11,
    -7.82394575359355297437491915705e-12,
    -8.44240072511090922609176843848e-14,
    2.41333276776166240844516922196e-14,
    2.02015531985181413114834031833e-16,
    -5.68171271075270422851146478874e-17,
    -3.80082421064644521052871349836e-19,
    1.05551739229841670238163200361e-19,
    5.7758422925275435667221605993e-22,
    -1.58774695838716531303310462626e-22,
    -7.24181766014636685673730787292e-25
};

static const double _s11i[]
{
    0.100116671557942715638078149123,
    0.0429600096728215971268270800599,
    -0.00799014859477407505770275088389,
    -0.000664114111384495427035329182866,
    0.000240714510952202000864758517061,
    9.89085259369337382687437812294e-6,
    -3.22040860178194578481012477174e-6,
    -8.08401148192350365282200249069e-8,
    2.48351290049260966544658921605e-8,
    4.24154988067028660399867468349e-10,
    -1.25611378629704490237955971836e-10,
    -1.56053077919196502557674988724e-12,
    4.50565044006801278137904597946e-13,
    4.2641179237225098728291226479e-15,
    -1.2084245714879456268965803807e-15,
    -9.01338537885038989528688031325e-18,
    2.5180796700698002962991581923e-18,
    1.51955263898294940481729370636e-20,
    -4.19737873024216866691628952458e-21,
    -2.092488792285595339755624521e-23,
    5.72708467031136321701747126611e-24
};

static const double _s12r[]
{
    -0.376145877558191778393359413441,
    0.0775244431850198578126067647425,
    0.0120396593748540634695397747695,
    -0.00385683684390247509721340352427,
    -0.000232359275790231209370627606991,
    0.0000697318379146209092637310696007,
    2.32354473986257272021507575389e-6,
    -6.71692140309360615694979580992e-7,
    -1.43946361256617673523038166877e-8,
    4.06087820907414336567714443732e-9,
    6.10183339004616075548375321861e-11,
    -1.69196418769523832825063863136e-11,
    -1.88669746820541798989965091628e-13,
    5.16473095452962111184823547686e-14,
    4.45066881692009291504139737861e-16,
    -1.20625107617859803735741992452e-16,
    -8.28193837331508300767103116139e-19,
    2.22680015825230528892642524445e-19,
    1.24755889505424049389100515561e-21,
    -3.33254971913153176741833960484e-22,
    -1.55307002839777371508497520751e-24
};

static const double _s12i[]
{
    0.0527472790869782317601048210983,
    0.00823962722148093961886198320927,
    -0.0205185842051817330153151013327,
    -0.00184683218270819613487368071941,
    0.000569681886932212757533488372406,
    0.0000248774530818801164177266528608,
    -7.31121019876580624171992432347e-6,
    -1.92744564223806538367454388776e-7,
    5.49794278719049727550379096876e-8,
    9.78237385539447442446850072421e-10,
    -2.7341624177723508216430132999e-10,
    -3.51839815887772323640101921381e-12,
    9.68934411607055794052256859665e-13,
    9.45703963505047353201918875825e-15,
    -2.57516976113400217760868402425e-15,
    -1.97419921753098238455550504742e-17,
    5.32820017906655555903355375475e-18,
    3.29581793797656865402793252539e-20,
    -8.83137325823594007269279476114e-21,
    -4.50279718100548728336329365981e-23,
    1.19941679774924468309434420379e-23
};

static const double _m12r[]
{
    0.148523151773238914750879360089,
    -0.0117856118001224048185631301904,
    -0.00248887208039014371691400683052,
    0.000250045060357076469386198883676,
    0.0000227217776065076434637230864113,
    -2.48764935230787745662127026799e-6,
    -1.32138506847814502856384193414e-7,
    1.50966754393693942843767293542e-8,
    5.3472999553162661403204445045e-10,
    -6.26136041009708550772228055719e-11,
    -1.59574066624737000616598104732e-12,
    1.89788785691219687197167013023e-13,
    3.66030609080549274006207730375e-15,
    -4.39955659500182569051978906011e-16,
    -6.65848768159000092224193226014e-18,
    8.06343127453005031535923212263e-19,
    9.84397490339224661524630997726e-21,
    -1.19869887155210161836484730378e-21,
    -1.20634550494837590549640883469e-23,
    1.47512193662595435067359954287e-24,
    1.24549093756962710863096766634e-26
};

static const double _m12i[]
{
    -0.0454399665519585306943416687117,
    -0.0210517666740874019203591488894,
    0.00194647501081621201871675259482,
    0.000253466068123907163346571754613,
    -0.0000268083453427538717591876419304,
    -1.82138740336918117478832696004e-6,
    2.04357511048425337951376869602e-7,
    8.75944656915074206478854298947e-9,
    -1.01466837126303146739791005703e-9,
    -3.02573132377805421636557302451e-11,
    3.57358222114420372764650037191e-12,
    7.88121312149152771558608913996e-14,
    -9.42758576193708862552405242331e-15,
    -1.60439904050827900099939709069e-16,
    1.93624791035947590366500765061e-17,
    2.62394448214143482490534256935e-19,
    -3.18700789496399461681365308408e-20,
    -3.52400207248027768109209530864e-22,
    4.30074555255053206057921088056e-23,
    3.95655079023456015736315286131e-25,
    -4.84642137915095135859812028886e-26
};

/// Evaluate Chebyshev series at a point
/// Adapted from fortran algo. TOMS446
/// x  : xcoord in [-1, 1]
/// nt : number of terms in the series
/// cs : Chebyshev series generated by csterm()
inline double cseval(const double x, const int nt, const double *cs)
{
    double b0 {}, b1 {}, b2 {};
    const double x2 = x * 2.;

    for (int i = 0; i < nt; ++i)
    {
        int j = nt - i - 1;
        b2 = b1;
        b1 = b0;
        b0 = x2*b1 - b2 + cs[j];
    }

    return (b0 - b2) * 0.5;
}

/// From the original fortran inits
/// April 1977 version.  W. Fullerton, c3, Los Alamos Scientific Lab.
///
/// Initialize the orthogonal series so that inits is the number of terms
/// needed to insure the error is no larger than eta. Ordinarily, eta
/// will be chosen to be one-tenth machine precision.
inline int csinit(int nt, const double *cs, const double eta)
{
    double err { };
    int i = nt - 1;

    for ( ; err <= eta && i >= 0; --i)
        err += fabs(cs[i]);

    return i + 1;
}

static int _ns11r {};
static int _ns11i {};
static int _ns12r {};
static int _ns12i {};
static int _nm12r {};
static int _nm12i {};

static void csinits()
{
    if (!_ns11r) _ns11r = csinit(sizeof(_s11r) / sizeof(*_s11r), _s11r, _mach[2] / 10);
    if (!_ns11i) _ns11i = csinit(sizeof(_s11i) / sizeof(*_s11i), _s11i, _mach[2] / 10);
    if (!_ns12r) _ns12r = csinit(sizeof(_s12r) / sizeof(*_s12r), _s12r, _mach[2] / 10);
    if (!_ns12i) _ns12i = csinit(sizeof(_s12i) / sizeof(*_s12i), _s12i, _mach[2] / 10);
    if (!_nm12r) _nm12r = csinit(sizeof(_m12r) / sizeof(*_m12r), _m12r, _mach[2] / 10);
    if (!_nm12i) _nm12i = csinit(sizeof(_m12i) / sizeof(*_m12i), _m12i, _mach[2] / 10);
}

inline double s11r(const double t) { return cseval(t, _ns11r, _s11r); }
inline double s11i(const double t) { return cseval(t, _ns11i, _s11i); }
inline double s12r(const double t) { return cseval(t, _ns12r, _s12r); }
inline double s12i(const double t) { return cseval(t, _ns12i, _s12i); }
inline double m12r(const double t) { return cseval(t, _nm12r, _m12r); }
inline double m12i(const double t) { return cseval(t, _nm12i, _m12i); }

inline Comx s11(const double t) { return { s11r(t), s11i(t) }; }
inline Comx s12(const double t) { return { s12r(t), s12i(t) }; }
inline Comx m12(const double t) { return { m12r(t), m12i(t) }; }

inline Comx e_i(const double s) { return { cos(s), sin(s) }; }

///            k   _
///        /  / \  \
///       1  /   \  0
///     |/  /     \  \
///        /_______\
///      i   --2->   j
///
/// e_0 = e_{jk},  e_1 = e_{ki},  e_2 = e_{ij}
/// g11 = |e_1|^2, g22 = |e_2|^2, g12 = <-e_1, e_2>
inline Comx _D_JK(const double s, const double g11, const double g22, const double g12)
{
    const Comx es = e_i(s);
    const double s2 = s*s, s3 = s2*s, s4 = s2*s2, s5 = s3*s2;
    const Comx f1 = 3. + im*s            + s4/24. - im*s5/60. + es*(-3. + im*s*2. + s2*.5);
    const Comx f2 = 4. + im*s - im*s3/6. - s4/12. + im*s5/30. + es*(-4. + im*s*3. + s2);
    return (f1 * (g11 + g22) + f2 * g12) / s4;
}

///            k   _
///        /  / \  \
///       1  /   \  0
///     |/  /     \  \
///        /_______\
///      i   --2->   j
///
/// e_0 = e_{jk},  e_1 = e_{ki},  e_2 = e_{ij}
/// g11 = |e_1|^2, g22 = |e_2|^2, g12 = <-e_1, e_2>
inline Comx D_JK(const double s, const double g11, const double g22, const double g12)
{
    if (abs(s) > kPi)
    {
        return _D_JK(s, g11, g22, g12);
    }
    else if (s > 0) // [0, pi]
    {
        const double t = s*2/kPi - 1; // [-1,1]
        return conj(s11(t)*(g11 + g22) + s12(t)*g12);
    }
    else // [-pi, 0]
    {
        const double t = -s*2/kPi - 1; // [-1,1]
        return s11(t)*(g11 + g22) + s12(t)*g12;
    }
}

///            k   _
///        /  / \  \
///       1  /   \  0
///     |/  /     \  \
///        /_______\
///      i   --2->   j
///
/// e_0 = e_{jk},  e_1 = e_{ki},  e_2 = e_{ij}
/// g11 = |e_1|^2, g22 = |e_2|^2, g12 = <-e_1, e_2>
/// g00 = |e_1 + e_2|^2 = |e_1|^2 + |e_2|^2 + 2*<e_1, e_2> = g11 + g22 - 2*g12
inline Comx D_II(const double s, const double g11, const double g22, const double g12)
{
    return { ((g11 + g22 - g12*2) + s*s*(g11 + g22 + g12)/90.) / 4., 0 };
}

inline Comx _M_JK(const double s)
{
    const Comx es = e_i(s);
    const double s2 = s*s;
    return (es*6. - 6. - im*s*6. + s2*3 + im*s2*s) / (s2*s2*3.);
}

inline Comx M_JK(const double s)
{
    if (abs(s) > kPi)
    {
        return _M_JK(s);
    }
    else if (s > 0) // [0, pi]
    {
        const double t = s*2/kPi - 1; // [-1,1]
        return conj(m12(t));
    }
    else // [-pi, 0]
    {
        const double t = -s*2/kPi - 1; // [-1,1]
        return m12(t);
    }
}

inline Comx M_II()
{
    return { 1 / 6., 0 };
}

/// As the halfedge of a vertex is always the boundary one in openmesh
/// convention, we explicitly use the most clockwise halfedge instead
/// which eases the calculation of the angle around a vertex.
inline Hh anchor_halfedge(const TriMesh &mesh, const Vh &vh)
{
    const auto hh = mesh.halfedge_handle(vh);
    return mesh.is_boundary(hh) ? mesh.ccw_rotated_halfedge_handle(hh) : hh;
}

inline double calc_corner_angle(const TriMesh &mesh, const Hh &hh)
{
    auto hdge = make_smart(hh, mesh);
    if (hdge.is_boundary()) return 0;
    const double a = mesh.calc_edge_length(hdge.next());
    const double b = mesh.calc_edge_length(hdge.prev());
    const double c = mesh.calc_edge_length(hdge);
    return acos(cosine(a, b, c));
}

inline double calc_rescaled_ratio(const TriMesh &mesh, const Vh &vh)
{
    auto vert = make_smart(vh, mesh);
    if (vert.is_boundary()) return 1;
    double s {};

    for (auto hdge : vert.outgoing_halfedges())
        s += calc_corner_angle(mesh, hdge.next());

    return kPi*2. / s;
}

inline double calc_rescaled_angle(const TriMesh &mesh, const Hh &hh)
{
    double s0 {}, s1 {}; bool is_met {};
    auto vert = make_smart(hh, mesh).from();
    auto hanc = make_smart(anchor_halfedge(mesh, vert), mesh);

    auto hdge = hanc; do
    {
        const double t = calc_corner_angle(mesh, hdge.next());
        if (hdge == hh) is_met = true;
        if (is_met) s1 += t;
        else        s0 += t;
    }
    while ((hdge = hdge.prev().opp()) != hanc);

    const double s = vert.is_boundary() ? 1 : kPi*2. / (s0 + s1);
    return s0 * s;
}

inline double calc_parallel_transport(const TriMesh &mesh, const Hh &hh)
{
    auto hdge = make_smart(hh, mesh);
    const double t0 = calc_rescaled_angle(mesh, hdge);
    const double t1 = calc_rescaled_angle(mesh, hdge.opp());
    return t1 - kPi - t0;
}

inline double calc_rescaled_curvature(const TriMesh &mesh, const Fh &fh)
{
    auto face = make_smart(fh, mesh);
    double s {};

    for (auto hdge : face.halfedges())
        s += calc_corner_angle(mesh, hdge)
           * calc_rescaled_ratio(mesh, hdge.next().to());

    return s - kPi;
}

inline double normalize_angle(double t)
{
    t = fmod(t, kPi*2.);
    return (t>=kPi) ? (t - kPi*2.) : (t<-kPi) ? (t + kPi*2.) : t; // [-pi, pi)
}

static const char *_var_vso { "vert:GODF:solution" };

static int setup_indices(TriMesh &mesh)
{
    auto v_i = getOrMakeProperty<Vh, int>(mesh, var_v_index());
    int nv {};

    for (auto vert : mesh.vertices())
    { v_i[vert] = nv++; }

    return nv;
}

#if 0
static int setup_connections(TriMesh &mesh, const int n)
{
    auto e_t = getOrMakeProperty<Eh, double>(mesh, _var_etr);
    auto f_k = getOrMakeProperty<Fh, double>(mesh, _var_fkv);

    for (auto edge : mesh.edges())
    { e_t[edge] = calc_parallel_transport(mesh, edge.h0()) * n; }

    for (auto face : mesh.faces())
    { f_k[face] = calc_rescaled_curvature(mesh, face) * n; }

    return 0;
}

static int calculate_singularities_deprecated(TriMesh &mesh)
{
    auto e_t = getProperty<Eh, double>(mesh, _var_etr);
    auto f_k = getProperty<Fh, double>(mesh, _var_fkv);
    auto v_u = getProperty<Vh, Comx>(mesh, _var_vso);
    int sum_idx {};

    for (auto face : mesh.faces())
        set_marked(mesh, face, false);

    for (auto face : mesh.faces()) if (!face.is_boundary())
    {
        auto hdge0 = face.halfedge(); // h_{12}
        auto hdge1 = hdge0.next();    // h_{20}
        auto hdge2 = hdge0.prev();    // h_{01}
        const auto &z0 = v_u[hdge1.to()];
        const auto &z1 = v_u[hdge2.to()];
        const auto &z2 = v_u[hdge0.to()];
        const double a  = mesh.calc_face_area(face);
        const double t0 = sync(mesh, hdge0, e_t[hdge0.edge()]); // rho_{12}
        const double t1 = sync(mesh, hdge1, e_t[hdge1.edge()]); // rho_{20}
        const double t2 = sync(mesh, hdge2, e_t[hdge2.edge()]); // rho_{01}
        const double tf = normalize_angle(f_k[face]); // Omega
        const double w0 = normalize_angle(arg(z2/z1) - t0);
        const double w1 = normalize_angle(arg(z0/z2) - t1);
        const double w2 = normalize_angle(arg(z1/z0) - t2);
        const double ph = (w0 + w1 + w2 + tf) / (kPi*2.);
        const int idx = (int)round(ph); // index can only be -1, 0, 1
        set_marked(mesh, face, idx != 0);
        sum_idx += idx;
    }

    return sum_idx;
}
#endif

///
///            2   _
///        /  / \  \
///       1  /   \  0
///     |/  /     \  \
///        /_______\
///      0   --2->   1
///
/// e_0 = e_{12},  e_1 = e_{20},  e_2 = e_{01}
///
static int setup_mass_matrix(const TriMesh &mesh, const int n, Eigen::SparseMatrix<Comx> &M)
{
    auto v_i = getProperty<Vh, int>(mesh, var_v_index());
    const int nv = (int)mesh.n_vertices();

    std::vector<Eigen::Triplet<Comx>> coef {};

    for (auto face : mesh.faces())
    {
        auto hdge0 = face.halfedge(); // h_{12}
        auto hdge1 = hdge0.next();    // h_{20}
        auto hdge2 = hdge0.prev();    // h_{01}
        const int i0 = v_i[hdge1.to()];
        const int i1 = v_i[hdge2.to()];
        const int i2 = v_i[hdge0.to()];
        const double a  = mesh.calc_face_area(face);
        const double t0 = calc_parallel_transport(mesh, hdge0) * n; // rho_{12}
        const double t1 = calc_parallel_transport(mesh, hdge1) * n; // rho_{20}
        const double t2 = calc_parallel_transport(mesh, hdge2) * n; // rho_{01}
        const double tf = calc_rescaled_curvature(mesh, face)  * n; // Omega: bundle curvature
        const auto m_ii = M_II();
        const auto m_jk = M_JK(tf);
        const auto r0 = conj(e_i(t0)); // r_{12}^H
        const auto r1 = conj(e_i(t1)); // r_{20}^H
        const auto r2 = conj(e_i(t2)); // r_{01}^H
        coef.emplace_back(i0, i0, m_ii * a);
        coef.emplace_back(i1, i1, m_ii * a);
        coef.emplace_back(i2, i2, m_ii * a);
        coef.emplace_back(i1, i2, m_jk * r0 * a);
        coef.emplace_back(i2, i0, m_jk * r1 * a);
        coef.emplace_back(i0, i1, m_jk * r2 * a);
        coef.emplace_back(i2, i1, conj(m_jk * r0 * a));
        coef.emplace_back(i0, i2, conj(m_jk * r1 * a));
        coef.emplace_back(i1, i0, conj(m_jk * r2 * a));
    }

    M.resize(nv, nv); M.setZero();
    M.setFromTriplets(coef.begin(), coef.end()); coef.clear();

    return 0;
}

static int setup_energy_matrix(const TriMesh &mesh, const int n, Eigen::SparseMatrix<Comx> &A, const double s)
{
    auto v_i = getProperty<Vh, int>(mesh, var_v_index());
    const int nv = (int)mesh.n_vertices();

    std::vector<Eigen::Triplet<Comx>> coef {};

    for (auto face : mesh.faces())
    {
        auto hdge0 = face.halfedge(); // h_{12}
        auto hdge1 = hdge0.next();    // h_{20}
        auto hdge2 = hdge0.prev();    // h_{01}
        const int i0 = v_i[hdge1.to()];
        const int i1 = v_i[hdge2.to()];
        const int i2 = v_i[hdge0.to()];
        const double a  = mesh.calc_face_area(face);
        const double t0 = calc_parallel_transport(mesh, hdge0) * n; // rho_{12}
        const double t1 = calc_parallel_transport(mesh, hdge1) * n; // rho_{20}
        const double t2 = calc_parallel_transport(mesh, hdge2) * n; // rho_{01}
        const double tf = calc_rescaled_curvature(mesh, face)  * n; // Omega: bundle curvature
        const auto m_ii = M_II();
        const auto m_jk = M_JK(tf);
        const auto r0 = conj(e_i(t0)); // r_{12}^H
        const auto r1 = conj(e_i(t1)); // r_{20}^H
        const auto r2 = conj(e_i(t2)); // r_{01}^H
        const auto e0 = mesh.calc_edge_vector(hdge0); // e_{12}
        const auto e1 = mesh.calc_edge_vector(hdge1); // e_{20}
        const auto e2 = mesh.calc_edge_vector(hdge2); // e_{01}
        const double g00 = dot(e0, e0);
        const double g11 = dot(e1, e1);
        const double g22 = dot(e2, e2);
        const double g12 = dot(e1,-e2);
        const double g20 = dot(e2,-e0);
        const double g01 = dot(e0,-e1);
        const auto d_00 = D_II(tf, g11, g22, g12) / a;
        const auto d_11 = D_II(tf, g22, g00, g20) / a;
        const auto d_22 = D_II(tf, g00, g11, g01) / a;
        const auto d_12 = D_JK(tf, g11, g22, g12) / a;
        const auto d_20 = D_JK(tf, g22, g00, g20) / a;
        const auto d_01 = D_JK(tf, g00, g11, g01) / a;
        const auto k_ii = (m_ii * tf) * s;
        const auto k_jk = (m_jk * tf - im*.5) * s;
        coef.emplace_back(i0, i0, d_00 - k_ii);
        coef.emplace_back(i1, i1, d_11 - k_ii);
        coef.emplace_back(i2, i2, d_22 - k_ii);
        coef.emplace_back(i1, i2, (d_12 - k_jk) * r0);
        coef.emplace_back(i2, i0, (d_20 - k_jk) * r1);
        coef.emplace_back(i0, i1, (d_01 - k_jk) * r2);
        coef.emplace_back(i2, i1, conj((d_12 - k_jk) * r0));
        coef.emplace_back(i0, i2, conj((d_20 - k_jk) * r1));
        coef.emplace_back(i1, i0, conj((d_01 - k_jk) * r2));
    }

    A.resize(nv, nv); A.setZero();
    A.setFromTriplets(coef.begin(), coef.end()); coef.clear();

    return 0;
}

static int setup_curvature_alignment(const TriMesh &mesh, const int n, Eigen::VectorX<Comx> &b)
{
    auto v_i = getProperty<Vh, int>(mesh, var_v_index());
    const int nv = (int)mesh.n_vertices();

    b.resize(nv); b.setZero();

    for (auto vert : mesh.vertices())
    {
        for (auto hdge : vert.outgoing_halfedges())
        {
            const double t = calc_rescaled_angle(mesh, hdge);
            const double td = mesh.calc_dihedral_angle(hdge);
            const double l = mesh.calc_edge_length(hdge);
            b(v_i[vert]) += e_i(t*2.) * -td * l;
        }
    }

    if (n == 4) b = b.array() * b.array(); // b_i = b_i^2

    return 0;
}

inline double calc_average_rescaled_sharp_angle(const TriMesh &mesh, const Vh &vh)
{
    double avg_theta {};
    int nl {};

    for (auto hdge : mesh.voh_range(vh)) if (is_marked(mesh, hdge.edge()))
    {
        avg_theta += calc_rescaled_angle(mesh, hdge);
        ++nl;
    }

    return avg_theta / (double)nl;
}

#if 0
inline double calc_one_of_rescaled_sharp_angles(const TriMesh &mesh, const Vh &vh)
{
    auto hh = anchor_halfedge(mesh, vh);
    double max_dha {};

    for (auto hdge : mesh.voh_range(vh)) if (is_marked(mesh, hdge.edge()))
    {
        double dha = abs(mesh.calc_dihedral_angle(hdge));
        if (max_dha < dha) { max_dha = dha; hh = hdge; }
    }

    return calc_rescaled_angle(mesh, hh);
}
#endif

static int setup_fixed_boundary(const TriMesh &mesh, const int n, Eigen::VectorX<Comx> &u, Eigen::VectorXi &C)
{
    auto v_i = getProperty<Vh, int>(mesh, var_v_index());
    const int nv = (int)mesh.n_vertices();
    int nc {};

    u.resize(nv); u.setZero();
    C.resize(nv); C.setZero();

    for (auto vert : mesh.vertices())
    {
        if (is_marked(mesh, vert))
        {
            u(v_i[vert]) = e_i(calc_average_rescaled_sharp_angle(mesh, vert) * n);
            C(v_i[vert]) = 1;
            ++nc;
        }
    }

    return nc;
}

static void populate_solution(TriMesh &mesh, const Eigen::VectorX<Comx> &u)
{
    auto v_i = getProperty<Vh, int>(mesh, var_v_index());
    auto v_u = getOrMakeProperty<Vh, Comx>(mesh, _var_vso);
    for (auto vert : mesh.vertices()) v_u[vert] = u(v_i[vert]);
}

int generate_n_rosy_free(TriMesh &mesh, const int n, const double s, const double lambda)
{
    constexpr double kEps = 1e-8;
    Eigen::SparseMatrix<Comx> M, A;
    Eigen::VectorX<Comx> u;

    csinits();
    setup_indices(mesh);
    setup_mass_matrix(mesh, n, M);
    setup_energy_matrix(mesh, n, A, s);

    A += M*(-lambda + kEps);
    u = Eigen::VectorX<Comx>::Ones(A.rows());
    int err = solve_inversed_power(A, M, u, 1e-6, 8192);

    populate_solution(mesh, u);

    //for (int i = 0; i < u.size(); ++i) printf("%lf\n", degree(arg(u(i))));
    //std::cout << "min eig = " << u.dot(A*u) / u.dot(M*u) << std::endl;

    return err;
}

#if 0
int generate_n_rosy_curvature_aligned(TriMesh &mesh, const int n, const double s, const double lambda)
{
    constexpr double kEps = 1e-8;
    Eigen::SparseMatrix<Comx> M, A;
    Eigen::VectorX<Comx> u, b;

    csinits();
    setup_indices(mesh);
    setup_mass_matrix(mesh, n, M);
    setup_energy_matrix(mesh, n, A, s);

    setup_curvature_alignment(mesh, n, b);
    const double b2 = (b.conjugate().transpose()*M*b).norm();
    if (b2 > kEps) b = M*b / sqrt(b2);

    A += M*(-lambda + kEps);

    int err = solve_simplical_LDLT(A, b, u);

    populate_solution(mesh, u);

    removeProperty<Vh, int>   (mesh, _var_vid);

    return err;
}
#endif

int generate_n_rosy_aligned(TriMesh &mesh, const int n, const double s, const double lambda)
{
    constexpr double kEps = 1e-8;
    Eigen::SparseMatrix<Comx> M, A;
    Eigen::VectorX<Comx> u, b;
    Eigen::VectorXi C;
    int err {};

    csinits();
    setup_indices(mesh);
    setup_mass_matrix(mesh, n, M);
    setup_energy_matrix(mesh, n, A, s);
    A += M*(-lambda + kEps);

    // prescribed directions as fixed conditions
    int nc = setup_fixed_boundary(mesh, n, u, C);

    // curvatures as guidance of the solution
    setup_curvature_alignment(mesh, n, b);
    const double b2 = (b.conjugate().transpose()*M*b).norm();
    if (b2 > kEps) b = M*b / sqrt(b2);

    if (nc) err = solve_simplical_LDLT(A, C, b, u);
    else    err = solve_simplical_LDLT(A, b, u);

    populate_solution(mesh, u);

    return err;
}

void release_n_rosy(TriMesh &mesh)
{
    removeProperty<Vh, Comx>(mesh, _var_vso);
}

int calculate_n_rosy_singularities(TriMesh &mesh, const int n)
{
    auto v_u = getProperty<Vh, Comx>(mesh, _var_vso);
    int sum_idx {};

    for (auto face : mesh.faces())
        set_marked(mesh, face, false);

    for (auto face : mesh.faces()) if (!face.is_boundary())
    {
        auto hdge0 = face.halfedge(); // h_{12}
        auto hdge1 = hdge0.next();    // h_{20}
        auto hdge2 = hdge0.prev();    // h_{01}
        const auto &z0 = v_u[hdge1.to()];
        const auto &z1 = v_u[hdge2.to()];
        const auto &z2 = v_u[hdge0.to()];
        const double a  = mesh.calc_face_area(face);
        const double t0 = calc_parallel_transport(mesh, hdge0) * n; // rho_{12}
        const double t1 = calc_parallel_transport(mesh, hdge1) * n; // rho_{20}
        const double t2 = calc_parallel_transport(mesh, hdge2) * n; // rho_{01}
        const double tf = normalize_angle(calc_rescaled_curvature(mesh, face) * n); // Omega
        const double w0 = normalize_angle(arg(z2/z1) - t0);
        const double w1 = normalize_angle(arg(z0/z2) - t1);
        const double w2 = normalize_angle(arg(z1/z0) - t2);
        const double ph = (w0 + w1 + w2 + tf) / (kPi*2.);
        const int idx = (int)round(ph); // index can only be -1, 0, 1
        set_marked(mesh, face, idx != 0);
        sum_idx += idx;
    }

    return sum_idx;
}

//inline Vec3 pull_back(const TriMesh &mesh, const Vh &vh, const Comx &u)
//{
//    const auto ex = mesh.calc_edge_vector(mesh.halfedge_handle(vh));
//    const auto nz = mesh.calc_normal(vh).normalized();
//    const auto nx = (ex - dot(ex,nz)*nz).normalized();
//    const auto ny = cross(nz, nx).normalized();
//    return nx*real(u) + ny*imag(u);
//}

void pull_back_vertex_space(TriMesh &mesh, const char *var_vvec, const int n)
{
    auto v_u = getProperty<Vh, Comx>(mesh, _var_vso);
    auto v_v = getOrMakeProperty<Vh, Vec3>(mesh, var_vvec);

    for (auto vert : mesh.vertices())
    {
        const auto ex = mesh.calc_edge_vector(anchor_halfedge(mesh, vert));
        const auto nz = mesh.calc_normal(vert).normalized();
        const auto nx = (ex - dot(ex,nz)*nz).normalized();
        const auto ny = cross(nz, nx).normalized();
        const auto u = e_i(arg(v_u[vert]) / n);
        const auto d = nx*real(u) + ny*imag(u);
        v_v[vert] = d.normalized();
    }
}

void pull_back_face_space(TriMesh &mesh, const char *var_fvec, const int n)
{
    auto v_u = getProperty<Vh, Comx>(mesh, _var_vso);
    auto f_v = getOrMakeProperty<Fh, Vec3>(mesh, var_fvec);

    ///     2
    ///    / \
    /// 1 /   \ 0
    ///  /     \
    /// 0-------1
    ///     2
    for (auto face : mesh.faces())
    {
        auto hdge0 = face.halfedge();
        auto hdge1 = hdge0.next();
        auto hdge2 = hdge1.next();
        auto vert0 = hdge1.to();
        auto vert1 = hdge2.to();
        auto vert2 = hdge0.to();
        const double c0 = calc_corner_angle(mesh, hdge0) * n;
        const double c1 = calc_corner_angle(mesh, hdge1) * n;
        const double c2 = calc_corner_angle(mesh, hdge2) * n;
        const double t0 = calc_rescaled_angle(mesh, hdge2) * n;
        const double t1 = calc_rescaled_angle(mesh, hdge0) * n;
        const double t2 = calc_rescaled_angle(mesh, hdge1) * n;
        const auto z0 = e_i(arg(v_u[vert0]) - t0);
        const auto z1 = e_i(arg(v_u[vert1]) - t1 + c0 + c2);
        const auto z2 = e_i(arg(v_u[vert2]) - t2 - c1 - c2);
        //std::cout << z0 << ", " << z1 << ", " << z2 << "\n";
        const auto z  = (z0 + z1 + z2) / 3.;

        const auto nz = mesh.calc_normal(face).normalized();
        const auto nx = mesh.calc_edge_vector(hdge2).normalized();
        const auto ny = cross(nz, nx).normalized();
        const auto u = e_i(arg(z) / n);
        const auto d = nx*real(u) + ny*imag(u);

        f_v[face] = d.normalized();
    }
}

inline Vec2 vt_to_uv(const TriMesh &mesh, const Fh &fh, const Vec3 &vb, const Vec3 &vt)
{
    auto b0 = vb.normalized();
    auto b2 = mesh.calc_normal(fh).normalized();
    auto b1 = cross(b2, b0).normalized();
    return { dot(vt, b0), dot(vt, b1) };
}

inline int calc_mismatch(const Vec2 &u0, const Vec2 &u1, const int n)
{
    const double dt = atan2(cross(u0, u1), dot(u0, u1)); // dt = [-180, 180)
    const int m = (int)floor(dt / (kPi*2./n) + 0.5); // (dt - (-45)) / 90
    return (m + n) % n;
}

/// TODO: singularities calculation of 1-ROSY(vector) needs its own way

int calculate_n_rosy_singularities(TriMesh &mesh, const char *var_fvec, const int n)
{
    auto f_v = getProperty<Fh, Vec3>(mesh, var_fvec);
    int sum_mism_index {};

    for (auto vert : mesh.vertices())
        set_marked(mesh, vert, false);

    for (auto vert : mesh.vertices()) if (!vert.is_boundary())
    {
        int m {};

        for (auto hdge : vert.outgoing_halfedges())
        {
            auto face0 = hdge.opp().face();
            auto face1 = hdge.face();
            const auto vb = mesh.calc_edge_vector(hdge);
            const auto u0 = vt_to_uv(mesh, face0, vb, f_v[face0]);
            const auto u1 = vt_to_uv(mesh, face1, vb, f_v[face1]);
            m += calc_mismatch(u0, u1, n);
        }

        m = ((m % n) + n) % n;
        sum_mism_index += m;
        if (m) set_marked(mesh, vert, true);
    }

    return sum_mism_index;
}
