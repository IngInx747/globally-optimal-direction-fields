#include <queue>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

////////////////////////////////////////////////////////////////
/// Utility
////////////////////////////////////////////////////////////////

template <class ScalarT>
int save_matrix(const Eigen::SparseMatrix<ScalarT> &mat, const char *path)
{
    using ColumnIterator = Eigen::SparseMatrix<ScalarT>::InnerIterator;

    std::ofstream out(path, std::ios::out);
    out << mat.rows() << " " << mat.cols() << "\n";

    for (int j = 0; j < mat.outerSize(); ++j)
        for (ColumnIterator it(mat, j); it; ++it)
            out << it.row() << " " << it.col() << " " << it.value() << "\n";

    return 0;
}

template <class ScalarT>
int save_matrix(const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor> &mat, const char *path)
{
    using RowIterator = Eigen::SparseMatrix<ScalarT, Eigen::RowMajor>::InnerIterator;

    std::ofstream out(path, std::ios::out);
    out << mat.rows() << " " << mat.cols() << "\n";

    for (int i = 0; i < mat.outerSize(); ++i)
        for (RowIterator it(mat, i); it; ++it)
            out << it.row() << " " << it.col() << " " << it.value() << "\n";

    return 0;
}

template
int save_matrix(const Eigen::SparseMatrix<double>&, const char*);

template
int save_matrix(const Eigen::SparseMatrix<std::complex<double>>&, const char*);

template
int save_matrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>&, const char*);

template
int save_matrix(const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>&, const char*);

////////////////////////////////////////////////////////////////
/// General solvers
////////////////////////////////////////////////////////////////

template <class SolverT, typename ScalarT, int N_row, int N_col>
inline int solve_directly(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &y,
          Eigen::Matrix<ScalarT, N_row, N_col> &x)
{
    SolverT solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    x = solver.solve(y);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    return Eigen::Success;
}

template <class SolverT, typename ScalarT, int N_row, int N_col>
inline int solve_iteratively(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &y,
          Eigen::Matrix<ScalarT, N_row, N_col> &x)
{
    SolverT solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    x = solver.solveWithGuess(y, x);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    return Eigen::Success;
}

template <typename ScalarT, int N_row, int N_col>
inline int solve_sparse_LU(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &y,
          Eigen::Matrix<ScalarT, N_row, N_col> &x)
{
    using SolverType = Eigen::SparseLU<Eigen::SparseMatrix<ScalarT>>;
    return solve_directly<SolverType, ScalarT, N_row, N_col>(A, y, x);
}

template <typename ScalarT, int N_row, int N_col>
inline int solve_simplical_LDLT(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &y,
          Eigen::Matrix<ScalarT, N_row, N_col> &x)
{
    using SolverType = Eigen::SimplicialLDLT<Eigen::SparseMatrix<ScalarT>>;
    return solve_directly<SolverType, ScalarT, N_row, N_col>(A, y, x);
}

template <typename ScalarT, int N_row, int N_col>
inline int solve_conjugate_gradient(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &y,
          Eigen::Matrix<ScalarT, N_row, N_col> &x)
{
    using SolverType = Eigen::ConjugateGradient<Eigen::SparseMatrix<ScalarT>>;
    return solve_iteratively<SolverType, ScalarT, N_row, N_col>(A, y, x);
}

//template
//int solve_sparse_LU(
//    const Eigen::SparseMatrix<double>&,
//    const Eigen::Matrix<double, -1, 1>&,
//          Eigen::Matrix<double, -1, 1>&);

//template
//int solve_sparse_LU(
//    const Eigen::SparseMatrix<std::complex<double>>&,
//    const Eigen::Matrix<std::complex<double>, -1, 1>&,
//          Eigen::Matrix<std::complex<double>, -1, 1>&);

//template
//int solve_simplical_LDLT(
//    const Eigen::SparseMatrix<double>&,
//    const Eigen::Matrix<double, -1, 1>&,
//          Eigen::Matrix<double, -1, 1>&);

template
int solve_simplical_LDLT(
    const Eigen::SparseMatrix<std::complex<double>>&,
    const Eigen::Matrix<std::complex<double>, -1, 1>&,
          Eigen::Matrix<std::complex<double>, -1, 1>&);

//template
//int solve_conjugate_gradient(
//    const Eigen::SparseMatrix<double>&,
//    const Eigen::Matrix<double, -1, 1>&,
//          Eigen::Matrix<double, -1, 1>&);

//template
//int solve_conjugate_gradient(
//    const Eigen::SparseMatrix<std::complex<double>>&,
//    const Eigen::Matrix<std::complex<double>, -1, 1>&,
//          Eigen::Matrix<std::complex<double>, -1, 1>&);

////////////////////////////////////////////////////////////////
/// Constraint solvers
////////////////////////////////////////////////////////////////

#if 0

template <typename ScalarT>
inline int reduce_fixed_constraints(
    const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor> &A_in,
    const Eigen::VectorX<ScalarT>                       &y_in,
    const Eigen::VectorX<ScalarT>                       &x_in,
    const Eigen::VectorXi                               &is_fixed,
          Eigen::SparseMatrix<ScalarT>                  &A_out,
          Eigen::VectorX<ScalarT>                       &y_out,
          Eigen::VectorX<ScalarT>                       &x_out,
          Eigen::VectorXi                               &indices)
{
    using RowIter = Eigen::SparseMatrix<ScalarT, Eigen::RowMajor>::InnerIterator;

    indices.resize(y_in.size()); indices.setConstant(-1);
    int nv {};

    // reindex variables
    for (int i = 0; i < y_in.size(); ++i)
        if (!is_fixed(i)) indices(i) = nv++;

    // reduced linear system
    A_out.resize(nv, nv); A_out.setZero();
    y_out.resize(nv);     y_out.setZero();
    x_out.resize(nv);     x_out.setZero();

    std::vector<Eigen::Triplet<ScalarT>> coef {};

    for (int i = 0; i < A_in.outerSize(); ++i) if (!is_fixed(i))
    {
        const int id = indices(i);

        y_out(id) = y_in(i);
        x_out(id) = x_in(i);

        for (RowIter iter(A_in, i); iter; ++iter)
        {
            const int j = (int)iter.col();
            const int jd = indices(j);

            if (!is_fixed(j))
            {
                coef.emplace_back(id, jd, A_in.coeff(i, j));
            }
            else
            {
                y_out(id) -= A_in.coeff(i, j) * x_in(j);
            }
        }
    }

    A_out.setFromTriplets(coef.begin(), coef.end()); coef.clear();

    return nv;
}

template <typename ScalarT>
inline int reduce_fixed_constraints(
    const Eigen::SparseMatrix<ScalarT> &A_in,
    const Eigen::VectorX<ScalarT>      &y_in,
    const Eigen::VectorX<ScalarT>      &x_in,
    const Eigen::VectorXi              &is_fixed,
          Eigen::SparseMatrix<ScalarT> &A_out,
          Eigen::VectorX<ScalarT>      &y_out,
          Eigen::VectorX<ScalarT>      &x_out,
          Eigen::VectorXi              &indices)
{
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor> A_rm(A_in);
    return reduce_fixed_constraints(A_rm, y_in, x_in, is_fixed, A_out, y_out, x_out, indices);
}

template <class SolverT, typename ScalarT>
inline int solve_directly(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    Eigen::VectorXi I;
    Eigen::VectorX<ScalarT> x1, y1;
    Eigen::SparseMatrix<ScalarT> A1;

    reduce_fixed_constraints(A, y, x, C, A1, y1, x1, I);

    int err = solve_directly<SolverT, ScalarT, -1, 1>(A1, y1, x1);

    // write back to the original array
    for (int i = 0; i < x.size(); ++i)
        if (!C(i)) x(i) = x1(I(i));

    return err;
}

template <class SolverT, typename ScalarT>
inline int solve_iteratively(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    Eigen::VectorXi I;
    Eigen::VectorX<ScalarT> x1, y1;
    Eigen::SparseMatrix<ScalarT> A1;

    reduce_fixed_constraints(A, y, x, C, A1, y1, x1, I);

    int err = solve_iteratively<SolverT, ScalarT, -1, 1>(A1, y1, x1);

    // write back to the original array
    for (int i = 0; i < x.size(); ++i)
        if (!C(i)) x(i) = x1(I(i));

    return err;
}

#endif

template <typename StorageIndex>
inline Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, StorageIndex> column_permutation(const Eigen::VectorXi &c)
{
    Eigen::VectorXi p(c.rows());
    int nv {}; // number of variables

    for (int i = 0; i < c.rows(); ++i)
        if (c(i)) { p(i) = nv++; }

    for (int i = 0; i < c.rows(); ++i)
        if (!c(i)) { p(i) = nv++; }

    // Default Permutation ctor receives an array of row permutation.
    // Transpose it to get a column permutation.
    return Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, StorageIndex>(p).transpose();
}

/// Transform variables of the original constraint optimal
///   problem to the ones of an reduced free optimal problem.
/// x = M*u + k
///   where u is the free variable vector
///   and x is the original variable vector.
template <typename ScalarT>
inline int variables_tramsform(
    const Eigen::VectorXi              &c,
    const Eigen::VectorX<ScalarT>      &x,
          Eigen::SparseMatrix<ScalarT> &M,
          Eigen::VectorX<ScalarT>      &k)
{
    eigen_assert(c.rows() == x.rows() && "Number of variables are inconsistent");

    const int nv = (int)c.rows(); // number of variables
    const int n0 = c.sum();       // number of fixed variables
    const int n1 = nv - n0;       // number of free variables

    // Construct permutation P that
    //   [...] P = [I | 0]
    //
    // Note that when permutation P is applied to the left of the masking vector
    //   P | I | = | . | = c
    //     | 0 |   | . |
    //   it actually does the inversed permutation as of P applied on the right.
    //
    // This fact can be illustrated by transposing any matrix applied by P on the right
    //   M'^T = (MP)^T = P^T M^T
    //   in which P^T on the left permutates the same as P on the right.
    //
    // Applying P on the left of M'^T = P P^T M^T = M^T
    //   shows P on the left does an inversed permutation of P on the right.
    const auto P = column_permutation<Eigen::SparseMatrix<ScalarT>::StorageIndex>(c);

    // construct k by masking
    k = x.array() * (c.cast<ScalarT>()).array();

    // construct M
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor> M_(nv, n1); M_.setZero();
    Eigen::SparseMatrix<ScalarT> M1(n1, n1); M1.setIdentity();
    M_.bottomRows(n1) = M1;
    M = P * M_;

    return Eigen::Success;
}

/// The same as above but to initialize the reduced vector meanwhile
/// u0 = (P^T x0) | truncated | fixed part
///               | retained  | free  part
template <typename ScalarT>
inline int variables_tramsform(
    const Eigen::VectorXi              &c,
    const Eigen::VectorX<ScalarT>      &x,
          Eigen::SparseMatrix<ScalarT> &M,
          Eigen::VectorX<ScalarT>      &k,
          Eigen::VectorX<ScalarT>      &u)
{
    eigen_assert(c.rows() == x.rows() && "Number of variables are inconsistent");

    const int nv = (int)c.rows(); // number of variables
    const int n0 = c.sum();       // number of fixed variables
    const int n1 = nv - n0;       // number of free variables

    const auto P = column_permutation<Eigen::SparseMatrix<ScalarT>::StorageIndex>(c);

    // construct k by masking
    k = x.array() * (c.cast<ScalarT>()).array();

    // construct M
    Eigen::SparseMatrix<ScalarT, Eigen::RowMajor> M_(nv, n1); M_.setZero();
    Eigen::SparseMatrix<ScalarT> M1(n1, n1); M1.setIdentity();
    M_.bottomRows(n1) = M1;
    M = P * M_;

    // initialize u with free variables in x
    u = (P.transpose() * x).bottomRows(n1);

    return Eigen::Success;
}

template <class SolverT, typename ScalarT>
inline int solve_directly(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    // Variables transform
    Eigen::SparseMatrix<ScalarT> M;
    Eigen::VectorX<ScalarT> k;
    variables_tramsform(C, x, M, k);

    // Reduce constraints
    // min |A x - b|^2 => min |A (M*u + k) - b|^2 =>
    // A x = b => M^T A M u = M^T (b - Ak) => A'u = b'
    // A' ::= M^T A M
    // b' ::= M^T (b - Ak)
    Eigen::SparseMatrix<ScalarT> G;
    Eigen::VectorX<ScalarT> h;
    G = M.transpose() * A * M;
    h = M.transpose() * (y - A*k);

    // Solve reduced linear system
    Eigen::VectorX<ScalarT> u;
    int err = solve_directly<SolverT, ScalarT, -1, 1>(G, h, u);

    // Write back
    x = M*u + k;

    return err;
}

template <class SolverT, typename ScalarT>
inline int solve_iteratively(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    // Variables transform
    Eigen::SparseMatrix<ScalarT> M;
    Eigen::VectorX<ScalarT> k, u;
    variables_tramsform(C, x, M, k, u);

    // Reduce constraints
    Eigen::SparseMatrix<ScalarT> G;
    Eigen::VectorX<ScalarT> h;
    G = M.transpose() * A * M;
    h = M.transpose() * (y - A*k);

    // Solve reduced linear system
    int err = solve_iteratively<SolverT, ScalarT, -1, 1>(G, h, u);

    // Write back
    x = M*u + k;

    return err;
}

template <typename ScalarT>
inline int solve_sparse_LU(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    using SolverType = Eigen::SparseLU<Eigen::SparseMatrix<ScalarT>>;
    return solve_directly<SolverType, ScalarT>(A, C, y, x);
}

template <typename ScalarT>
inline int solve_simplical_LDLT(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    using SolverType = Eigen::SimplicialLDLT<Eigen::SparseMatrix<ScalarT>>;
    return solve_directly<SolverType, ScalarT>(A, C, y, x);
}

template <typename ScalarT>
inline int solve_conjugate_gradient(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x)
{
    using SolverType = Eigen::ConjugateGradient<Eigen::SparseMatrix<ScalarT>>;
    return solve_iteratively<SolverType, ScalarT>(A, C, y, x);
}

template
int solve_simplical_LDLT(
    const Eigen::SparseMatrix<std::complex<double>> &A,
    const Eigen::VectorXi                           &C,
    const Eigen::VectorX<std::complex<double>>      &y,
          Eigen::VectorX<std::complex<double>>      &x);

////////////////////////////////////////////////////////////////
/// EigV solver
////////////////////////////////////////////////////////////////

template <class ScalarT>
inline ScalarT rayleigh_quotient(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorX<ScalarT> &x)
{
    return x.dot(A*x) / x.dot(x);
};

template <class ScalarT>
inline double rayleigh_residual(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorX<ScalarT> &x)
{
    const ScalarT ev = rayleigh_quotient(A, x);
    return (A*x - x*ev).norm() / x.norm();
};

template <class ScalarT>
inline ScalarT rayleigh_quotient(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::SparseMatrix<ScalarT> &B,
    const Eigen::VectorX<ScalarT> &x)
{
    return x.dot(A*x) / x.dot(B*x);
};

template <class ScalarT>
inline double rayleigh_residual(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::SparseMatrix<ScalarT> &B,
    const Eigen::VectorX<ScalarT> &x)
{
    const ScalarT ev = rayleigh_quotient(A, B, x);
    return (A*x - B*x*ev).norm() / x.norm();
};

template <class SolverT, class ScalarT>
inline int solve_inversed_power(
    const Eigen::SparseMatrix<ScalarT> &A,
          Eigen::VectorX<ScalarT> &x,
    const double tol,
    const int n_iter)
{
    SolverT solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    int iter {};

    for (iter = 0; iter < n_iter; ++iter)
    { // begin of iteration

    x = x.array() - x.mean();
    x = solver.solve(x);
    if (solver.info() != Eigen::Success) { return solver.info(); }
    x.normalize();

    if (rayleigh_residual(A, x) < tol) break;

    } // end of iteration

    return iter < n_iter ? Eigen::Success : Eigen::NoConvergence;
}

template <class SolverT, class ScalarT>
inline int solve_inversed_power(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::SparseMatrix<ScalarT> &B,
          Eigen::VectorX<ScalarT> &x,
    const double tol,
    const int n_iter)
{
    SolverT solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    int iter {};

    for (iter = 0; iter < n_iter; ++iter)
    { // begin of iteration

    x = solver.solve(B*x);
    if (solver.info() != Eigen::Success) { return solver.info(); }
    x.normalize();

    if (rayleigh_residual(A, B, x) < tol) break;

    } // end of iteration

    return iter < n_iter ? Eigen::Success : Eigen::NoConvergence;
}

template <class ScalarT>
inline int solve_inversed_power(
    const Eigen::SparseMatrix<ScalarT> &A,
          Eigen::VectorX<ScalarT> &x,
    const double tol,
    const int n_iter)
{
    using SolverType = Eigen::SimplicialLDLT<Eigen::SparseMatrix<ScalarT>>;
    return solve_inversed_power<SolverType, ScalarT>(A, x, tol, n_iter);
}

template <class ScalarT>
inline int solve_inversed_power(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::SparseMatrix<ScalarT> &B,
          Eigen::VectorX<ScalarT> &x,
    const double tol,
    const int n_iter)
{
    using SolverType = Eigen::SimplicialLDLT<Eigen::SparseMatrix<ScalarT>>;
    return solve_inversed_power<SolverType, ScalarT>(A, B, x, tol, n_iter);
}

//template
//int solve_inversed_power(
//    const Eigen::SparseMatrix<double> &A,
//    const Eigen::SparseMatrix<double> &B,
//          Eigen::VectorX<double> &x,
//    const double tol,
//    const int n_iter);

template
int solve_inversed_power(
    const Eigen::SparseMatrix<std::complex<double>> &A,
    const Eigen::SparseMatrix<std::complex<double>> &B,
          Eigen::VectorX<std::complex<double>> &x,
    const double tol,
    const int n_iter);
