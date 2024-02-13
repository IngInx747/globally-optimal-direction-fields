#include <queue>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

////////////////////////////////////////////////////////////////
/// Utility
////////////////////////////////////////////////////////////////

template <class ScalarT>
int dump_sparse_matrix(const Eigen::SparseMatrix<ScalarT> &mat, const char *path)
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
int dump_sparse_matrix(const Eigen::SparseMatrix<ScalarT, Eigen::RowMajor> &mat, const char *path)
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
int dump_sparse_matrix(const Eigen::SparseMatrix<double>&, const char*);

template
int dump_sparse_matrix(const Eigen::SparseMatrix<std::complex<double>>&, const char*);

template
int dump_sparse_matrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>&, const char*);

template
int dump_sparse_matrix(const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>&, const char*);

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