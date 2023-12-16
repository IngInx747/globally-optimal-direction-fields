#include <queue>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

////////////////////////////////////////////////////////////////
/// Utility
////////////////////////////////////////////////////////////////

template <class T>
int dump_sparse_matrix(const Eigen::SparseMatrix<T> &mat, const char *path)
{
    using ColumnIterator = Eigen::SparseMatrix<T>::InnerIterator;

    std::ofstream out(path, std::ios::out);
    out << mat.rows() << " " << mat.cols() << "\n";

    for (int j = 0; j < mat.outerSize(); ++j)
        for (ColumnIterator it(mat, j); it; ++it)
            out << it.row() << " " << it.col() << " " << it.value() << "\n";

    return 0;
}

template <class T>
int dump_sparse_matrix(const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat, const char *path)
{
    using RowIterator = Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator;

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
/// General solver
////////////////////////////////////////////////////////////////

template <class T, size_t N_row, size_t N_col>
int solve_sparse_LU(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::Matrix<T, N_row, N_col> &y,
          Eigen::Matrix<T, N_row, N_col> &x)
{
    Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    x = solver.solve(y);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    return Eigen::Success;
}

template <class T, size_t N_row, size_t N_col>
int solve_simplical_LDLT(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::Matrix<T, N_row, N_col> &y,
          Eigen::Matrix<T, N_row, N_col> &x)
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    x = solver.solve(y);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    return Eigen::Success;
}

template <class T, size_t N_row, size_t N_col>
int solve_conjugate_gradient(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::Matrix<T, N_row, N_col> &y,
          Eigen::Matrix<T, N_row, N_col> &x)
{
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T>> solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    //x = solver.solve(y);
    x = solver.solveWithGuess(y, x);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    return Eigen::Success;
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
/// Minumum eigvec solver
////////////////////////////////////////////////////////////////

template <class T>
inline T rayleigh_quotient(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::SparseMatrix<T> &B,
    const Eigen::VectorX<T> &x)
{
    return x.dot(A*x) / x.dot(B*x);
};

template <class T>
inline double rayleigh_residual(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::SparseMatrix<T> &B,
    const Eigen::VectorX<T> &x)
{
    const T ev = rayleigh_quotient(A, B, x);
    return (A*x - B*x*ev).norm() / x.norm();
};

template <class T>
int solve_inversed_power(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::SparseMatrix<T> &B,
          Eigen::VectorX<T> &x,
    const double tol,
    const int n_iter)
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;
    //Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    int iter {};

    for (iter = 0; iter < n_iter; ++iter)
    { // begin of iteration

    x = solver.solve(B*x);
    if (solver.info() != Eigen::Success) { return solver.info(); }
    x.normalize();

    //printf("iter = %d\n", iter);
    //printf("res = %lf\n", rayleigh_residual(A, B, x));

    if (rayleigh_residual(A, B, x) < tol) break;

    } // end of iteration

    return iter < n_iter ? Eigen::Success : Eigen::NoConvergence;
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
