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
    const Eigen::VectorX<T> &x)
{
    return x.dot(A*x) / x.dot(x);
};

template <class T>
inline double rayleigh_residual(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::VectorX<T> &x)
{
    const T ev = rayleigh_quotient(A, x);
    return (A*x - x*ev).norm() / x.norm();
};

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

    x = x.array() - x.mean();
    x = solver.solve(x);
    if (solver.info() != Eigen::Success) { return solver.info(); }
    x.normalize();

    if (rayleigh_residual(A, x) < tol) break;

    } // end of iteration

    return iter < n_iter ? Eigen::Success : Eigen::NoConvergence;
}

template <class T>
int solve_inversed_power(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::SparseMatrix<T> &B,
          Eigen::VectorX<T> &x,
    const double tol,
    const int n_iter)
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>> solver;

    solver.compute(A);
    if (solver.info() != Eigen::Success) { return solver.info(); }

    int iter {};

    for (iter = 0; iter < n_iter; ++iter)
    { // begin of iteration

    x = solver.solve(B*x);
    if (solver.info() != Eigen::Success) { return solver.info(); }
    x.normalize();

    //printf("iter[%d], res = %lf\n", iter, rayleigh_residual(A, B, x));

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

////////////////////////////////////////////////////////////////
/// Constraint
////////////////////////////////////////////////////////////////

template <class T>
inline int reduce_fixed_constraints(
    const Eigen::SparseMatrix<T, Eigen::RowMajor> &A_in,
    const Eigen::VectorX<T>                       &y_in,
    const Eigen::VectorX<T>                       &x_in,
    const Eigen::VectorXi                         &is_fixed,
          Eigen::SparseMatrix<T>                  &A_out,
          Eigen::VectorX<T>                       &y_out,
          Eigen::VectorX<T>                       &x_out,
          Eigen::VectorXi                         &indices)
{
    using RowIter = Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator;

    indices.resize(y_in.size()); indices.setConstant(-1);
    int nv {};

    // reindex variables
    for (int i = 0; i < y_in.size(); ++i)
        if (!is_fixed(i)) indices(i) = nv++;

    // reduced linear system
    A_out.resize(nv, nv); A_out.setZero();
    y_out.resize(nv);     y_out.setZero();
    x_out.resize(nv);     x_out.setZero();

    std::vector<Eigen::Triplet<T>> coef {};

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

template <class T>
inline int reduce_fixed_constraints(
    const Eigen::SparseMatrix<T> &A_in,
    const Eigen::VectorX<T>      &y_in,
    const Eigen::VectorX<T>      &x_in,
    const Eigen::VectorXi        &is_fixed,
          Eigen::SparseMatrix<T> &A_out,
          Eigen::VectorX<T>      &y_out,
          Eigen::VectorX<T>      &x_out,
          Eigen::VectorXi        &indices)
{
    Eigen::SparseMatrix<T, Eigen::RowMajor> A_rm(A_in);
    return reduce_fixed_constraints(A_rm, y_in, x_in, is_fixed, A_out, y_out, x_out, indices);
}

template <class T>
inline int solve_sparse_LU(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::VectorXi        &C,
    const Eigen::VectorX<T>      &y,
          Eigen::VectorX<T>      &x)
{
    Eigen::VectorXi I;
    Eigen::VectorX<T> x1, y1;
    Eigen::SparseMatrix<T> A1;

    reduce_fixed_constraints(A, y, x, C, A1, y1, x1, I);

    int err = solve_sparse_LU(A1, y1, x1);

    // write back to the original array
    for (int i = 0; i < x.size(); ++i)
        if (!C(i)) x(i) = x1(I(i));

    return err;
}

template <class T>
inline int solve_simplical_LDLT(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::VectorXi        &C,
    const Eigen::VectorX<T>      &y,
          Eigen::VectorX<T>      &x)
{
    Eigen::VectorXi I;
    Eigen::VectorX<T> x1, y1;
    Eigen::SparseMatrix<T> A1;

    reduce_fixed_constraints(A, y, x, C, A1, y1, x1, I);

    int err = solve_simplical_LDLT(A1, y1, x1);

    // write back to the original array
    for (int i = 0; i < x.size(); ++i)
        if (!C(i)) x(i) = x1(I(i));

    return err;
}

template
int solve_simplical_LDLT(
    const Eigen::SparseMatrix<std::complex<double>> &A,
    const Eigen::VectorXi                           &C,
    const Eigen::VectorX<std::complex<double>>      &y,
          Eigen::VectorX<std::complex<double>>      &x);
