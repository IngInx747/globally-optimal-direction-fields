#ifndef SPARSE_SOLVER_HH
#define SPARSE_SOLVER_HH

#include <Eigen/SparseCore>

template <typename ScalarT>
inline int save_matrix(
    const Eigen::MatrixX<ScalarT> &mat,
    const char *filename,
    const std::streamsize prec = 0);

template <typename ScalarT, int OptionsN>
inline int save_matrix(
    const Eigen::SparseMatrix<ScalarT, OptionsN>&,
    const char *filename,
    const std::streamsize prec = 0);

template <typename ScalarT, int N_row, int N_col>
inline int solve_simplical_LDLT(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &b,
          Eigen::Matrix<ScalarT, N_row, N_col> &x);

template <class ScalarT>
inline int solve_simplical_LDLT(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorX<ScalarT>      &b,
    const Eigen::VectorXi              &c,
          Eigen::VectorX<ScalarT>      &x);

template <class T>
inline int solve_inversed_power(
    const Eigen::SparseMatrix<T> &A,
          Eigen::VectorX<T> &x,
    const double tol,
    const int n_iter);

template <class T>
inline int solve_inversed_power(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::SparseMatrix<T> &B,
          Eigen::VectorX<T> &x,
    const double tol,
    const int n_iter);

#endif