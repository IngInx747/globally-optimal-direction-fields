#ifndef LINEAR_SYSTEM_SOLVER_HH
#define LINEAR_SYSTEM_SOLVER_HH

#include <Eigen/SparseCore>

template <class ScalarT>
int dump_sparse_matrix(const Eigen::SparseMatrix<ScalarT>&, const char*);

template <typename ScalarT, int N_row, int N_col>
int solve_simplical_LDLT(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::Matrix<ScalarT, N_row, N_col> &y,
          Eigen::Matrix<ScalarT, N_row, N_col> &x);

template <class ScalarT>
inline int solve_simplical_LDLT(
    const Eigen::SparseMatrix<ScalarT> &A,
    const Eigen::VectorXi              &C,
    const Eigen::VectorX<ScalarT>      &y,
          Eigen::VectorX<ScalarT>      &x);

template <class T>
int solve_inversed_power(
    const Eigen::SparseMatrix<T> &A,
          Eigen::VectorX<T> &x,
    const double tol,
    const int n_iter);

template <class T>
int solve_inversed_power(
    const Eigen::SparseMatrix<T> &A,
    const Eigen::SparseMatrix<T> &B,
          Eigen::VectorX<T> &x,
    const double tol,
    const int n_iter);

#endif