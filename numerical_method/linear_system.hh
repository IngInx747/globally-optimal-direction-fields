#ifndef LINEAR_SYSTEM_SOLVER_HH
#define LINEAR_SYSTEM_SOLVER_HH

#include <Eigen/SparseCore>

template <class T>
int dump_sparse_matrix(const Eigen::SparseMatrix<T>&, const char*);

template <class T, size_t N_row, size_t N_col>
int solve_sparse_LU(const Eigen::SparseMatrix<T>&, const Eigen::Matrix<T, N_row, N_col>&, Eigen::Matrix<T, N_row, N_col>&);

template <class T, size_t N_row, size_t N_col>
int solve_simplical_LDLT(const Eigen::SparseMatrix<T>&, const Eigen::Matrix<T, N_row, N_col>&, Eigen::Matrix<T, N_row, N_col>&);

template <class T, size_t N_row, size_t N_col>
int solve_conjugate_gradient(const Eigen::SparseMatrix<T>&, const Eigen::Matrix<T, N_row, N_col>&, Eigen::Matrix<T, N_row, N_col>&);

template <class T>
int solve_inversed_power(const Eigen::SparseMatrix<T>&, const Eigen::SparseMatrix<T>&, Eigen::VectorX<T>&, const double tol, const int n_iter);

#endif