#include "smooth.h"
#include<Eigen/SparseCholesky>
#include <igl/edge_lengths.h>
#include "dual_laplacian.h"
#include <iostream>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Find dual laplacian matrix
  Eigen::SparseMatrix<double> L, M;
  dual_laplacian(V, F, L, M);

  // Find Ax = b;
  Eigen::SparseMatrix<double> A = (-1.0)*lambda * L;
  for (int i = 0; i < A.rows(); i++){
    // Need to add the items of M onto the diagonal of A
    // Why can't we just use a sparse matrix for M?
  	A.coeffRef(i, i) += M.diagonal()[i];
  }
  Eigen::MatrixXd b = M*G;  // Why is b not sparse?

  // People suggested SimplicialLLT -- Bicc doesn't work?
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  U = solver.solve(b);
}
