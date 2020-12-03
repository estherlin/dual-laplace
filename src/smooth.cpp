#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include<Eigen/SparseCholesky>
#include <igl/edge_lengths.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <iostream>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  // U = G;

  // Find edge lengths
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  // Find Mass matrix
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);

  // Find cotangent matrix
  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);
  //igl::cotmatrix(V, F, L);

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
  std::cout << "completed 1 itr" << std::endl;
}
