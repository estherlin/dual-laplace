#include "dual_laplacian.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <igl/copyleft/tetgen/tetrahedralize.h>

typedef Eigen::Triplet<double> t;

void dual_laplacian(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& T, 
  Eigen::SparseMatrix<double>& L, 
  Eigen::SparseMatrix<double>& M){

  // The matrix is num verticies x num vertices
  L.resize(T.maxCoeff()+1, T.maxCoeff()+1);
  L.setZero(); // Always set to zero!

  std::vector<t> tripletList;
  tripletList.reserve(16*T.rows()); // 4x4 matrix for each tetrahedron


  double v_ijkl, v_ijlk;
}



