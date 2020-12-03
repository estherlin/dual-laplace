#include "dual_laplacian.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <igl/copyleft/tetgen/tetrahedralize.h>

typedef Eigen::Triplet<double> T;

void dual_laplacian(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& T, 
  Eigen::SparseMatrix<double>& L, 
  Eigen::SparseMatrix<double>& M){

  // First call tetrahedralize to mesh the interior of a surface mesh (V,F)
}



