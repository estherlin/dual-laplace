#include "massmatrix.h"
#include <vector>
#include <iostream>
#include <igl/doublearea.h>

typedef Eigen::Triplet<double> T;

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // The mass matrix is num verticies x num vertices
  M.resize(F.maxCoeff() + 1);
  M.setZero();  // Set everything to zero so we can just add on

  // Get a vector of DOUBLE areas to help later
  Eigen::VectorXd dbla;
  igl::doublearea(l, 0.0, dbla);

  // Each face contributes to an area
  for (int i = 0; i < F.rows(); i++){
    // Add area contributions for each vertex
    M.diagonal()[F(i,0)] += dbla(i) / 6;
    M.diagonal()[F(i,1)] += dbla(i) / 6;
    M.diagonal()[F(i,2)] += dbla(i) / 6;
  }
}

