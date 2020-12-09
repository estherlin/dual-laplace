#ifndef CIRCUMCENTRE3D_H
#define CIRCUMCENTRE3D_H
#include <Eigen/Core>

// Determines the circumcentre in 3D of triangles and tetrahedrons
// given a matrix of vertex points (A).
//
// Inputs:
//   A  #A by 3 list of vertex positions in 3D - each row is a vertex
// Outputs:
//   c  #c by 3 vertex of circumcentre position in 3D
void circumcentre3d(
  const Eigen::MatrixXd& A, 
  Eigen::Vector3d& c);
#endif