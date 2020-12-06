#ifndef CIRCUMCENTRE3D_H
#define CIRCUMCENTRE3D_H
#include <Eigen/Core>

// Inputs:
//   A  #A by d list of vertex positions in 3D - each row is a point
// Outputs:
//   c  #c by d vertex of circumcentre position in 3D
void circumcentre3d(
  const Eigen::MatrixXd& A, 
  Eigen::Vector3d& c);
#endif