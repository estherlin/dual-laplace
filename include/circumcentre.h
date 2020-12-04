#ifndef CIRCUMCENTRE3D_H
#define CIRCUMCENTRE3D_H
#include <Eigen/Core>

// Inputs:
//   A  #A by d list of vertex positions - each row is a point
// Outputs:
//   c  #c by d vertex of circumcentre position
void circumcentre3d(
  const Eigen::MatrixXd& A, 
  const Eigen::VectorXd& c);
#endif