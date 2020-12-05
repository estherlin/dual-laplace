#ifndef TETVOLUME_H
#define TETVOLUME_H
#include <Eigen/Core>

// Inputs:
//   a  #1 by 3 vertex position of point in 3D
//   b  #1 by 3 vertex position of point in 3D
//   c  #1 by 3 vertex position of point in 3D
//   d  #1 by 3 vertex position of point in 3D
// Outputs:
//   v  #volume of tetrahedron
void tet_volume(
  const Eigen::Vector3d& a, 
  const Eigen::Vector3d& b,
  const Eigen::Vector3d& c,
  const Eigen::Vector3d& d,
  double& vol);
#endif