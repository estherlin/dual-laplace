#ifndef TETVOLUME_H
#define TETVOLUME_H
#include <Eigen/Core>

// Determines the signed volume of a tetrahedron given the four 
// vertex points of the tetrahedron (a,b,c,d) in 3D.
//
// Inputs:
//   a  #1 by 3 vertex position of point in 3D
//   b  #1 by 3 vertex position of point in 3D
//   c  #1 by 3 vertex position of point in 3D
//   d  #1 by 3 vertex position of point in 3D
// Outputs:
//   v  #double signed volume of tetrahedron
//
// Note: This implementation assume that the user chooses to 
// centre the origin of the coordinate system with vertex d
// as described here: https://en.wikipedia.org/wiki/Tetrahedron#Volume
void tet_volume(
  const Eigen::Vector3d& a, 
  const Eigen::Vector3d& b,
  const Eigen::Vector3d& c,
  const Eigen::Vector3d& d,
  double& vol);
#endif