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
//
// Notes: 
// Math for the triangle circumcentre in 3D was taken from 
// https://math.stackexchange.com/a/2130522
//
// An excellent discussion on all things circumcentre 
// and C implementation examples from Jonathan Shewchuk here:
// https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
//
// Math for the tetrahedron circumcentre was taken from 
// https://en.wikipedia.org/wiki/Tetrahedron#Circumcenter
void circumcentre3d(
  const Eigen::MatrixXd& A, 
  Eigen::Vector3d& c);
#endif