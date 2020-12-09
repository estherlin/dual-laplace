#include "tet_volume.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

void tet_volume(
  const Eigen::Vector3d& a, 
  const Eigen::Vector3d& b,
  const Eigen::Vector3d& c,
  const Eigen::Vector3d& d,
  double& vol)
{  
  // Math from: https://en.wikipedia.org/wiki/Tetrahedron#Volume
  vol = (1.0/6.0) * (a-d).dot((b-d).cross(c-d));
}