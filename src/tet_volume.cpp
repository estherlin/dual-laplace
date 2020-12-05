#include "tet_volume.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

double tet_volume(
  const Eigen::Vector3d& a, 
  const Eigen::Vector3d& b,
  const Eigen::Vector3d& c,
  const Eigen::Vector3d& d)
{  
  return (1.0/6.0) * ((a-d).cross(b-d)).dot(c-d);
}