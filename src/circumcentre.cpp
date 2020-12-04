#include "cotmatrix.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>


void circumcentre3d(
  const Eigen::MatrixXd& A, 
  Eigen::Vector3d& c){

  // Find number of data points
  int n = A.rows();
  c.setZero();

  // If number of data points is 1, the centre is the point
  if (n == 1){
    c = A.row(0);
  } else if (n == 2)
  {
    c = 0.5 * (A.row(0) + A.row(1));

  } else if (n == 3)
  { // Math from: https://gamedev.stackexchange.com/a/60631
    Eigen::Vector3d a = A.row(0);
    Eigen::Vector3d b = A.row(1);
    Eigen::Vector3d c = A.row(2);

    Eigen::Vector3d area = (b-a).cross(c-a);
    Eigen::Vector3d num = (c-a).squaredNorm()*area.cross(b-a) - (b-a).squaredNorm()*area.cross(c-a);
    double den = 2*area.squaredNorm();

    c = a + (1.0/den)*num;

  } else if (n == 4)
  { // Math from: https://en.wikipedia.org/wiki/Tetrahedron#Circumcenter
    Eigen::Matrix3d Q(3,3);
    Q.row(0) = A.row(1) - A.row(0);
    Q.row(1) = A.row(2) - A.row(0);
    Q.row(2) = A.row(3) - A.row(0);

    Eigen::Vector3d b;
    b(0) = A.row(1).squaredNorm() - A.row(0).squaredNorm();
    b(1) = A.row(2).squaredNorm() - A.row(0).squaredNorm();
    b(2) = A.row(3).squaredNorm() - A.row(0).squaredNorm();

    c = Q.colPivHouseholderQr().solve(b);

  } else {
      std::cout << "Error: Can't handle circumcentre of more than 4 points" << std::endl;
  }

}