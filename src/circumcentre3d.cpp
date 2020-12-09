#include "circumcentre3d.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

void circumcentre3d(
  const Eigen::MatrixXd& A, 
  Eigen::Vector3d& c){

  // Find number of vertex points
  int n = A.rows();

  if (n == 1)
  {
    // If number of vertex points is 1, the centre is the point
    c = A.row(0);

  } else if (n == 2)
  {
    // If number of data points is 2, the centre is the midpoint
    c = 0.5 * (A.row(0) + A.row(1));

  } else if (n == 3)
  {
    // If number of vertex points is 3 we have a triangle:
    // Solutions like https://gamedev.stackexchange.com/a/60631
    // aren't good because they use cross products and cross products 
    // caused overflow problems for me. So I use Barycentric coordinates:
    // Math from: https://math.stackexchange.com/a/2130522 and some 
    // adaptations from: https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    double a_sq = (A.row(0)-A.row(1)).squaredNorm();
    double b_sq = (A.row(1)-A.row(2)).squaredNorm();
    double c_sq = (A.row(2)-A.row(0)).squaredNorm();

    Eigen::Vector3d O(3,1);
    O(0) = a_sq * (b_sq + c_sq - a_sq);
    O(1) = b_sq * (a_sq + c_sq - b_sq);
    O(2) = c_sq * (b_sq + a_sq - c_sq);
    double sum = O(0) + O(1) + O(2);

    c = (O(0)/sum)*A.row(2) + (O(1)/sum)*A.row(0) + (O(2)/sum)*A.row(1);

  } else if (n == 4)
  { 
    // If number of vertex points is 4 we have a tetrahedron:
    // Math from: https://en.wikipedia.org/wiki/Tetrahedron#Circumcenter
    Eigen::Matrix3d Q(3,3);
    Q.row(0) = A.row(1) - A.row(0);
    Q.row(1) = A.row(2) - A.row(0);
    Q.row(2) = A.row(3) - A.row(0);

    Eigen::Vector3d b;
    b(0) = A.row(1).squaredNorm() - A.row(0).squaredNorm();
    b(1) = A.row(2).squaredNorm() - A.row(0).squaredNorm();
    b(2) = A.row(3).squaredNorm() - A.row(0).squaredNorm();

    c = (0.5) * Q.colPivHouseholderQr().solve(b);

  } else {
    // If we have more than 4 vertex points, not applicable
    std::cout << "Error: Can't handle circumcentre of more than 4 points" << std::endl;
  }

}