#include "cotmatrix.h"
#include <vector>
#include <iostream>
#include <cmath>

typedef Eigen::Triplet<double> T;

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // The mass matrix is num verticies x num vertices
  L.resize(F.maxCoeff()+1, F.maxCoeff()+1);
  L.setZero(); // Always set to zero!

  std::vector<T> tripletList;
  tripletList.reserve(9*F.rows()); // 3x3 matrix for each triangle

  // For each triangular face we constuct a 3x3 matrix
  for (int i = 0; i < F.rows(); i++){
    // Get edge lengths
    double a = l(i,0); // edge opposite angle A
    double b = l(i,1); // edge opposite angle B
    double c = l(i,2); // edge opposite angle C

    // calculate area of trangle face
    double s = (a+b+c)/2;
    double area = sqrt(s*(s-a)*(s-b)*(s-c));

    // Calculate the cotangents already divided by 2
    double cotA = (b*b+c*c-a*a)/(2*4*area); // vertices: 1 <-> 2
    double cotB = (a*a+c*c-b*b)/(2*4*area); // vertices: 2 <-> 0
    double cotC = (a*a+b*b-c*c)/(2*4*area); // vertices: 0 <-> 1

    // Non diagonals, we only need to add one alpha for i->j because the 
    // beta in other j-> i will be taken care of in another loop
    tripletList.push_back(T(F(i, 0), F(i, 1), cotC));
    tripletList.push_back(T(F(i, 1), F(i, 0), cotC));
    tripletList.push_back(T(F(i, 1), F(i, 2), cotA));
    tripletList.push_back(T(F(i, 2), F(i, 1), cotA));
    tripletList.push_back(T(F(i, 2), F(i, 0), cotB));
    tripletList.push_back(T(F(i, 0), F(i, 2), cotB));
    
    // Diagonals: Add the other components
    tripletList.push_back(T(F(i, 0), F(i, 0), (-1.0)*(cotB+cotC)));
    tripletList.push_back(T(F(i, 1), F(i, 1), (-1.0)*(cotA+cotC)));
    tripletList.push_back(T(F(i, 2), F(i, 2), (-1.0)*(cotA+cotB)));
  }

  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

