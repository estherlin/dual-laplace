#include "dual_laplacian.h"
#include "tet_volume.h"
#include "circumcentre3d.h"
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <vector>
#include <iostream>


typedef Eigen::Triplet<double> t;

void dual_laplacian(
  const Eigen::MatrixXd& V, 
  const Eigen::MatrixXi& T, 
  Eigen::SparseMatrix<double>& L, 
  Eigen::SparseMatrix<double>& M){

  // For each vertex of tet, you have 3 possible orientations, 2 items for each
  std::vector<t> tripletListL;
  tripletListL.reserve(4*3*4*T.rows()); 

  // Also make the mass matrix in the same loop
  std::vector<t> tripletListM;
  tripletListM.reserve(4*3*2*T.rows()); 

  // List of 12 triangle face orientations per tetrahedron
  // Use this list of faces instead of 2 for loops over indices
  Eigen::MatrixXi faces(12,3);
  faces << 0,1,2,
           0,2,3,
           0,3,1,
           1,0,3,
           1,2,0,
           1,3,2,
           2,0,1,
           2,1,3,
           2,3,0,
           3,0,2,
           3,1,0,
           3,2,1;

  // Loop over each tetrahedron in mesh
  for (int i = 0; i < T.rows(); i++){

    // Get our tetrahedron
    Eigen::MatrixXd tet(4,3);
    for (int ind = 0; ind < 4; ind++){
      tet.row(ind) = V.row(T(i, ind));
    }

    // Get circumcentre of our tetrahedron
    Eigen::Vector3d tet_cc;
    circumcentre3d(tet, tet_cc);

    // Loop over each of our triangular faces
    for (int j = 0; j < 12; j++){

      // Get our triangle face (normal out)
      Eigen::MatrixXd tri(3,3);
      for (int ind = 0; ind < 3; ind++){
        tri.row(ind) = V.row(T(i, faces(j, ind)));
      }

      // Get circumcentre of our triangle
      Eigen::Vector3d tri_cc;
      circumcentre3d(tri, tri_cc);

      // Get circumcentre of our edge
      Eigen::Vector3d edge_cc;
      edge_cc = (0.5) * (V.row(T(i, faces(j, 0))) + V.row(T(i, faces(j, 1))));

      // Get volume of tet
      double vol;
      tet_volume(edge_cc, tri_cc, tet_cc, V.row(T(i, faces(j, 0))), vol);

      // Place volume into our M operator
      tripletListM.push_back(t(T(i, faces(j, 0)), T(i, faces(j, 0)), vol)); // i != j
      tripletListM.push_back(t(T(i, faces(j, 1)), T(i, faces(j, 1)), vol)); // i != j

      // Get edge norm and weight of this tetrahedron
      double edge_sq = (V.row(T(i, faces(j, 0))) - V.row(T(i, faces(j, 1)))).squaredNorm();
      double w = (6.0) * vol/edge_sq;

      // Place them into our L operator
      tripletListL.push_back(t(T(i, faces(j, 0)), T(i, faces(j, 1)), w)); // i != j
      tripletListL.push_back(t(T(i, faces(j, 1)), T(i, faces(j, 0)), w)); // i != j
      tripletListL.push_back(t(T(i, faces(j, 0)), T(i, faces(j, 0)), (-1.0)*w)); // i == i
      tripletListL.push_back(t(T(i, faces(j, 1)), T(i, faces(j, 1)), (-1.0)*w)); // j == j

    }
  }

  // The matrix is num verticies x num vertices
  L.resize(V.rows(), V.rows());
  M.resize(V.rows(), V.rows());

  L.setFromTriplets(tripletListL.begin(), tripletListL.end());
  M.setFromTriplets(tripletListM.begin(), tripletListM.end());
}



