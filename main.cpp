#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/barycenter.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/parula.h>
#include <igl/readOBJ.h>
#include <igl/parula.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/colon.h>
#include <igl/faces_first.h>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include "dual_laplacian.h"
#include "smooth.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;
Eigen::MatrixXd N;

// Tetrahedralized interior
Eigen::MatrixXd TV, tempTV;
Eigen::MatrixXi TT, tempTT;
Eigen::MatrixXi TF, tempTF;
Eigen::MatrixXd C;

// Eigen Problem
Eigen::VectorXd Z, Z_in;
Eigen::MatrixXi bound_faces, bound_inds, bound_across;
Eigen::MatrixXd bound_ver;

// Smoothing Problem
double lambda = 0.1;
Eigen::MatrixXd G, U;;


// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;
  std::cout<<key<<std::endl;

  if (key >= '1' && key <= '9')
  {
    double t = double((key - '1')+1) / 9.0;

    VectorXd v = B.col(2).array() - B.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    vector<int> s;

    for (unsigned i=0; i<v.size();++i)
      if (v(i) < t)
        s.push_back(i);

    MatrixXd V_temp(s.size()*4,3);
    MatrixXi F_temp(s.size()*4,3);
    VectorXd Z_temp(s.size()*4);

    for (unsigned i=0; i<s.size();++i)
    {
      V_temp.row(i*4+0) = TV.row(TT(s[i],0));
      V_temp.row(i*4+1) = TV.row(TT(s[i],1));
      V_temp.row(i*4+2) = TV.row(TT(s[i],2));
      V_temp.row(i*4+3) = TV.row(TT(s[i],3));

      Z_temp(i*4+0) = Z(TT(s[i],0));
      Z_temp(i*4+1) = Z(TT(s[i],1));
      Z_temp(i*4+2) = Z(TT(s[i],2));
      Z_temp(i*4+3) = Z(TT(s[i],3));

      F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
      F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
      F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
      F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
    }

    viewer.data().clear();
    viewer.data().set_mesh(V_temp,F_temp);
    viewer.data().set_data(Z_temp);
    viewer.data().set_colors(C);
    viewer.data().set_face_based(true);
  } else if (key == 'b' || key == 'B') {
    viewer.data().clear();
    viewer.data().set_colors(C);
    viewer.data().set_mesh(V, F);
    viewer.data().set_data(Eigen::VectorXd::Ones(V.rows(),1));
    viewer.data().set_face_based(true);
  } else if (key == 'n' || key == 'N'){
    smooth(TV, TT, U, lambda, U);
    viewer.data().set_mesh(TV,TF);
    viewer.data().compute_normals();
    Eigen::MatrixXd CC;
    igl::parula(U,G.minCoeff(),G.maxCoeff(),CC);
    viewer.data().set_colors(CC);
  }

  return false;
}

void solve_eigen(){
  // Compute barycenters
  igl::barycenter(TV,TT,B);

  // Find boundary faces
  igl::boundary_facets(TT, bound_faces, bound_inds, bound_across);

  // Find boundary vertices
  Eigen::VectorXi b, IA, IC;
  igl::unique(bound_faces,b,IA,IC);

  // List of all vertex indices
  Eigen::VectorXi all,in;
  igl::colon<int>(0,TV.rows()-1,all);

  // List of interior indices
  igl::setdiff(all,b,in,IA);

  // Construct and slice up Laplacian
  Eigen::SparseMatrix<double> L,L_in_in,L_in_b, M;
  // Swap to see primal implmenetation
  //igl::cotmatrix(TV, TT, L);
  //igl::massmatrix(TV, TT, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
  dual_laplacian(TV, TT, L, M);

  igl::slice(L,in,in,L_in_in);
  igl::slice(L,in,b,L_in_b);

  // Dirichlet boundary conditions from z-coordinate
  Eigen::VectorXd bc;
  Z = TV.col(2);
  igl::slice(Z,b,bc);

  // Solve PDE
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(L_in_in);
  Z_in = solver.solve(L_in_b*bc);
  // slice into solution
  igl::slice_into(Z_in,in,Z);

  // Alternative, short hand
  igl::min_quad_with_fixed_data<double> mqwf;
  // Linear term is 0
  Eigen::VectorXd B = Eigen::VectorXd::Zero(TV.rows(),1);
  // Empty constraints
  Eigen::VectorXd Beq;
  Eigen::SparseMatrix<double> Aeq;


  // Our L may be indefinite
  igl::min_quad_with_fixed_precompute(L.eval(),b,Aeq,true,mqwf);
  igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z);
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  igl::readOFF(argc>1?argv[1]:"../data/cube.off",V,F);

  std::cout<<V.rows()<<V.cols()<<std::endl;
  std::cout<<F.rows()<<F.cols()<<std::endl;

  // Tetrahedralize the interior  "pq1.1a0.05aA"
  Eigen::VectorXi IM;
  igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414a0.01", TV,TT,TF);
  //igl::faces_first(TV,F,IM);
  //TT = TT.unaryExpr(bind1st(mem_fun( static_cast<VectorXi::Scalar& (VectorXi::*)(VectorXi::Index)>(&VectorXi::operator())),&IM)).eval();

  // Solve the eigenvalue problem
  solve_eigen();

  // Smoothing
  G = TV.col(1);
  G += 0.1*(G.maxCoeff()-G.minCoeff())*
        Eigen::MatrixXd::Random(G.rows(),G.cols());
  U = G;

  //Set color scale
  igl::parula(Z,true,C);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;

  // Print command keys
  std::cout<<R"(
    1-9  cross-sectional views
    b    see full shape
    f    full shape
    n    smoothing
  )"<<std::endl;

  viewer.callback_key_down = &key_down;
  key_down(viewer,'5',0);
  viewer.data().show_lines = false;
  viewer.launch();
}