#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/readOFF.h>
#include <igl/barycenter.h>
#include <igl/massmatrix.h>
#include <igl/readSTL.h>
#include <igl/parula.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/colon.h>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include "dual_laplacian.h"

// Input polygon
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd B;
Eigen::MatrixXd N;

// Tetrahedralized interior
Eigen::MatrixXd TV;
Eigen::MatrixXi TT;
Eigen::MatrixXi TF;

// Solved
Eigen::VectorXd Z, Z_in;
Eigen::MatrixXi bound_faces, bound_inds, bound_across;



// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

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
    viewer.data().set_face_based(true);
  } else if (key == 'b') {
    std::cout<<"boundary faces"<<std::endl;
    viewer.data().clear();
    viewer.data().set_mesh(b, bound_faces);
    viewer.data().set_face_based(true);
  }

  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  if ((argc > 1) && (strcmp(".sl", argv[1]) == 0)) {
    std::cout << "Reading .stl file..." << std::endl;
    FILE * fp = fopen(argv[1], "r");
    igl::readSTL(fp,V,F,N);
    fclose(fp);
  } else {
    std::cout << "Reading .off file..." << std::endl;
    igl::readOFF(argc>1?argv[1]:"../data/3holes.off",V,F);
  }

  // Tetrahedralize the interior  
  igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414Y", TV,TT,TF);

  // Compute barycenters
  igl::barycenter(TV,TT,B);

  // Find boundary faces
  igl::boundary_facets(TT, bound_faces, bound_inds, bound_across);

  // Find boundary vertices
  Eigen::VectorXi b, IA, IC;
  igl::unique(bound_faces,b,IA,IC);

  // List of all vertex indices
  VectorXi all,in;
  igl::colon<int>(0,TV.rows()-1,all);

  // List of interior indices
  igl::setdiff(all,b,in,IA);

  // Construct and slice up Laplacian
  SparseMatrix<double> L,L_in_in,L_in_b, M;
  dual_laplacian(TV, TT, L, M);

  igl::slice(L,in,in,L_in_in);
  igl::slice(L,in,b,L_in_b);

  // Dirichlet boundary conditions from z-coordinate
  VectorXd bc;
  Z = TV.col(2);
  igl::slice(Z,b,bc);

  // Solve PDE
  SimplicialLLT<SparseMatrix<double > > solver(-L_in_in);
  Z_in = solver.solve(L_in_b*bc);
  // slice into solution
  igl::slice_into(Z_in,in,Z);

  // Alternative, short hand
  igl::min_quad_with_fixed_data<double> mqwf;
  // Linear term is 0
  VectorXd B = VectorXd::Zero(TV.rows(),1);
  // Empty constraints
  VectorXd Beq;
  SparseMatrix<double> Aeq;


  // Our L may be indefinite
  igl::min_quad_with_fixed_precompute(L.eval(),b,Aeq,true,mqwf);
  igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z);

  //std::cout<<Z.rows()<<Z.cols()<<std::endl;
  //std::cout<<Z_in.rows()<<Z_in.cols()<<std::endl;
  //std::cout<<TV.rows()<<TV.cols()<<std::endl;
  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;

  // Print command keys
  std::cout<<R"(
    1-9  cross-sectional views
    b    boundary conditions surface
    r  reset mesh geometry and data
  )"<<std::endl;

  viewer.callback_key_down = &key_down;
  key_down(viewer,'5',0);
  viewer.data().show_lines = true;
  viewer.launch();
}