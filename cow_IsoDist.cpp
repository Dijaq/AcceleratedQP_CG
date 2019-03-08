#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include "code/mathelpers.h"
#include "class/Param_State.h"
#include "code/utils.h"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "libs/writeOBJ.h"
#include "code/OptimProblemIsoDist.h"
#include "code/OptimSolverAcclQuadProx.h"
#include <chrono>

using namespace std;
using namespace cv;
using namespace Eigen;

void exportOptimProblem(OptimProblemIsoDist optimProblem);

int main()
{
    	//Example of gecko deformation
	int num_iter =40;
	double TolX = 1e-10;
	double TolFun = 1e-6;

	Param_State mesh;
//Seccion of read a mesh
    read_mesh_3D("data_cow/V.csv", "data_cow/F.csv","data_cow/Vt.csv", "data_cow/Ft.csv", mesh);
    update_F(mesh.F);

    
//    SparseMatrix<double>  eq_lhs = kron_sparse(MatrixXd::Identity(2,2).sparseView(), sparse_to_parameterization(mesh.inds_bF, mesh.V.rows()));
    mesh.eq_lhs = kron_sparse(MatrixXd::Identity(2,2).sparseView(), sparse_to_parameterization(mesh.inds_bF, mesh.V.rows()));
    mesh.eq_rhs = MatrixXd::Zero(2,1);

    cout << mesh.eq_lhs.rows() << endl;
    cout << mesh.eq_lhs.cols() << endl;
    cout << MatrixXd::Identity(2,2) << endl;

    cout << "Init Optim Problem" << endl;
    OptimProblemIsoDist optimProblem(mesh, mesh.Vt, 1); 

    cout << "Init OptimAccelQuad prox" << endl; 
    OptimSolverAcclQuadProx optimProblemAQP("AQP", optimProblem, true, true, true);
    optimProblemAQP.setKappa(1000);

    vector<OptimSolverAcclQuadProx> listOptimSolverAccQuadProx;
    listOptimSolverAccQuadProx.push_back(optimProblemAQP);


    int n_solvers = listOptimSolverAccQuadProx.size();

    cout << "Start Iterations: " << endl;
    for(int i=0; i<n_solvers; i++)
    {
        listOptimSolverAccQuadProx[i].solveTol(TolX, TolFun, num_iter);  
    }

    cout << "Finish Program" << endl;

    /*print_dimensions("V: ", mesh.V);
    //print_dimensions("F: ", mesh.F);
    print_dimensions("Vt: ", mesh.Vt);
    print_dimensions("Ft: ", mesh.Ft);*/
    /*cout << mesh.F.rows() << " - " << mesh.F.cols() << endl;
    cout << mesh.V.rows() << " - " << mesh.V.cols() << endl;*/

/*    MatrixXd V0; 
    //Setup optimization problem
    cout << "Optim Problem" << endl;
   
    OptimProblemIsoDist optimProblem(mesh, V0, 25);
    //exportOptimProblem(optimProblem);   
  
    OptimSolverAcclQuadProx optimProblemAQP("AQP", optimProblem, true, true, true);
    
    optimProblemAQP.setKappa(1000);

    vector<OptimSolverAcclQuadProx> listOptimSolverAccQuadProx;
    listOptimSolverAccQuadProx.push_back(optimProblemAQP);

    int n_solvers = listOptimSolverAccQuadProx.size();

    for(int i=0; i<n_solvers; i++)
    {
        cout << "Start Iterations: " << endl;
        
        listOptimSolverAccQuadProx[i].solveTol(TolX, TolFun, num_iter);
       
    }

    cout << "Finish Program" << endl;
     
    MatrixXd exportV = optimProblem.x0;
    //matrix_reshape(exportV, exportV.rows()/2, 2);
    MatrixXd nV(exportV.rows(), exportV.cols()+1);
    create_column_zeros(exportV, nV);
    igl::writeOBJ("IsoDist_cpp.obj", nV, mesh.F);

    print_dimensions("Fx: ", optimProblemAQP.x);
    //MatrixXd exportVAQP = optimProblemAQP.x;
    MatrixXd exportVAQP = listOptimSolverAccQuadProx[0].x;
    matrix_reshape(exportVAQP, exportVAQP.rows()/2, 2);
    MatrixXd nAQPV(exportVAQP.rows(), exportVAQP.cols()+1);
    create_column_zeros(exportVAQP, nAQPV);*/

    //igl::writeOBJ("cow.obj", mesh.V, mesh.F);

    //Create init parameterization
   
    //igl::writeOBJ("cow_param.obj", mesh.V, mesh.F, Eigen::MatrixXd(), Eigen::MatrixXi(), mesh.Vt, mesh.F);

    MatrixXi Fi(mesh.F.rows(), mesh.F.cols());
    for(int i=0; i<mesh.F.rows(); i++)
    {
        for(int j=0; j<mesh.F.cols(); j++)
        {
            Fi(i,j) = mesh.F(i,j);
        }
    }
    //mesh.Vt = mesh.Vt*10;
    igl::writeOBJ("presentacion_cow/cow_init_0.obj", mesh.V, Fi, Eigen::MatrixXd(), Eigen::MatrixXi(), mesh.Vt, Fi);

    //Create final parameterization
    print_dimensions("xx: ", listOptimSolverAccQuadProx[0].x);
    MatrixXd VerticeTexture = listOptimSolverAccQuadProx[0].x;
    matrix_reshape(VerticeTexture, VerticeTexture.rows()/2, 2);

    igl::writeOBJ("cow_param_x_final.obj", mesh.V, Fi, Eigen::MatrixXd(), Eigen::MatrixXi(), VerticeTexture, Fi);

    /*MatrixXd VerticeText = listOptimSolverAccQuadProx[0].p;
    matrix_reshape(VerticeText, VerticeText.rows()/2, 2);*/
    MatrixXd VerticeText = VerticeTexture/100;

    igl::writeOBJ("cow_param_final.obj", mesh.V, Fi, Eigen::MatrixXd(), Eigen::MatrixXi(), VerticeText, Fi);

	return 0;
}

void exportOptimProblem(OptimProblemIsoDist optimProblem)
{
    export_mat_to_excel(optimProblem.V, "V"); 
    export_mat_to_excel(optimProblem.F, "F"); 
    export_mat_to_excel(optimProblem.areas, "areas"); 
    export_mat_to_excel(optimProblem.Tx, "Tx"); 
    export_mat_to_excel(optimProblem.R, "R"); 
    export_mat_to_excel(optimProblem.Tx_grad, "Tx_grad"); 

    cout << "LOG" << endl;
    cout << "dim: "<< optimProblem.dim << endl;
    cout << "n_vert: "<< optimProblem.n_vert << endl;
    cout << "n_tri: "<< optimProblem.n_tri << endl;
    cout << "f_val: "<< optimProblem.f_val << endl;
    cout << "flips: "<< optimProblem.flips << endl;
    cout << "localHess: "<< optimProblem.localHess << endl;
    cout << "initArapIter: "<< optimProblem.initArapIter << endl;
    cout << "ENDLOG" << endl;

}
