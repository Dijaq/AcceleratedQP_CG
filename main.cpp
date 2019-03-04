#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
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
	int num_iter = 20;
	double TolX = 1e-10;
	double TolFun = 1e-6;

	Param_State mesh;

//Seccion of read a mesh
    read_mesh_2D("data_gecko/V.csv", "data_gecko/F.csv","data_gecko/eq_lhs.csv", "data_gecko/eq_rhs.csv", mesh);
    print_dimensions("V: ", mesh.V);

    /*cout << "F: " << mesh.F.rows() << " - " << mesh.F.cols() << endl;
    cout << "V: " << mesh.V.rows() << " - " << mesh.V.cols() << endl;
    cout << "eq_lhs: " << mesh.eq_lhs.rows() << " - " << mesh.eq_lhs.cols() << endl;
    cout << "eq_rhs: " << mesh.eq_rhs.rows() << " - " << mesh.eq_rhs.cols() << endl;*/

    update_F(mesh.F);

    /*MatrixXd mV(mesh.V.rows(), mesh.V.cols()+1);
    create_column_zeros(mesh.V, mV);
    igl::writeOBJ("Isodist_gecko.obj", mV, mesh.F);*/

    MatrixXd V0; 
    //Setup optimization problem
    cout << "Optim Problem" << endl;
    auto t01 = std::chrono::high_resolution_clock::now();
    OptimProblemIsoDist optimProblem(mesh, V0, 25);
    //exportOptimProblem(optimProblem);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto durationOP = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count();
    cout << "End Optim Problem: " << durationOP << endl;

    cout << optimProblem.T.rows() << " - "<< optimProblem.T.cols()<<endl;
    cout << optimProblem.areas.rows()<< " - "<< optimProblem.areas.cols() << endl;

    //Setup solver
    auto t11 = std::chrono::high_resolution_clock::now();
    OptimSolverAcclQuadProx optimProblemAQP("AQP", optimProblem, true, true, true);
    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Optim Problem Acc Quad Prox: " << duration << endl;
    optimProblemAQP.setKappa(1000);

    vector<OptimSolverAcclQuadProx> listOptimSolverAccQuadProx;
    listOptimSolverAccQuadProx.push_back(optimProblemAQP);

    int n_solvers = listOptimSolverAccQuadProx.size();

    for(int i=0; i<n_solvers; i++)
    {
        cout << "Start Iterations: " << endl;
        auto t21 = std::chrono::high_resolution_clock::now();
        listOptimSolverAccQuadProx[i].solveTol(TolX, TolFun, num_iter);
        auto t22 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t22 - t21).count();
        cout << "Optim Solver Iter Acc Quad Prox: " << duration2 << endl;
    }

    print_dimensions("x", optimProblemAQP.x);
    export_mat_to_excel(optimProblemAQP.x, "xFinal"); 

    cout << "Finish Program" << endl;
    /*cout << "Start Practice Section" << endl;

    MatrixXd practice(3,3);
    practice(0,0) = 3;
    practice(0,2) = 1;
    practice(1,0) = -1;
    practice(1,1) = 2;
    practice(2,1) = 1;
    practice(2,2) = -5;

    cout << "Original Matrix" << endl;
    cout << practice << endl;

    FullPivLU<MatrixXd> lu(practice);
    cout << lu.matrixLU() << endl;
    cout << "print matrix L " << endl;
    MatrixXd l = MatrixXd::Identity(3,3);
    l.block<3,3>(0,0).triangularView<StrictlyLower>() = lu.matrixLU();
    cout << l << endl;

    cout << "print matrix U " << endl;
    MatrixXd u = lu.matrixLU().triangularView<Upper>();
    cout << u << endl;
    //cout << lu.permutationQ() << endl;
    cout << "Original Matrix: " << endl;
    cout << lu.permutationP().inverse()*l*u*lu.permutationQ().inverse() << endl;

    cout << lu.permutationP().inverse() << endl;*/

    

    /*for(int i=0; i<practice.rows(); i++)
    {
        for(int j=0; j<practice.cols(); j++)
        {
            practice(i,j) = i+j;
        }
    }*/

    //practice.block<100,100>(0,0) = practice.block<100,100>(3,3);

    //cout << practice.block(0,0,10,10) << endl;

    //cout << practice.row(0)*2 << endl;

    /*SparseMatrix<double> practiceSparce(1000,1000);

    for(int i=0; i<practiceSparce.rows(); i++)
    {
        for(int j=0; j<practiceSparce.cols(); j++)
        {
            practiceSparce.insert(i,j) = i+j;
        }
    }

    MatrixXd spars = practiceSparce;

    spars.block<100,100>(0,0) = spars.block<100,100>(3,3);

    auto t11 = std::chrono::high_resolution_clock::now();
    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Time: " << duration << endl;
    */

    print_dimensions("Fx0: ", optimProblem.x0);
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
    create_column_zeros(exportVAQP, nAQPV);
    igl::writeOBJ("AQP_cpp.obj", nAQPV, mesh.F);
  
	return 0;
}

void exportOptimProblem(OptimProblemIsoDist optimProblem)
{
    /*export_sparsemat_to_excel(optimProblem.T, "T");
    export_sparsemat_to_excel(optimProblem.eq_lhs, "eq_lhs");
    export_sparsemat_to_excel(optimProblem.H, "H");
    export_mat_to_excel(optimProblem.eq_rhs, "eq_rhs");
    export_mat_to_excel(optimProblem.x0, "x0"); 
    cout << "LOG" << endl;
    cout << "n_vars: "<< optimProblem.n_vars << endl;
    cout << "n_eq: "<< optimProblem.n_eq << endl;
    cout << "ENDLOG" << endl;*/

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
  /*int n = 3;
    VectorXd x(n), b(n);
    SparseMatrix<double> A(3,3);
    A.insert(0,0) = 2;
    A.insert(0,1) = -1;
    A.insert(0,2) = -8;
    A.insert(1,0) = 1;
    A.insert(1,1) = 3;
    A.insert(1,2) = -3;
    A.insert(2,0) = 7;
    A.insert(2,1) = 4;
    A.insert(2,2) = -6;

    SparseLU<SparseMatrix<double>> solver;

    solver.analyzePattern(A);
    solver.factorize(A);
    x = solver.solve(b);

    cout << solver(0,0) << endl;
    cout << solver.cols() << endl;*/

    /*VectorXd y = b;
    cout << solver.m_mapU;*/

    /*cout << b << endl;
    cout << x << endl;*/

    /*auto t11 = std::chrono::high_resolution_clock::now();
    SparseMatrix<double> T(mesh.eq_lhs.rows(),mesh.eq_lhs.cols());
    for(int i=0; i<T.rows(); i++)
    {
        for(int j=0; j<T.cols(); j++)
        {
            if(mesh.eq_lhs(i,j) != 0)
                T.insert(i,j) = mesh.eq_lhs(i,j);
        }

    }

    auto t12 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Time: " << duration << endl;*/



   
    /*MatrixXd nV(mesh.V.rows(), mesh.V.cols()+1);
    create_column_zeros(mesh.V, nV);
    mesh.V = nV;
    
    //Meshlab Render->shader->gooch
    igl::writeOBJ("prueba.obj", mesh.V, mesh.F);
    */
   

   /*MatrixXd A(3,3);
    A << 1,4,-3,-2, 8, 5, 3, 4, 7;
    MatrixXd desLU = A;
    FullPivLU<Ref<MatrixXd> > lu(desLU);
    
    cout << A << endl;
    cout << desLU << endl;
    MatrixXd p = lu.permutationP();
    MatrixXd q = lu.permutationQ();
    //cout << A << endl;  */