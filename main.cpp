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

int main()
{
	//Example of gecko deformation
	int num_iter = 2500;
	double TolX = 1e-10;
	double TolFun = 1e-6;

	Param_State mesh;

//Seccion of read a mesh
    read_mesh_2D("data_gecko/Vx0.csv", "data_gecko/F.csv","data_gecko/eq_lhs.csv", "data_gecko/eq_rhs.csv", mesh);
    print_dimensions("V: ", mesh.V);

    /*cout << "F: " << mesh.F.rows() << " - " << mesh.F.cols() << endl;
    cout << "V: " << mesh.V.rows() << " - " << mesh.V.cols() << endl;
    cout << "eq_lhs: " << mesh.eq_lhs.rows() << " - " << mesh.eq_lhs.cols() << endl;
    cout << "eq_rhs: " << mesh.eq_rhs.rows() << " - " << mesh.eq_rhs.cols() << endl;*/

    update_F(mesh.F);

    MatrixXd mV(mesh.V.rows(), mesh.V.cols()+1);
    create_column_zeros(mesh.V, mV);
    igl::writeOBJ("Isodist_gecko.obj", mV, mesh.F);

    MatrixXd V0; 
    //Setup optimization problem
    OptimProblemIsoDist optimProblem(mesh, V0, 25);

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
        auto t21 = std::chrono::high_resolution_clock::now();
        listOptimSolverAccQuadProx[i].solveTol(TolX, TolFun, num_iter);
        auto t22 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(t22 - t21).count();
        cout << "Optim Solver Iter Acc Quad Prox: " << duration2 << endl;
    }

    print_dimensions("x0", optimProblemAQP.x);

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

    MatrixXd exportV = optimProblemAQP.x;
    matrix_reshape(exportV, exportV.rows()/2, 2);

    MatrixXd nV(exportV.rows(), exportV.cols()+1);
    create_column_zeros(exportV, nV);
    igl::writeOBJ("prueba_it2.obj", nV, mesh.F);
  
	return 0;
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