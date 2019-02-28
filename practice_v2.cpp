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

void multiplicar_matrices(MatrixXd &result, MatrixXd matrixA, MatrixXd matrixB);

int main()
{
	//Example of gecko deformation
	int num_iter = 2500;
	double TolX = 1e-10;
	double TolFun = 1e-6;

	Param_State mesh;

//Seccion of read a mesh
    read_mesh_2D("data_gecko/V.csv", "data_gecko/F.csv","data_gecko/eq_lhs.csv", "data_gecko/eq_rhs.csv", mesh);
    print_dimensions("V: ", mesh.V);
    update_F(mesh.F);

    int cols = 200;
    int rows = 200;

    MatrixXd matrix(rows, cols);
    for(int i=0; i<matrix.rows(); i++)
    {
        for(int j=0; j<matrix.cols(); j++)
        {
            matrix(i,j) = i+j;
        }
    }

    MatrixXd mul = MatrixXd::Zero(rows,cols);
    SparseMatrix<double> mulS(rows, cols);

    auto t11 = std::chrono::high_resolution_clock::now();
    //mulS = (matrix.transpose()*matrix).sparseView();
    mul = matrix*matrix;
    //multiplicar_matrices(mul, matrix, matrix*matrix);
    auto t12 = std::chrono::high_resolution_clock::now();


    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Time Mul Matrix: " << duration << endl;

  
	return 0;
}

/*
Al multiplicar
*/
void multiplicar_matrices(MatrixXd &result, MatrixXd matrixA, MatrixXd matrixB)
{
    /*for(int i=0; i<result.rows(); i++)  m
    {
        for(int j=0; j<result.cols(); j++)
        {
            result(i,j) = matrixA(i,j)*matrixB(i,j);
        }
    }
}


 