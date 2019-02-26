#ifndef OPTIMSOLVERACCLQUADPROX_H
#define OPTIMSOLVERACCLQUADPROX_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "math.h"
#include <chrono>
#include "../class/DSparseLU.h"
#include "../class/Param_State.h"
#include "../mex/computeMeshTranformationCoeffsMex.h"
#include "OptimSolverIterative.h"

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimSolverAcclQuadProx : public OptimSolverIterative
{
public:

    //Parameters
    MatrixXd x_prev;
    MatrixXd y;
    int p;
    //SparseMatrix<double> KKT;
    DSparseLU KKT;
    MatrixXd KKT_rhs;
    int p_lambda;
    int y_f;
    int y_fgrad;
    int y_start;
    int t_init; //Store initialization time
    int f_count = 0; //Store function evaluation count

    //Solver parameters
    bool useAcceleration;
    bool useQuadProxy;
    bool useLineSearch;
    float theta; // = [];

    //Step size limiting
    bool useAccelerationStepSizeLimit = true;
    float accelarationStepSizeLimitFactor = 0.5;
    float accelerationStepSize;
    bool useLineSearchStepSizeLimit = true;
    float lineSearchStepSizeLimitFactor = 0.5;
    bool useLineSearchStepSizeMemory = true;
    float lineSearchStepSizeMemoryFactor = pow(2,5.5);

    //Line serch parameters
    float ls_alpha = 0.2;
    float ls_beta = 0.5;


    OptimSolverAcclQuadProx();
    OptimSolverAcclQuadProx(string tag, OptimProblemIsoDist optimProblem, bool useAccelaration, bool useQuadProxy, bool useLineSearch);
    ~OptimSolverAcclQuadProx();
    void setKappa(float kappa);
    void solveTol(float TolX, float TolFun, int max_iter);
    void iterate();
    double min(double value1, double value2);
};

OptimSolverAcclQuadProx::OptimSolverAcclQuadProx()
{}

OptimSolverAcclQuadProx::OptimSolverAcclQuadProx(string tag, OptimProblemIsoDist optimProblem, bool useAcceleration, bool useQuadProxy, bool useLineSearch)
{
    //t_init_start = tic

    this->useAcceleration = useAcceleration;
    this->useQuadProxy = useQuadProxy;
    this->useLineSearch = useLineSearch;

    //Init Solver
    initSolver(tag, optimProblem);

    SparseMatrix<double> KKT_mat;

    if(this->useQuadProxy)
    {
        SparseMatrix<double> spar(optimProblem.n_eq, optimProblem.n_eq);
        //auto t11 = std::chrono::high_resolution_clock::now();
        KKT_mat = join_matrices(optimProblem.H, optimProblem.eq_lhs.transpose(), optimProblem.eq_lhs, spar);
        //auto t12 = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
        //cout << "Time: " << duration << endl;
        /*cout << optimProblem.H.rows() << " - " << optimProblem.H.cols()<<endl;
        cout << optimProblem.eq_lhs.rows() << " - " << optimProblem.eq_lhs.cols() <<endl;
        cout << optimProblem.n_eq << endl;*/
    }
    else
    {
        SparseMatrix<double> spar(optimProblem.n_eq, optimProblem.n_eq);
        SparseMatrix<double> ones = create_SparseMatrix_ones(optimProblem.H.rows(), optimProblem.H.cols());
        KKT_mat = join_matrices(ones, optimProblem.eq_lhs.transpose(), optimProblem.eq_lhs, spar);
    }

    /*SparseMatrix<double> smatrix(3,3);
    smatrix.insert(2,1) = -5;
    smatrix.insert(0,0) = 3;
    smatrix.insert(0,2) = 2;
    smatrix.insert(1,0) = 1;
    smatrix.insert(1,1) = -2;

    MatrixXd matrix = smatrix;
    DSparseLU sLU(matrix);
    cout << smatrix << endl;
    cout << sLU.U << endl;
    cout << sLU.L << endl;
    cout << sLU.P << endl;
    cout << sLU.Q << endl;*/    

    MatrixXd matrix = KKT_mat;
    DSparseLU sparseLU(matrix);
    this->KKT = sparseLU;
    //this->KKT = SparseLU(matrix);

    //init internal variables
    MatrixXd initZeros(this->optimProblem.n_vars+this->optimProblem.n_eq,1);
    this->KKT_rhs = initZeros;
    this->x_prev = this->x;
    this->t = 1/this->lineSearchStepSizeMemoryFactor;
    //Store init time
    //t_init = toc(t_init_start)
}

OptimSolverAcclQuadProx::~OptimSolverAcclQuadProx(){}

void OptimSolverAcclQuadProx::setKappa(float kappa)
{
    this->theta = (1-sqrt(1/kappa))/(1+sqrt(1/kappa));
}

void OptimSolverAcclQuadProx::solveTol(float TolX, float TolFun, int max_iter)
{
    //This is only logs is not necesarialy for the implementation
    //logInit();
    //logState();

    //run num_iter iteration
    //for(int i=0; i<max_iter; i++)
    for(int i=0; i<1; i++)
    {
        //this is the time to start no necesarialy
        //t_iter_start = tic

        iterate();
        //Stop criteria
    }
}

void OptimSolverAcclQuadProx::iterate()
{
    //Use Acceleration
    if(this->useAcceleration)
    {
        if(this->useAccelerationStepSizeLimit)
        {
            this->accelerationStepSize = min(this->theta, this->accelarationStepSizeLimitFactor*this->optimProblem.getMaxStep(this->x, (this->x-this->x_prev)));
        }
        else
        {
            this->accelerationStepSize = this->theta;
            this->y = this->x+this->accelerationStepSize*(this->x-this->x_prev);
        }
    }
    else
    {
        this->y = this->x;
    }

    //Quadratic proxy minimization
    if(this->useLineSearch)
    {
        cout << "useLineSearch" << endl;
    }
    else
    {

    }

    //Initialize step size
    if(this->useLineSearchStepSizeMemory)
    {
        cout << "useLineSearchStepSizeMemory" << endl;
    }
    else
    {

    }

    if(this->useLineSearchStepSizeLimit)
    {
        cout << "useLineSearchStepLimit" << endl;
    }
    else
    {

    }

    //Line search
    if(this->useLineSearch)
    {
        cout << "useLineSearch" << endl;
    }

    //Update values

}

double OptimSolverAcclQuadProx::min(double value1, double value2)
{
    if(value1 < value2)
        return value1;
    else
        return value2;
}

#endif