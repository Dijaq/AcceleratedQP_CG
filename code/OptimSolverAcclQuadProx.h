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
#include "../class/Splu.h"
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
    VectorXd p;
    //SparseMatrix<double> KKT;
    //DSparseLU KKT;
    Splu KKT;
    MatrixXd KKT_rhs;
    VectorXd p_lambda;
    double y_f;
    VectorXd y_fgrad;
    float t_start;
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
    void computeLineSearchCond(double &linesearch_cond_lhs, double &linesearch_cond_rhs);
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
        //cout << "Time Sparse: " << duration << endl;
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

    /*auto t31 = std::chrono::high_resolution_clock::now();
    DSparseLU ssparseLU(matrix, true);
    auto t32 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t32 - t31).count();
    cout << "Time AQP1: " << duration << endl;*/

    MatrixXd matrix = KKT_mat;
    auto t11 = std::chrono::high_resolution_clock::now();
    //SparseLUPQ sparseLUPQ(KKT_mat);
    //this->KKT = SparseLUPQ(KKT_mat);
    //this->KKT = DSparseLU(matrix);
    this->KKT = Splu(KKT_mat);

    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Time AQP2: " << duration1 << endl;



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
    for(int i=0; i<100; i++)
    {
        //this is the time to start no necesarialy
        //t_iter_start = tic

        iterate();
        //Falta Stop criteria

    }
}

void OptimSolverAcclQuadProx::iterate()
{
    //Use Acceleration
    if(this->useAcceleration)
    {
        //cout << "useAcceleration" << endl;
        if(this->useAccelerationStepSizeLimit)
        {
            MatrixXd tempX = this->x;
            MatrixXd tempX_prev = this->x_prev;
            matrix_reshape(tempX, tempX.rows()/this->optimProblem.dim, this->optimProblem.dim);
            matrix_reshape(tempX_prev,tempX_prev.rows()/this->optimProblem.dim, this->optimProblem.dim);
            
            this->accelerationStepSize = min(this->theta, this->accelarationStepSizeLimitFactor*this->optimProblem.getMaxStep(tempX, (tempX-tempX_prev)));
        }
        else
        {
            this->accelerationStepSize = this->theta;
        }
        this->y = this->x+this->accelerationStepSize*(this->x-this->x_prev);
    }
    else
    {
        this->y = this->x;
    }

    //Quadratic proxy minimization
    if(this->useLineSearch)
    {
        //cout << "useLineSearch" << endl;
        this->optimProblem.evaluateValueGrad(this->y, this->y_f, this->y_fgrad);
        this->f_count++;
    }
    else
    {
        this->optimProblem.evaluateGrad(this->y, this->y_f, this->y_fgrad);
        this->f_count++;
    }

    //print_dimensions("->", this->KKT_rhs);//Dimension of 1728
    //print_dimensions("->", this->y_fgrad);
    //cout << this->optimProblem.n_vars << endl;
    for(int i=0; i<this->optimProblem.n_vars; i++)
    {
        this->KKT_rhs(i,0) = -this->y_fgrad(i,0);
    }

    /*Prueba matrix por descompistion LU and solve
    MatrixXd smatrix(3,3);
    smatrix(2,1) = -5;
    smatrix(0,0) = 3;
    smatrix(0,2) = 2;
    smatrix(1,0) = 1;
    smatrix(1,1) = -2;
    cout << smatrix << endl;
    DSparseLU sLU(smatrix);
    cout << sLU.U << endl;
    cout << sLU.L << endl;
    cout << sLU.P.transpose()*sLU.L*sLU.U << endl;
    VectorXd vec(3,1);
    vec(0,0) = -2;
    vec(1,0) = 3;
    vec(2,0) = 6;
    sLU.solve(vec);*/

    /*print_dimensions("->", this->KKT_rhs);
    print_dimensions("->", this->KKT.L);
    */

//    print_dimensions("->", this->KKT.U);
    auto t11 = std::chrono::high_resolution_clock::now();

    this->p_lambda = this->KKT.solve(this->KKT_rhs);
    auto t12 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Time solve LU: " << duration << endl;
    this->p = VectorXd(this->optimProblem.n_vars,1);

    for(int i=0; i<this->optimProblem.n_vars; i++)
    {
        this->p(i,0) = this->p_lambda(i,0);
    }

    //Initialize step size
    if(this->useLineSearchStepSizeMemory)
    {
        //cout << "useLineSearchStepSizeMemory" << endl;
        this->t_start = min(this->t*this->lineSearchStepSizeMemoryFactor,1);
    }
    else
    {
        this->t_start = 1;
    }

    if(this->useLineSearchStepSizeLimit)
    {
        MatrixXd tempY = this->y;
        MatrixXd tempP = this->p;
        matrix_reshape(tempY, tempY.rows()/this->optimProblem.dim, this->optimProblem.dim);
        matrix_reshape(tempP,tempP.rows()/this->optimProblem.dim, this->optimProblem.dim);
        //cout << "useLineSearchStepLimit" << endl;
        this->t = min(this->t_start, this->lineSearchStepSizeLimitFactor*this->optimProblem.getMaxStep(tempY, tempP));
    }
    else
    {
        this->t = this->t_start;
    }

    //Line search
    if(this->useLineSearch)
    {
        //cout << "useLineSearch" << endl;
        double linesearch_cond_lhs, linesearch_cond_rhs;
        computeLineSearchCond(linesearch_cond_lhs, linesearch_cond_rhs);

        while(linesearch_cond_lhs > linesearch_cond_rhs)
        {
          //  cout << linesearch_cond_lhs << " > " << linesearch_cond_rhs << endl;
            this->t = this->ls_beta*this->t;
            computeLineSearchCond(linesearch_cond_lhs, linesearch_cond_rhs);
        }
    }

    this->x_prev = this->x;
    this->x = this->y+this->t*this->p;

    //Update values

}

double OptimSolverAcclQuadProx::min(double value1, double value2)
{
    if(value1 < value2)
        return value1;
    else
        return value2;
}

void OptimSolverAcclQuadProx::computeLineSearchCond(double &linesearch_cond_lhs, double &linesearch_cond_rhs)
{
    MatrixXd tempY = this->y;
    VectorXd tempP = this->p;
    for(int i=0; i<tempP.rows(); i++)
    {
        for(int j=0; j<tempP.cols(); j++)
        {
            tempP(i,j) = tempP(i,j)*this->t;
        }
    }

//Add elemento from TempY and TempP
    for(int i=0; i<tempY.rows(); i++)
    {
        for(int j=0; j<tempY.cols(); j++)
        {
            tempY(i,j) = tempY(i,j)+tempP(i,j);
        }
    }

//    matrix_reshape(tempP, tempP.rows()/this->optimProblem.dim, this->optimProblem.dim);
    //cout << "Evaluate Value" << endl;
//    cout << "->> evaluateValue" << endl;
    this->optimProblem.evaluateValue(tempY, linesearch_cond_lhs);
    this->f_count++;
    linesearch_cond_rhs = this->y_f+this->ls_alpha*this->t*this->y_fgrad.transpose()*tempP;
}

#endif