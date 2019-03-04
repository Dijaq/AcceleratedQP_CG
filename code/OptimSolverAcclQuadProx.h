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
//#include "../class/Splu.h"
#include "../class/Param_State.h"
#include "../mex/computeMeshTranformationCoeffsMex.h"
#include "OptimSolverIterative.h"

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimSolverAcclQuadProx : public OptimSolverIterative
{
public:

    //SparseLU<SparseMatrix<double>> AQPLU;
    //Parameters
    MatrixXd x_prev;
    MatrixXd y;
    VectorXd p;
    //SparseMatrix<double> KKT;
    DSparseLU KKT_Class;
    //Splu KKT;
    SparseMatrix<double> KKT_mat;
    MatrixXd KKT_rhs;
    VectorXd p_lambda;
    double y_f = 0;
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

    //SparseMatrix<double> KKT_mat;

   
    //export_sparsemat_to_excel(optimProblem.H, "H_spars");
    //export_sparsemat_to_excel(this->optimProblem.H, "PH_spars");

    if(this->useQuadProxy)
    {
        SparseMatrix<double> spar(optimProblem.n_eq, optimProblem.n_eq);
        //auto t11 = std::chrono::high_resolution_clock::now();
        export_sparsemat_to_excel(optimProblem.eq_lhs, "ValidarDatos/eq_lhs_load");
        this->KKT_mat = join_matrices(optimProblem.H, optimProblem.eq_lhs.transpose(), optimProblem.eq_lhs, spar);
        MatrixXd m = this->KKT_mat;
        export_mat_to_excel(m, "ValidarDatos/KKT_matmm");
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
        this->KKT_mat = join_matrices(ones, optimProblem.eq_lhs.transpose(), optimProblem.eq_lhs, spar);
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

    /*SparseMatrix<double> pract(3,3);
    pract.insert(0,0) = 3;
    pract.insert(1,1) = 7;
    pract.insert(2,0) = 1;
    pract.insert(2,1) = 2;
    pract.insert(2,2) = 1;*/


    MatrixXd matrix = this->KKT_mat;
    export_mat_to_excel(matrix, "c_KKT_matm");
    auto t11 = std::chrono::high_resolution_clock::now();
    //SparseLUPQ sparseLUPQ(KKT_mat);
    //this->KKT = SparseLUPQ(KKT_mat);
    //DSparseLU(matrix);
    cout << "init matrix" << endl;
    MatrixXd m = MatrixXd::Zero(4,4);
    m(0,1) = 1;
    m(0,2) = 3;
    m(0,3) = 4;
    m(1,0) = 3;
    m(1,2) = 2;
    m(2,0) = 1;
    m(2,1) = 1;
    m(3,1) = 4;
    m(3,2) = 4;
    //this->KKT_Class = DSparseLU(m, true);
    this->KKT_Class = DSparseLU(matrix, true);
    //export_mat_to_excel(this->KKT_Class.P.transpose()*this->KKT_Class.L*this->KKT_Class.U, "c_LU");
    //this->KKT = DSparseLU(matrix);
    //this->KKT = Splu(KKT_mat);
    //cout << "***" << this->KKT.LU.matrixL().rows() << endl;

    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "Time AQP2: " << duration1 << endl;
    
    //cout << "1****" << (*this->KKT.LU).matrixL().rows() << endl;

    //init internal variables
    //MatrixXd initZeros(this->optimProblem.n_vars+this->optimProblem.n_eq,1);
    //this->KKT_rhs = initZeros;
    this->KKT_rhs = MatrixXd::Zero(this->optimProblem.n_vars+this->optimProblem.n_eq,1);
    //cout << "2*****" << (*this->KKT.LU).matrixL().rows() << endl;

    this->x_prev = this->x;
    //cout << "******" << (*this->KKT.LU).matrixL().rows() << endl;
    this->t = 1/this->lineSearchStepSizeMemoryFactor;
    //Store init time
    //t_init = toc(t_init_start)
    //cout << "*******" << (*this->KKT.LU).matrixL().rows() << endl;
}

OptimSolverAcclQuadProx::~OptimSolverAcclQuadProx(){}

void OptimSolverAcclQuadProx::setKappa(float kappa)
{
    this->theta = (1-sqrt(1/kappa))/(1+sqrt(1/kappa));
}

void OptimSolverAcclQuadProx::solveTol(float TolX, float TolFun, int max_iter)
{
    /*Secion of create the matriz L and U
    */
    //cout << "varsssss: "<< this->optimProblem.n_vars << endl;
    this->p = VectorXd::Zero(this->optimProblem.n_vars,1);
    this->y_fgrad = VectorXd::Zero(this->optimProblem.n_vars,1);

    export_sparsemat_to_excel(this->KKT_mat, "KKT_mat");
    SparseLU<SparseMatrix<double>> LU;
    LU.analyzePattern(this->KKT_mat);
    LU.factorize(this->KKT_mat);

/*    cout << "-->"<<LU.matrixL().cols() << endl;
    cout << "-->"<<LU.matrixU().cols() << endl;
    cout << "-->"<<LU.rowsPermutation().rows() << endl;
    cout << "-->"<<LU.colsPermutation().rows() << endl;
*/

    //This is only logs is not necesarialy for the implementation
    //logInit();
    //logState();

    //run num_iter iteration
    //for(int i=0; i<100; i++)
    for(int i=0; i<max_iter; i++)
    {
        auto t11 = std::chrono::high_resolution_clock::now();
        //this is the time to start no necesarialy
        //t_iter_start = tic

        //iterate(LU);
        //Falta Stop criteria
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

        export_mat_to_excel(this->y, "ValidarDatos/c_y"+to_string(i));
        export_mat_to_excel(this->x, "ValidarDatos/c_x"+to_string(i));
        export_mat_to_excel(this->x_prev, "ValidarDatos/c_x_prev"+to_string(i));
        //cout << "y_fff: " << this->y_f << endl;

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

        export_mat_to_excel(this->y_fgrad, "ValidarDatos/c_y_fgrad"+to_string(i));
        //export_mat_to_excel(this->y_fgrad, "y_fgrad");
        for(int i=0; i<this->optimProblem.n_vars; i++)
        {
            this->KKT_rhs(i,0) = -this->y_fgrad(i,0); 
            
        }

        export_mat_to_excel(this->KKT_rhs, "ValidarDatos/c_KKT_rhs"+to_string(i));

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
        //print_dimensions("KKT_rhs", this->KKT_rhs);

        //export_mat_to_excel(this->KKT_rhs, "KKT_rhs");

        /*VectorXd b = (LU.rowsPermutation().transpose())*this->KKT_rhs;
        export_mat_to_excel(b, "b_perR");
        
        LU.matrixL().solveInPlace(b);
        LU.matrixU().solveInPlace(b);
        this->p_lambda = (LU.colsPermutation().transpose())*b;*/


        //this->p_lambda = b;

//        VectorXd prueba_lambda = this->KKT_Class.solve(this->KKT_rhs);
        this->p_lambda = this->KKT_Class.solve(this->KKT_rhs);
        //export_mat_to_excel(prueba_lambda, "prueba_lambda");
        export_mat_to_excel(this->p_lambda, "ValidarDatos/c_p_lambda"+to_string(i));

        for(int i=0; i<this->optimProblem.n_vars; i++)
        {
            this->p(i,0) = this->p_lambda(i,0);
        }

        export_mat_to_excel(this->p, "ValidarDatos/c_p"+to_string(i+1));

        //Initialize step size
        cout << "t1: " << this->t << endl;
        if(this->useLineSearchStepSizeMemory)
        {
            //cout << "useLineSearchStepSizeMemory" << endl;
            this->t_start = min(this->t*this->lineSearchStepSizeMemoryFactor,1);
        }
        else
        {
            this->t_start = 1;
        }

        cout << "t_start: " << this->t_start << endl;

        //cout << "start: " << this->t_start << endl;

        if(this->useLineSearchStepSizeLimit)
        {
            //export_sparsemat_to_excel(this->y.sparseView(), "1y");
            //export_mat_to_excel(this->p, "1p");
            MatrixXd tempY = this->y;
            MatrixXd tempP = this->p;
            matrix_reshape(tempY, tempY.rows()/this->optimProblem.dim, this->optimProblem.dim);
            matrix_reshape(tempP,tempP.rows()/this->optimProblem.dim, this->optimProblem.dim);
            //cout << "useLineSearchStepLimit" << endl;
            //cout << "getstep: " << this->optimProblem.getMaxStep(tempY, tempP) << endl;
            cout << "min: " << this->optimProblem.getMaxStep(tempY, tempP) << endl;
            this->t = min(this->t_start, this->lineSearchStepSizeLimitFactor*this->optimProblem.getMaxStep(tempY, tempP));
        }
        else
        {
            this->t = this->t_start;
        }

        cout << "t2: " << this->t << endl;

        //Line search
        if(this->useLineSearch)
        {
            //cout << "useLineSearch" << endl;
            double linesearch_cond_lhs, linesearch_cond_rhs;
            computeLineSearchCond(linesearch_cond_lhs, linesearch_cond_rhs);

            cout << "---------->"<<linesearch_cond_lhs/100000 << " - " << linesearch_cond_rhs/100000 << endl;

            while(linesearch_cond_lhs > linesearch_cond_rhs)
            {
                cout << linesearch_cond_lhs << " <------www------> " << linesearch_cond_rhs << endl;
                //cout << "beta: " << this->ls_beta<< endl;
                this->t = this->ls_beta*this->t;
                computeLineSearchCond(linesearch_cond_lhs, linesearch_cond_rhs);
            }
        }
        cout << "t3: " << this->t << endl;

        this->x_prev = this->x;
        this->x = this->y+this->t*this->p;

        auto t12 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
        cout << "Iteration: " << i << " time: " <<duration << endl;

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
    print_dimensions("KKT_rhs", this->KKT_rhs);

    //this->p_lambda = this->KKT.solve(this->KKT_rhs);
    
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
    linesearch_cond_rhs = this->y_f+this->ls_alpha*this->t*this->y_fgrad.transpose()*this->p;
}

#endif