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
    int x_prev;
    int y;
    int p;
    SparseMatrix<double> KKT;
    int KKT_rhs;
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
    int theta; // = [];

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
        KKT_mat = join_matrices(optimProblem.H, optimProblem.eq_lhs.transpose(), optimProblem.eq_lhs, spar);
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

    //this->KKT = SparseLU(KKT_mat);

    //init internal variables

    //Store init time
    //t_init = toc(t_init_start)
}

OptimSolverAcclQuadProx::~OptimSolverAcclQuadProx()
{}

#endif