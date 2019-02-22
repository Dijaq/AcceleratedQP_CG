#ifndef OPTIMSOLVERACCLQUADPROX_H
#define OPTIMSOLVERACCLQUADPROX_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "../class/Param_State.h"
#include "../mex/computeMeshTranformationCoeffsMex.h"

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimSolverAcclQuadProx
{
public:

    //Parameters
    int x_prev;
    int y;
    int p;
    int KKT;
    int KKT_rhs;
    int p_lambda

    OptimSolverAcclQuadProx();
    OptimSolverAcclQuadProx(Param_State mesh, MatrixXd V0, int initArapIter);
    ~OptimSolverAcclQuadProx();
};

OptimSolverAcclQuadProx::OptimSolverAcclQuadProx()
{
    /*this->maxStepTol = 1e-8;
    this->proj_bd_K = 10;
    this->proj_bd_lb = -1;
    this->proj_bd_ub = -1;
    this->proj_bd_iter_max = 1000;
    this->proj_bd_tol_err = 1e-10;
    this->proj_bd_verbose = true;*/
}

OptimSolverAcclQuadProx::OptimSolverAcclQuadProx(Param_State mesh, MatrixXd V0, int initArapIter)
{
    this->V = mesh.V;
    this->F = mesh.F;
    this->eq_lhs = mesh.eq_lhs;
    this->eq_rhs = mesh.eq_rhs;
    this->dim = mesh.F.cols()-1;
    this->n_vert = mesh.V.rows();
    this->n_tri = mesh.F.rows();

    /*VectorXd areas(this->F.rows());
    this->areas = areas;*/

    //See the real implementation
    this->initArapIter = initArapIter;

    computeMeshTranformationCoeffsFullDim(this->F, this->V, this->T, this->areas);
}

OptimSolverAcclQuadProx::~OptimSolverAcclQuadProx()
{}

#endif