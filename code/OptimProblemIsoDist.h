#ifndef OPTIMPROBLEMISODIST_H
#define OPTIMPROBLEMISODIST_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "OptimProblem.h"
#include "../class/Param_State.h"
#include "../mex/computeMeshTranformationCoeffsMex.h"

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimProblemIsoDist : public OptimProblem
{
public:

    //Parameter mesh
    MatrixXd V;
    MatrixXd F;
    int dim;
    int n_vert;
    int n_tri;
    VectorXd areas;

    //Internal Variables
    int Tx;
    int R;
    int f_val;
    int Tx_grad;
    int flips;
    int localHess;

    //Parameters
    float maxStepTol = 1e-10;
    int initArapIter;

    //bd projection parameters
    int proj_bd_K = 10;
    int proj_bd_lb = -1;
    int proj_bd_ub = -1;
    int proj_bd_iter_max = 1000;
    int proj_bd_tol_err = 1e-10;
    bool proj_bd_verbose = true;

    OptimProblemIsoDist();
    OptimProblemIsoDist(Param_State mesh, MatrixXd V0, int initArapIter);
    ~OptimProblemIsoDist();
};

OptimProblemIsoDist::OptimProblemIsoDist()
{
    /*this->maxStepTol = 1e-8;
    this->proj_bd_K = 10;
    this->proj_bd_lb = -1;
    this->proj_bd_ub = -1;
    this->proj_bd_iter_max = 1000;
    this->proj_bd_tol_err = 1e-10;
    this->proj_bd_verbose = true;*/
}

OptimProblemIsoDist::OptimProblemIsoDist(Param_State mesh, MatrixXd V0, int initArapIter)
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

    //Compute transformations
    computeMeshTranformationCoeffsFullDim(this->F, this->V, this->T, this->areas);
    //set initial configuration
}

OptimProblemIsoDist::~OptimProblemIsoDist()
{}

#endif