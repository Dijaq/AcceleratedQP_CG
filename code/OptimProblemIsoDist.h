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
#include "mathelpers.h"
#include "../class/Param_State.h"
#include "../mex/computeMeshTranformationCoeffsMex.h"
#include "../mex/computeFunctionalIsoDistMex.h"

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
    VectorXd  Tx;
    int R;
    int f_val;
    int Tx_grad;
    bool flips;
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
    void initVertices(MatrixXd V0);
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

OptimProblemIsoDist::~OptimProblemIsoDist()
{}

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
    initVertices(V0);
}

void OptimProblemIsoDist::initVertices(MatrixXd v0)
{
    MatrixXd x0;
    if(!MatrixXd_isempty(v0))
    {
        x0 = v0;
    }
    else
    {
        //Here in our example
        /*Falta implementar revisar el codigo fuente*/
        x0 = colStack(this->V);
    }

/*    cout << this->T.rows() << "-" << this->T.cols() << endl;
    cout << x0.rows() << "-" << x0.cols() << endl;*/

    //cout << "->"<<(this->T*x0).rows() <<" - " <<(this->T*x0).cols()<<endl;
    this->Tx = this->T*x0;

    double val;
    helperFcuntionalIsoDist2x2(this->Tx, this->areas, this->dim, val, this->flips);
    //Fix if there are flips
    if(this->flips)
    {

    }

    //ProgBar is only for matlab


    //Check if the solution is orientation preserving
    //this->Tx = this->T*colStack(x0);
}

#endif