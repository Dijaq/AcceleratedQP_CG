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
#include <limits>
#include "utils.h"
#include "../class/Param_State.h"
#include "../mex/computeMeshTranformationCoeffsMex.h"
#include "../mex/computeFunctionalIsoDistMex.h"
#include "../mex/computeInjectiveStepSizeMex.h"
#include "../mex/projectRotationMexFast.h"

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
    VectorXd R;
    double f_val;
    VectorXd Tx_grad;
    bool flips;
    int localHess;

    //Parameters
    float maxStepTol = 1e-8;
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
    void setQuadraticProxy();
    double getMaxStep(MatrixXd x, MatrixXd p);
    void evaluateValue(MatrixXd x, double &f);
    void evaluateGrad(MatrixXd x, double &f, VectorXd &f_grad);
    void evaluateValueGrad(MatrixXd x, double &f, VectorXd &f_grad);
    void evaluateFunctional(MatrixXd x, bool doVal, bool doGrad, bool doHess, double &f, VectorXd &f_grad);
    MatrixXd sparseMul(SparseMatrix<double> A, MatrixXd B);
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
    if(mesh.V.cols() == 2)
        this->V = mesh.V;
    else
        this->V = V0;
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
    cout << "first step" << endl;
    computeMeshTranformationCoeffsFullDim(this->F, this->V, this->T, this->areas);//Finished
    //set initial configuration
    cout << "second step" << endl;
    initVertices(V0);//Falta
    //set quadratic proxy
    cout << "third step" << endl;
    setQuadraticProxy();//Finished
    //Finish construction
    cout << "four step" << endl;
    initProblem();//Finished
}

void OptimProblemIsoDist::initVertices(MatrixXd v0)
{ 
    print_dimensions("V0: ", v0);
    //MatrixXd x0;
    if((v0.rows() != 0) || (v0.cols() != 0))
    {
        cout << "no es empty" << endl;
        this->x0 = v0;
    }
    else
    {
        //Here in our example
        this->x0 = colStack(this->V);
 //       export_mat_to_excel(this->x0, "stackx0");
        //SparseMatrix<double> x0 = this->x0.sparseView();
        /*Falta implementar revisar el codigo fuente*/
        //cout << "Tx: " << this->T.adjoint().rows() << " - " << this->T.adjoint().cols() << endl;
        //cout << "Tx0: " << x0.rows() << " - " << x0.cols() << endl;
        //export_mat_to_excel(sparseMul(this->T, x0), "mul");
        sparseMul(this->T,x0);
        for(int arapIter = 0; arapIter<this->initArapIter; arapIter++)
        {

            //this->Tx = this->T*this->x0;
            this->Tx = this->T*this->x0;
            //export_mat_to_excel(this->Tx, "Tx"+to_string(arapIter)); 
            this->R = this->Tx;

            if(this->dim == 2)
                projBlockRotation2x2(this->R, this->dim);
            else
                projBlockRotation3x3(this->R, this->dim);
            //this->x0 = solveConstrainedLS(this->T, this->R, this->eq_lhs, this->eq_rhs);

            /*if(arapIter == 0)
            {
                print_dimensions("T: ", this->T);
                print_dimensions("R: ", this->R);
                print_dimensions("EQ_LHS: ", this->eq_lhs);
                print_dimensions("EQ_RHS: ", this->eq_rhs);
            }*/
            //This is too low, falta mejorar
            this->x0 = solveConstrainedLS(this->T, this->R, this->eq_lhs, this->eq_rhs);
        }


        matrix_reshape(this->x0, this->n_vert, this->dim);
        //cout << "size: "<<this->x0.rows() << "-"<<this->x0.cols() << endl;
    }

    print_dimensions("T: ", this->T);

/*    cout << this->T.rows() << "-" << this->T.cols() << endl;
    cout << x0.rows() << "-" << x0.cols() << endl;*/

    //cout << "->"<<(this->T).rows() <<" - " <<(this->T).cols()<<endl;
    this->Tx = this->T*colStack(this->x0);

    cout << "dim: " << this->dim << endl;

    double val;
    if(this->dim == 2)
        helperFunctionalIsoDist2x2(this->Tx, this->areas, this->dim, val, this->flips);
    else
        helperFunctionalIsoDist3x3(this->Tx, this->areas, this->dim, val, this->flips);
    //Fix if there are flips
    if(this->flips)
    {
        cout << "Hay flips" << endl;
        //Project onto BD
        /*Falta implementar SOLVER PROJECT BD*/
    }
    else
    {
        //No hace nada
        cout << "No hay flips" << endl;
    }

    //ProgBar is only for matlab


    //Check if the solution is orientation preserving
    //this->Tx = this->T*colStack(x0);
}

MatrixXd OptimProblemIsoDist::sparseMul(SparseMatrix<double> A, MatrixXd B)
{
    MatrixXd sparseMul = MatrixXd::Zero(A.rows(), 1);

    for(int i=0; i< A.outerSize(); i++)
    {
        for(SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
        {
            sparseMul(it.row(),0) += it.value()*B(it.col(),0); 
            /*if(it.row() == 1)
            {
                cout << "-->"<< "("<< it.row()<<","<< it.col()<<")"<< it.value() << " : " << B(it.col(),0) <<" = "<< it.value()*B(it.col(),0) << endl;//sparseMul(it.row(),0)
            }*/
            //sparseMul.insert(it.row(), it.col()) = it.value();
        }
    }

    return sparseMul;

}

void OptimProblemIsoDist::setQuadraticProxy()
{
    //export_sparsemat_to_excel(this->T, "T");
    //cout << "T: " << this->T.rows() << endl;
    //cout << "T: " << this->T.cols() << endl;
    VectorXd ar = this->areas;
    VectorXd ones = VectorXd::Ones(pow(this->dim,2),1);
    VectorXd_sqrt(ar);

    //export_mat_to_excel(kron(ar, ones), "kron_areas");

    SparseMatrix<double> wT = spdiag(kron(ar, ones))*this->T;
    this->H = 2*(wT.transpose()*wT);
    //export_sparsemat_to_excel(H, "H");
    //cout << "k: "<<this->H.cols() << endl;

}

double OptimProblemIsoDist::getMaxStep(MatrixXd x, MatrixXd p)
{

    double t_max;
    //cout <<"x: " <<x.cols()<<endl;
    if(this->dim == 2)
        computeInjectiveStepSize_2d(this->F, x, p, this->maxStepTol, t_max);
    else
        computeInjectiveStepSize_3d(this->F, x, p, this->maxStepTol, t_max);
    return t_max;
}

void OptimProblemIsoDist::evaluateValue(MatrixXd x, double &f)
{
    VectorXd f_grad;
    evaluateFunctional(x, true, false, false, f, f_grad);
}

void OptimProblemIsoDist::evaluateGrad(MatrixXd x, double &f, VectorXd &f_grad)
{
    evaluateFunctional(x, false, true, false, f, f_grad);
}

void OptimProblemIsoDist::evaluateValueGrad(MatrixXd x, double &f, VectorXd &f_grad)
{
    evaluateFunctional(x, true, true, false, f, f_grad);
}

void OptimProblemIsoDist::evaluateFunctional(MatrixXd x, bool doVal, bool doGrad, bool doHess, double &f, VectorXd &f_grad)
{
    //print_dimensions("T: ", this->T);
    //print_dimensions("x: ", x);
    this->Tx = this->T*x;
    //export_mat_to_excel(this->Tx, "Tx_funn");
    //export_mat_to_excel(x, "x_funn");// regular
    this->Tx_grad = this->Tx;
    //export_mat_to_excel(Tx_grad, "Tx");//Este esta medio raro
    if(doVal || doGrad)
    {
        //cout << "Inpute do Val" << endl;
        if(this->dim == 2)
            helperFunctionalIsoDist2x2(this->Tx_grad, this->areas, this->dim, this->f_val, this->flips);
        else
            helperFunctionalIsoDist3x3(this->Tx_grad, this->areas, this->dim, this->f_val, this->flips);

        if(this->flips)
        {
            this->f_val = numeric_limits<int>::max();
        }
    }

    if(doVal)
    {
        f = f_val;
    }

    if(doGrad)
    {
        //export_mat_to_excel(this->Tx_grad, "Tx_gradd");
        //print_dimensions("Tx->", this->Tx_grad);
        //print_dimensions("T->", this->T);
        f_grad = (this->Tx_grad.transpose()*this->T).transpose();
    }

    if(doHess)
    {
        cout << "The function computeHessianIsoDistMex is not implement by the author" << endl;
        /*n_arg = n_arg + 1;
        obj.localHess = computeHessianIsoDistMex(obj.Tx, obj.areas, obj.dim);
        varargout{n_arg} = obj.T'*obj.localHess*obj.T;*/
    }

    /*
    % return
    n_arg = 0;
    if doVal
        n_arg = n_arg + 1;
        varargout{n_arg} = obj.f_val;
    end
    if doGrad
        n_arg = n_arg + 1;
        varargout{n_arg} = (obj.Tx_grad'*obj.T)';
    end
    if doHess
        n_arg = n_arg + 1;
        obj.localHess = computeHessianIsoDistMex(obj.Tx, obj.areas, obj.dim);
        varargout{n_arg} = obj.T'*obj.localHess*obj.T;
    end
    */
}

#endif