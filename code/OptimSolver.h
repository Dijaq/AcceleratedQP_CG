#ifndef OPTIM_SOLVER_H
#define OPTIM_SOLVER_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "mathelpers.h"
#include "utils.h"
#include "OptimProblemIsoDist.h"

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimSolver
{
public:

	//properties
	OptimProblemIsoDist optimProblem;
	MatrixXd x; //This is a matrix of 1 column

	//Log
	//char tag = '';
	string tag = "";
	int verbose = 2;

	//default solver parameters
	int default_max_iter = 500;
	float default_TolX = 1e-10;
	float default_TolFun = 1e-6;



	OptimSolver();
	~OptimSolver();
	void initSolver(string tag, OptimProblemIsoDist optimProblem);
};

OptimSolver::OptimSolver()
{}

OptimSolver::~OptimSolver()
{}

void OptimSolver::initSolver(string tag, OptimProblemIsoDist optimProblem)
{
	//Copy
	this->tag = tag;
	this->optimProblem = optimProblem;
	this->x = colStack(optimProblem.x0);
	//Report no implemented

}

#endif