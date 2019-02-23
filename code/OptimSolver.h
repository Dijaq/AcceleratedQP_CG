#ifndef OPTIM_SOLVER_H
#define OPTIM_SOLVER_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimSolver
{
public:

	//properties
	int optimProblem;
	int x;

	//Log
	char tag = '';
	int verbose = 2;

	//default solver parameters
	int default_max_iter = 500;
	float default_TolX = 1e-10;
	float default_TolFun = 1e-6;



	OptimSolver();
	~OptimSolver();
};

OptimSolver::OptimSolver()
{}

OptimSolver::~OptimSolver()
{}

#endif