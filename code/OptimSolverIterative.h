#ifndef OPTIM_SOLVER_ITERATIVE_H
#define OPTIM_SOLVER_ITERATIVE_H

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "OptimSolver.h"

using namespace std;
using namespace cv;
using namespace Eigen;

class OptimSolverIterative : public OptimSolver
{
public:

	//value logs
	int iter;
	int t_iter;

	//solver vars
	//int t = nan;
	float t;

	int stopCntAccept = 5;
	int tolXCnt = 0;
	int tolFunCnt = 0;

	OptimSolverIterative();
	~OptimSolverIterative();
	void logInit();
	void logState();
};

OptimSolverIterative::OptimSolverIterative()
{}

OptimSolverIterative::~OptimSolverIterative()
{}


void OptimSolverIterative::logInit()
{
	/*this->iter = 0;
	this->t_iter = 0;
	this-tag = */
}

void OptimSolverIterative::logState()
{

}

#endif