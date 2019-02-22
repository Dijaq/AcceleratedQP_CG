#ifndef OPTIMPROBLEM_H
#define OPTIMPROBLEM_H

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

class OptimProblem
{
public:
	SparseMatrix<double> T;
	MatrixXd eq_lhs;
	MatrixXd eq_rhs;
	int x0;
	int n_vars;
	int n_eq;
	MatrixXd H;

	//Parameters
	int verbose = 2;

	OptimProblem();
	~OptimProblem();
	void initProblem();
};

OptimProblem::OptimProblem()
{}

void OptimProblem::initProblem()
{
	this->n_vars = this->T.cols();
	this->n_eq = this->eq_lhs.rows();
}

OptimProblem::~OptimProblem()
{}

#endif