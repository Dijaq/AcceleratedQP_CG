#ifndef SPLU_H
#define SPLU_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>

class Splu{
public:

	Eigen::SparseLU<SparseMatrix<double>> *LU;

	Splu();
	Splu(SparseMatrix<double> A);
	VectorXd solve(VectorXd LHS);

};

Splu::Splu(){}

Splu::Splu(SparseMatrix<double> A)
{
	SparseLU<SparseMatrix<double>> solveLU;
	solveLU.analyzePattern(A);
	solveLU.factorize(A);
	LU = &solveLU;
}

VectorXd Splu::solve(VectorXd LHS)
{
	VectorXd b = (*LU).rowsPermutation()*LHS;
	(*LU).matrixL().solveInPlace(b);
	(*LU).matrixU().solveInPlace(b);

	return b;
}


#endif