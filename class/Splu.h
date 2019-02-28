#ifndef SPLU_H
#define SPLU_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>
#include "../code/utils.h"

class Splu{
public:

	SparseLU<SparseMatrix<double>> *LU;

//	Splu();
	Splu(SparseMatrix<double> A);
	VectorXd solve(VectorXd LHS);

};

/*Splu::Splu(){
	//this->LU = 
}
*/
Splu::Splu(SparseMatrix<double> A)
{
	//cout << "<A>: " << A.rows() << " - " << A.cols() << endl;
	SparseLU<SparseMatrix<double>> solveLU;
	solveLU.analyzePattern(A);
	solveLU.factorize(A);
	//this->LU = solveLU;
	//cout << "**" << solveLU.matrixL().rows() << " - "<< solveLU.matrixL().cols() <<endl;
}

VectorXd Splu::solve(VectorXd LHS)
{
	print_dimensions("LHS: ", LHS);
	cout << "*" << (*this->LU).matrixL().rows() << endl;
	/*cout << "row: " <<  (*LU).rowsPermutation().cols() << endl;
	VectorXd b = ((*LU).rowsPermutation()).transpose()*LHS;*/
	VectorXd b = LHS;
	(*LU).matrixL().solveInPlace(b);
	(*LU).matrixU().solveInPlace(b);

	return b;
}


#endif