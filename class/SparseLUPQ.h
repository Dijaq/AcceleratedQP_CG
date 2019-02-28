#ifndef SPARSE_LUPQ_H
#define SPARSE_LUPQ_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>

class SparseLUPQ{
public:

	Eigen::SparseLU<SparseMatrix<double>> LU;

	SparseLUPQ();
	SparseLUPQ(SparseMatrix<double> A);
//	VectorXd solve(VectorXd LHS);

};

SparseLUPQ::SparseLUPQ(){}

SparseLUPQ::SparseLUPQ(SparseMatrix<double> A)
{
	LU.analyzePattern(A);
	LU.factorize(A);

}

/*VectorXd SparseLUPQ::solve(VectorXd LHS)
{
	VectorXd b = LU.rowsPermutation()*LHS;
	LU.matrixL().solveInPlace(b);
	LU.matrixU().solveInPlace(b);

	return b;
}*/


#endif