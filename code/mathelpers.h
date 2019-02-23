#ifndef MATHELPERS_H
#define MATHELPERS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <math.h>

using namespace std;
using namespace Eigen;

MatrixXd colStack(MatrixXd matrix)
{
	MatrixXd y(matrix.rows()*matrix.cols(), 1);
	for(int j=0; j<matrix.cols(); j++)
	{
		for(int i=0; i<matrix.rows(); i++)
		{
			y(j*matrix.rows()+i,0) = matrix(i,j);
		}
	}
	//y = ;

	return y;
}

bool MatrixXd_isempty(MatrixXd matrix)
{
	int sum =0;
	for(int i=0; i<matrix.rows(); i++)
	{
		for(int j=0; j<matrix.cols(); j++)
		{
			sum += matrix(i,j);
		}
	}

	if(sum == 0)
		return true;
	else 
		return false;
	//y = ;

}

void VectorXd_sqrt(VectorXd &vec)
{
	for(int i=0; i<vec.rows(); i++)
	{
		vec(i,0) = sqrt(vec(i,0));
	}
}

VectorXd kron(VectorXd first, VectorXd second)
{
	VectorXd result(first.rows()*second.rows(),1);
	for(int i=0; i<second.rows(); i++)
	{
		for(int j=0; j<first.rows(); j++)
		{
			result(j*(i+1),0) = first(j,0)*second(i,0);
		}
	}

	return result;
}

SparseMatrix<double> spdiag(VectorXd vect)
{
	int n = vect.rows();
	SparseMatrix<double> T(n,n);

	for(int i=0; i<n; i++)
	{
		T.insert(i,i) = vect(i,0);
	}

	return T;
}

#endif // MATHELPERS
