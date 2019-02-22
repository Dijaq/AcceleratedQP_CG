#ifndef MATHELPERS_H
#define MATHELPERS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

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

#endif // MATHELPERS
