#ifndef MATHELPERS_H
#define MATHELPERS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

using namespace std;
using namespace Eigen;

MatrixXd colStack(MatrixXd x)
{
	MatrixXd y(x.rows()*x.cols(), 1);
	for(int j=0; j<x.cols(); j++)
	{
		for(int i=0; i<x.rows(); i++)
		{
			y(j*x.rows()+i,0) = x(i,j);
		}
	}
	//y = ;

	return y;
}

#endif // MATHELPERS
