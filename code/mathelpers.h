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

//Join matrices libe [matriz1 matriz2; matriz3 matriz4]
SparseMatrix<double> join_matrices(SparseMatrix<double> H, SparseMatrix<double> eq_lhs_transpose, SparseMatrix<double> eq_lhs, SparseMatrix<double> sparse)
{
	SparseMatrix<double> KKT_mat(H.rows()+eq_lhs.rows(),H.cols()+eq_lhs_transpose.cols());

	int j=0;
	for(int i=0; i<H.rows(); i++)
	{
		for(j=0; j<H.cols(); j++)
		{
			if(H.coeffRef(i,j) != 0)
				KKT_mat.insert(i,j) = H.coeffRef(i,j);
		}
	}

	for(int i=0; i<eq_lhs_transpose.rows(); i++)
	{
		for(int j=0; j<eq_lhs_transpose.cols(); j++)
		{
			if(eq_lhs_transpose.coeffRef(i,j) != 0)
				KKT_mat.insert(i,j+H.cols()) = eq_lhs_transpose.coeffRef(i,j);
		}
	}

	for(int i=0; i<eq_lhs.rows(); i++)
	{
		for(j=0; j<eq_lhs.cols(); j++)
		{
			if(eq_lhs.coeffRef(i,j) != 0)
				KKT_mat.insert(i+H.rows(),j) = eq_lhs.coeffRef(i,j);
		}
	}

	for(int i=0; i<sparse.rows(); i++)
	{
		for(int j=0; j<sparse.cols(); j++)
		{
			if(sparse.coeffRef(i,j) != 0)
				KKT_mat.insert(i+H.rows(),j+H.cols()) = sparse.coeffRef(i,j);
			//sparse(i,j);
		}
	}

	return KKT_mat;

}

SparseMatrix<double> create_SparseMatrix_ones(int rows, int cols)
{
	SparseMatrix<double> ones(rows, cols);
	int menor = 0;
	if(rows < cols)
		menor = rows;
	else
		menor = cols;
	for(int i=0; i<menor; i++)
	{
		ones.insert(i,i) = 1;
	}

	return ones;
}

#endif // MATHELPERS
