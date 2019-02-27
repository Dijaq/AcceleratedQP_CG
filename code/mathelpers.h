#ifndef MATHELPERS_H
#define MATHELPERS_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <math.h>
#include "utils.h"
#include <chrono>

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

	for(int i=0; i<H.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(H, i); it; ++it)
		{
			KKT_mat.insert(it.row(), it.col()) = it.value();
			/*cout << "value: " << it.value() << endl;
			cout << "row: "<< it.row() << endl;
			cout << "col: "<< it.col() << endl;
			cout << it.index() << endl;*/
		}
	}

	for(int i=0; i<eq_lhs_transpose.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(eq_lhs_transpose, i); it; ++it)
		{
			KKT_mat.insert(it.row(), it.col()+H.cols()) = it.value();
		}
	}

	for(int i=0; i<eq_lhs.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(eq_lhs, i); it; ++it)
		{
			KKT_mat.insert(it.row()+H.rows(), it.col()) = it.value();
		}
	}

	for(int i=0; i<sparse.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(sparse, i); it; ++it)
		{
			KKT_mat.insert(it.row()+H.rows(), it.col()+H.cols()) = it.value();
		}
	}

	/*int j=0;
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
	}*/

	return KKT_mat;

}

SparseMatrix<double> join_matrices_2x1(SparseMatrix<double> A, SparseMatrix<double> B)
{
	SparseMatrix<double> KKT_mat(A.rows()+B.rows(),A.cols());

	for(int i=0; i<A.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
		{
			KKT_mat.insert(it.row(), it.col()) = it.value();
		}
	}

	for(int i=0; i<B.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(B, i); it; ++it)
		{
			KKT_mat.insert(it.row()+A.rows(), it.col()) = it.value();
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

VectorXd solve_Lx(MatrixXd &L, VectorXd X)
{
	for(int k=0; k<L.rows(); k++)
	{
		//cout << smatrix << endl;
		for(int i=k+1; i<L.rows(); i++)
		{
			float first = L(i,k);
			if(L(i,k) != 0)
			{
				L.row(i) = L.row(i)-(first/L(k,k))*L.row(k);
			}

			X(i,0) = X(i,0)-(first/L(k,k))*X(k,0);
		}
	}
	return X;
}

VectorXd solve_Ux(MatrixXd &U, VectorXd X)
{
	for(int k=U.rows()-1; k>=0; k--)
	{
		//cout << smatrix << endl;
		for(int i=k-1; i>=0; i--)
		{
			float first = U(i,k);
			if(U(i,k) != 0)
			{
				U.row(i) = U.row(i)-(first/U(k,k))*U.row(k);
			}

			X(i,0) = X(i,0)-(first/U(k,k))*X(k,0);
		}
	}

	for(int i=0; i< U.rows(); i++)
	{
		float div = U(i,i);
		U(i,i) /= div;
		X(i,0) /= div;
	}
	return X;
}

MatrixXd solve_Ax(SparseMatrix<double> sA, SparseMatrix<double> sb)
{
	MatrixXd A = sA;
	MatrixXd b = sb;

	//A.colPivHouseholderQr().solve(b);
	cout << "A*A'" << endl;
	auto t11 = std::chrono::high_resolution_clock::now();
	MatrixXd A_t = A.transpose()*A;
    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    cout << "A_t: " << duration << endl;
	VectorXd B_t = A.transpose()*b;

	cout << "Solve LX" << endl;
	solve_Lx(A_t, B_t);
	cout << "End LX" << endl;
	//(A.transpose() * A).ldlt().solve(A.transpose() * b);

	print_dimensions("A: ", A);
	print_dimensions("b: ", b);
}

void solveConstrainedLS(SparseMatrix<double> T,MatrixXd R, SparseMatrix<double> eq_lhs, MatrixXd eq_rhs)
{
	int n_vars = eq_lhs.cols();
	int n_eq = eq_lhs.rows();
	
	SparseMatrix<double> sR = R.sparseView();
	SparseMatrix<double> sEq_rhs = eq_rhs.sparseView();
	/*SparseMatrix<double> sT = T;
	SparseMatrix<double> sEq_lhs = eq_lhs;*/
	SparseMatrix<double> sp(n_eq, n_eq);

	solve_Ax(join_matrices((T.transpose()*T), sR.transpose(), eq_lhs,sp),
		join_matrices_2x1(T.transpose()*sR, sEq_rhs));
	//T.transpose()*sR;
	//Solve Matrix
}

#endif // MATHELPERS
