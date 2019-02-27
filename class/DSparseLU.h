#ifndef DSPARSE_LU_H
#define DSPARSE_LU_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>

class DSparseLU{
public:

	Eigen::MatrixXd L;
	Eigen::MatrixXd U;
	Eigen::MatrixXd P;//This is a sparse matrix
	Eigen::MatrixXd Q;

	DSparseLU();
	DSparseLU(MatrixXd smatrix);
	DSparseLU(SparseMatrix<double> smatrix);
	MatrixXd matriz_diagonal_ones(int rows, int cols);
	VectorXd solve(VectorXd LHS);
	VectorXd solve_Lx(MatrixXd L, VectorXd X);
	VectorXd solve_Ux(MatrixXd U, VectorXd X);
};

DSparseLU::DSparseLU(){}

DSparseLU::DSparseLU(MatrixXd smatrix)
{
	this->U = smatrix;
	
	auto t11 = std::chrono::high_resolution_clock::now();
	MatrixXd ones = matriz_diagonal_ones(smatrix.rows(), smatrix.cols());
    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();

    cout << "Create Identity: " << duration << endl;
	this->L = ones;
	this->P = ones;
	this->Q = ones;


	for(int k=0; k<this->U.rows(); k++)
	{
		//Pivote Parcial, Calculate change
		double mayor = abs(this->U(k,k));
		int fila = k;
		for(int i=k+1; i<this->U.rows(); i++)
		{
			if(mayor < abs(this->U(i,k)))
			{
				fila = i;
				mayor = abs(this->U(i,k));
			}
		}

		//Change in U and L
		for(int j=0; j<this->U.cols(); j++)
		{
			double tempU = this->U(k,j);
			this->U(k,j) = this->U(fila, j);
			this->U(fila, j) = tempU;

			double tempP = this->P(k,j);
			this->P(k,j) = this->P(fila, j);
			this->P(fila, j) = tempP;

			if(k > j)
			{
				double tempL = this->L(k,j);
				this->L(k,j) = this->L(fila, j);
				this->L(fila, j) = tempL;
			}
		}

		//cout << smatrix << endl;
		for(int i=k+1; i<this->U.rows(); i++)
		{
			float first = this->U(i,k);
			if(this->U(i,k) != 0)
				for(int j=0; j<this->U.cols(); j++)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					this->U(i,j) = this->U(i,j)-(first/this->U(k,k))*this->U(k,j);
					
				}
			this->L(i,k) = (first/this->U(k,k));
		}
	}

}

DSparseLU::DSparseLU(SparseMatrix<double> smatrix)
{
	for(int i=0; i<smatrix.outerSize(); i++)
	{
		for(SparseMatrix<double>::InnerIterator it(smatrix, i); it; ++it)
		{
			cout << it.row() << " - "<< it.col() << " value: " << it.value() << endl;
		}
	}
}

VectorXd DSparseLU::solve(VectorXd LHS)
{
	return solve_Ux(this->U,solve_Lx(this->L, (this->P.transpose()*LHS)));
}

VectorXd DSparseLU::solve_Lx(MatrixXd L, VectorXd X)
{
	for(int k=0; k<L.rows(); k++)
	{
		//cout << smatrix << endl;
		for(int i=k+1; i<L.rows(); i++)
		{
			float first = L(i,k);
			if(L(i,k) != 0)
			{
				for(int j=0; j<i; j++)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					L(i,j) = L(i,j)-(first/L(k,k))*L(k,j);
				}
			}

			//cout << "Num/Den: " << first << "/" << L(k,k) << " val: " << X(i,0) << endl;
			X(i,0) = X(i,0)-(first/L(k,k))*X(k,0);
		}
	}
	//cout << "L: " <<L << endl;
	//cout << "X: " << X << endl;
	return X;
}

VectorXd DSparseLU::solve_Ux(MatrixXd U, VectorXd X)
{
	for(int k=U.rows()-1; k>=0; k--)
	{
		//cout << smatrix << endl;
		for(int i=k-1; i>=0; i--)
		{
			float first = U(i,k);
			if(U(i,k) != 0)
				for(int j=U.cols()-1; j>i; j--)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					U(i,j) = U(i,j)-(first/U(k,k))*U(k,j);
				}

			//cout << "Num/Den: " << first << "/" << L(k,k) << " val: " << X(i,0) << endl;
			X(i,0) = X(i,0)-(first/U(k,k))*X(k,0);
		}
	}

	for(int i=0; i< U.rows(); i++)
	{
		float div = U(i,i);
		U(i,i) /= div;
		X(i,0) /= div;
	}
	//cout << "U: " <<U << endl;
	//cout << "X: " << X << endl;
	return X;
}


MatrixXd DSparseLU::matriz_diagonal_ones(int rows, int cols)
{
	MatrixXd ones(rows, cols);
	/*int menor = rows;
	if(cols < menor)
		menor = cols;

	for(int i=0; i<menor; i++)
	{
		ones(i,i) = 1;
	}*/

	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
		{
			if(i==j)
				ones(i,j) = 1;
			else
				ones(i,j) = 0;
		}
	}

	return ones;
}

#endif