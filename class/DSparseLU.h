#ifndef DSPARSE_LU_H
#define DSPARSE_LU_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

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
};

DSparseLU::DSparseLU(){}

DSparseLU::DSparseLU(MatrixXd smatrix)
{
	this->U = smatrix;
	MatrixXd ones = matriz_diagonal_ones(smatrix.rows(), smatrix.cols());
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
			int first = this->U(i,k);
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