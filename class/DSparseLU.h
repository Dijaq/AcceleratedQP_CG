#ifndef DSPARSE_LU_H
#define DSPARSE_LU_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>

class DSparseLU{
public:

	Eigen::MatrixXd L;
	Eigen::MatrixXd U;
	Eigen::SparseMatrix<double> sL;
	Eigen::SparseMatrix<double> sU;
	Eigen::MatrixXd P;//This is a sparse matrix
	Eigen::MatrixXd Q;
	FullPivLU<MatrixXd> pivlu;
//	SparseLU<SparseMatrix<double>> LU;

	DSparseLU();
	DSparseLU(MatrixXd smatrix, bool s);
	DSparseLU(MatrixXd smatrix);
	DSparseLU(SparseMatrix<double> A);
	MatrixXd matriz_diagonal_ones(int rows, int cols);
	VectorXd solve(VectorXd LHS);
	VectorXd solve_Lx(MatrixXd L, VectorXd X);
	VectorXd solve_Ux(MatrixXd U, VectorXd X);
};

DSparseLU::DSparseLU(){}


/*LU implementation of Eigen*/
//Toma 8 veces mas tiempo que la implementacion propia
DSparseLU::DSparseLU(MatrixXd smatrix, bool s)
{
	pivlu = FullPivLU<MatrixXd>(smatrix);
	/*cout << "Original Matrix" << endl;
    cout << practice << endl;

    FullPivLU<MatrixXd> lu(practice);
    cout << lu.matrixLU() << endl;
    cout << "print matrix L " << endl;
    MatrixXd l = MatrixXd::Identity(3,3);
    l.block<3,3>(0,0).triangularView<StrictlyLower>() = lu.matrixLU();
    cout << l << endl;

    cout << "print matrix U " << endl;
    MatrixXd u = lu.matrixLU().triangularView<Upper>();
    cout << u << endl;
    //cout << lu.permutationQ() << endl;
    cout << "Original Matrix: " << endl;
    cout << lu.permutationP().inverse()*l*u*lu.permutationQ().inverse() << endl;

    cout << lu.permutationP().inverse() << endl;*/
}

DSparseLU::DSparseLU(MatrixXd smatrix)
{
	this->U = smatrix;
	
	//auto t11 = std::chrono::high_resolution_clock::now();
	MatrixXd ones = MatrixXd::Identity(smatrix.rows(), smatrix.cols());
    //auto t12 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();

    //cout << "Create Identity: " << duration << endl;
	this->L = ones;
	this->P = ones;
	this->Q = ones;

	auto duration = 0;


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

		auto t11 = std::chrono::high_resolution_clock::now();
		MatrixXd tempU = this->U.row(k);
		this->U.row(k) = this->U.row(fila);
		this->U.row(fila) = tempU;

		MatrixXd tempP = this->P.row(k);
		this->P.row(k) = this->P.row(fila);
		this->P.row(fila) = tempP;


		for(int j=0; j<this->U.cols(); j++)
		{
			/*double tempU = this->U(k,j);
			this->U(k,j) = this->U(fila, j);
			this->U(fila, j) = tempU;

			double tempP = this->P(k,j);
			this->P(k,j) = this->P(fila, j);
			this->P(fila, j) = tempP;*/

			if(k > j)
			{
				double tempL = this->L(k,j);
				this->L(k,j) = this->L(fila, j);
				this->L(fila, j) = tempL;
			}
		}
    	auto t12 = std::chrono::high_resolution_clock::now();
    	duration += std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();

		//cout << smatrix << endl;
		for(int i=k+1; i<this->U.rows(); i++)
		{
			float first = this->U(i,k);
			if(this->U(i,k) != 0)
			{
				this->U.row(i) = this->U.row(i)-(first/this->U(k,k))*this->U.row(k);
				/*for(int j=0; j<this->U.cols(); j++)
				{
					this->U(i,j) = this->U(i,j)-(first/this->U(k,k))*this->U(k,j);
					
				}*/
			}
			this->L(i,k) = (first/this->U(k,k));
		}
	}

	this->sL = L.sparseView();
	this->sU = U.sparseView();	
    //cout << "Time Permutations: " << duration << endl;

}

/*DSparseLU::DSparseLU(SparseMatrix<double> A)
{
	LU.analyzePattern(A);
	LU.factorize(A);

}

VectorXd DSparseLU::solve(VectorXd LHS)
{
	VectorXd b = LU.rowsPermutation()*LHS;
	LU.matrixL().solveInPlace(b);
	LU.matrixU().solveInPlace(b);

	return b;
}*/

VectorXd DSparseLU::solve(VectorXd LHS)
{
	/*VectorXd bL = (this->P.transpose()*LHS);

	auto t11 = std::chrono::high_resolution_clock::now();
	SparseLU<SparseMatrix<double>> solverL;
	solverL.analyzePattern(this->sL);
	solverL.factorize(this->sL);
    auto t12 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    //cout << "L to solve: " << duration << endl;
	VectorXd xL = solverL.solve(bL);

	SparseLU<SparseMatrix<double>> solverU;
	solverU.analyzePattern(this->sU);
	solverU.factorize(this->sU);
	VectorXd xU = solverU.solve(xL);

	return xU;*/
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
				L.row(i) = L.row(i)-(first/L(k,k))*L.row(k);
				/*for(int j=0; j<i; j++)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					L(i,j) = L(i,j)-(first/L(k,k))*L(k,j);
				}*/
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
			{
				U.row(i) = U.row(i)-(first/U(k,k))*U.row(k);
				/*
				for(int j=U.cols()-1; j>i; j--)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					U(i,j) = U(i,j)-(first/U(k,k))*U(k,j);
				}*/
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