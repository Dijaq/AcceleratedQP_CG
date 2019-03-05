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

/*this constructor is not used*/

DSparseLU::DSparseLU(MatrixXd smatrix, bool s)
{
	cout << "constructor" << endl;
	
	this->P = MatrixXd::Identity(smatrix.rows(), smatrix.cols());
	this->Q = MatrixXd::Identity(smatrix.rows(), smatrix.cols());

	this->L = MatrixXd::Identity(smatrix.rows(), smatrix.cols());

	/*for(int i=smatrix.rows()-1; i>=0; i--)
	{
		int contador = 1;
		for(int j=smatrix.cols()-1; j>=0; j--)
		{
			if(smatrix(i,j) != 0)
			{
				permutationP(i,j) = permutationP(i,j)+contador;
				contador++;
			}
		}

	}

	for(int j=smatrix.rows()-1; j>=0; j--)
	{
		int contador = 1;
		for(int i=smatrix.cols()-1; i>=0; i--)
		{
			if(smatrix(i,j) != 0)
			{
				permutationQ(i,j) = permutationQ(i,j)+contador;
				contador++;
			}
		}

	}*/

//Create matriz de permutation of LU
	this->U = smatrix;

	/*cout << "Permutation P" << endl;
	cout << permutationP << endl;

	cout << "Permutation Q" << endl;
	cout << permutationQ << endl;*/

	for(int k=0; k<smatrix.rows(); k++)
	{
		//This matrix contain the number of elements at rigth of each row
		MatrixXd permutationP = MatrixXd::Zero(smatrix.rows(), smatrix.cols());
		//This matrix contain the number of elements at down of each col
		MatrixXd permutationQ = MatrixXd::Zero(smatrix.rows(), smatrix.cols());
		
		/*Create the Matrix of Permutation of count rows and cols*/
		for(int i=this->U.rows()-1; i>=k; i--)
		{
			int contador = 1;
			for(int j=this->U.cols()-1; j>=k; j--)
			{
				if(this->U(i,j) != 0)
				{
					permutationP(i,j) = permutationP(i,j)+contador;
					contador++;
				}
			}

		}
		
		/*End of the Matrix of Permutation of count rows and cols*/

		int select_row = k;
		int value_row = 0;

		for(int i=k; i<permutationP.rows(); i++)
		{
			if(value_row < permutationP(i,k))
			{
				value_row = permutationP(i,k);
				select_row = i;
			}	
		}

		MatrixXd tempU_row;
		tempU_row = this->U.row(k);
		this->U.row(k) = this->U.row(select_row);
		this->U.row(select_row) = tempU_row;

		/*MatrixXd tempPermuP;
		tempPermuP = permutationP.row(k);
		permutationP.row(k) = permutationP.row(select_row);
		permutationP.row(select_row) = tempPermuP;*/

		MatrixXd tempP;
		tempP = this->P.row(k);
		this->P.row(k) = this->P.row(select_row);
		this->P.row(select_row) = tempP;


		/*Create the Matrix of Permutation of count rows and cols*/
		for(int j=this->U.rows()-1; j>=k; j--)
		{
			int contador = 1;
			for(int i=this->U.cols()-1; i>=k; i--)
			{
				if(this->U(i,j) != 0)
				{
					permutationQ(i,j) = permutationQ(i,j)+contador;
					contador++;
				}
			}

		}
		/*End of the Matrix of Permutation of count rows and cols*/

		int select_col = k;
		int value_col = smatrix.cols()+1;

		for(int j=k; j<smatrix.cols(); j++)
		{
			if(permutationQ(k,j) < value_col && permutationQ(k,j) > 0)
			{
				value_col = permutationQ(k,j);
				select_col = j;
			}	
		}		

		MatrixXd tempU_col;
		tempU_col = this->U.col(k);
		this->U.col(k) = this->U.col(select_col);
		this->U.col(select_col) = tempU_col;

		/*MatrixXd tempPermuQ;
		tempPermuQ = permutationQ.col(k);
		permutationQ.col(k) = permutationQ.col(select_col);
		permutationQ.col(select_col) = tempPermuQ;*/

		MatrixXd tempQ;
		tempQ = this->Q.col(k);
		this->Q.col(k) = this->Q.col(select_col);
		this->Q.col(select_col) = tempQ;

		//Permutar values in L//Falta verificar
		for(int j=0; j<this->U.cols(); j++)
		{
			if(k > j)
			{
				double tempL = this->L(k,j);
				this->L(k,j) = this->L(select_row, j);
				this->L(select_row, j) = tempL;
			}
		}

		/*Eliminacion de Gaussiana*/
		for(int i=k+1; i<this->U.rows(); i++)
		{
			if(this->U(i,k) != 0)
			{
				double first_element = this->U(i,k);
				this->U.row(i) = this->U.row(i)-(first_element/this->U(k,k))*this->U.row(k);

				this->L(i,k) = (first_element/this->U(k,k));
			}
		}
	}

	//export_mat_to_excel(this->U, "ValidarDatos/U_prueba_in");


	/*cout << "Matrix U" << endl;
	cout << this->U << endl;

	cout << "Matrix L" << endl;
	cout << this->L << endl;

	cout << "Init Matrix" << endl;
	cout << this->P.transpose()*this->L*this->U*this->Q.transpose();*/

	/*cout << "Matrix P" << endl;
	cout << this->P << endl;

	cout << "Matrix Q" << endl;
	cout << this->Q << endl;

	cout << permutationP << endl;
	cout << permutationQ << endl;

	cout << "imprimir " << endl;*/
}
/*LU implementation of Eigen
//Toma 8 veces mas tiempo que la implementacion propia
DSparseLU::DSparseLU(MatrixXd smatrix, bool s)
{
	pivlu = FullPivLU<MatrixXd>(smatrix);
	cout << "Original Matrix" << endl;
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

    cout << lu.permutationP().inverse() << endl;
}*/

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

	//export_mat_to_excel(this->L, "ValidarDatos/LL_in");
	//export_mat_to_excel(this->U, "ValidarDatos/UU_in");	

	this->sL = this->L.sparseView();
	this->sU = this->U.sparseView();	
    //cout << "Time Permutations: " << duration << endl;

}

DSparseLU::DSparseLU(SparseMatrix<double> A)
{
	//cout << "------:::::::>>>>>" << A.nonZeros();
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
	VectorXd bL = (this->P*LHS);

	//auto t11 = std::chrono::high_resolution_clock::now();
	SparseLU<SparseMatrix<double>> solverL;
	solverL.analyzePattern(this->sL);
	solverL.factorize(this->sL);
    //auto t12 = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
    //cout << "L to solve: " << duration << endl;
	VectorXd xL = solverL.solve(bL);

	SparseLU<SparseMatrix<double>> solverU;
	solverU.analyzePattern(this->sU);
	solverU.factorize(this->sU);
	VectorXd xU = solverU.solve(xL);

	return xU;

	//print_dimensions( "PPP: ", this->P);
	//print_dimensions( "QQQ: ", this->Q);

	//return this->Q.transpose()*(solve_Ux(this->U,solve_Lx(this->L, (this->P.transpose()*LHS))));
	//return (solve_Ux(this->U,solve_Lx(this->L, (this->P*LHS))));
}

VectorXd DSparseLU::solve_Lx(MatrixXd L, VectorXd X)
{
	//export_mat_to_excel(L, "ValidarDatos/LLL_in");
	for(int k=0; k<L.rows(); k++)
	{
		//cout << smatrix << endl;
		for(int i=k+1; i<L.rows(); i++)
		{
			if(L(i,k) != 0)
			{
				float first = L(i,k);
				L.row(i) = L.row(i)-(first/L(k,k))*L.row(k);
				/*for(int j=0; j<i; j++)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					L(i,j) = L(i,j)-(first/L(k,k))*L(k,j);
				}*/
				X(i,0) = X(i,0)-(first/L(k,k))*X(k,0);
			}

			//cout << "Num/Den: " << first << "/" << L(k,k) << " val: " << X(i,0) << endl;
		}
	}

	//export_mat_to_excel(L, "ValidarDatos/LLL_out");
	//cout << "L: " <<L << endl;
	//cout << "X: " << X << endl;
	return X;
}

VectorXd DSparseLU::solve_Ux(MatrixXd U, VectorXd X)
{
	//export_mat_to_excel(U, "ValidarDatos/UUU_in");
	for(int k=U.rows()-1; k>=0; k--)
	{
		//cout << smatrix << endl;
		for(int i=k-1; i>=0; i--)
		{
			if(U(i,k) != 0)
			{
				float first = U(i,k);
				U.row(i) = U.row(i)-(first/U(k,k))*U.row(k);
				/*
				for(int j=U.cols()-1; j>i; j--)
				{
					//cout << "i: " << i << " j: " << j << " values: " << smatrix(i,k) <<" - " <<smatrix(k,k) << " - " <<smatrix(k,j) << endl;
					U(i,j) = U(i,j)-(first/U(k,k))*U(k,j);
				}*/
				X(i,0) = X(i,0)-(first/U(k,k))*X(k,0);
			}

			//cout << "Num/Den: " << first << "/" << L(k,k) << " val: " << X(i,0) << endl;
		}
	}

	for(int i=0; i< U.rows(); i++)
	{
		float div = U(i,i);
		U(i,i) /= div;
		X(i,0) /= div;
	}

	//export_mat_to_excel(U, "ValidarDatos/UUU_out");
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