#include <iostream>
#include <time.h>
#include <string.h>
#include <chrono>

using namespace std;

#define Threads 32

//This method return the matrix L and U of a LU factorization
__global__ void cuda_LU_factorization(float *U, float *L, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil > selected_row && col > (selected_col-1)))
	{
		int index_selected_row = (columnas*selected_row)+col;
		int index_selected_col = (columnas*fil)+selected_col;
		int index_kk = (columnas*selected_row)+selected_col;
		
		int index = (columnas*fil)+col;

		if(col == selected_col)
		{
			L[index] = U[index_selected_col]/U[index_kk];
		}

		U[index] = U[index]-U[index_selected_row]*U[index_selected_col]/U[index_kk];
		//L[index] = U[index]+L[index];		
	}

}

__global__ void cuda_permutationP(float *A, float *P, int filas, int columnas, int selected_row, int selected_col, int *dev_per)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	int index_kk = (columnas*selected_row)+selected_col;

	if((A[index_kk] == 0 && col == selected_col) && (col < columnas && fil < filas))
	{
		if(fil == selected_row)
			*dev_per = fil;	
	}

}

__global__ void cuda_solve_Lx(float *L, float *B, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	for(int selected=0; selected<filas-1; selected++)
    {
    	selected_row = selected;
    	selected_col = selected;
    	if((col < columnas && fil < filas) && (fil > selected_row && col > (selected_col-1)))
		{
			//Index for Matrix L
			int index_selected_row = (columnas*selected_row)+col;
			int index_selected_col = (columnas*fil)+selected_col;
			int index_kk = (columnas*selected_row)+selected_col;
			int index = (columnas*fil)+col;

			//Index for Matrix B
			if(fil > selected_row && col == selected_col)
			{
				
				int indexB = fil;
				B[indexB] =B [indexB]-B[col]*L[index_selected_col]/L[index_kk];
			}

			L[index] = L[index]-L[index_selected_row]*L[index_selected_col]/L[index_kk];
		}
    }

	

}

/*__global__ void cuda_solve_Lx(float *L, float *B, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil > selected_row && col > (selected_col-1)))
	{
		//Index for Matrix L
		int index_selected_row = (columnas*selected_row)+col;
		int index_selected_col = (columnas*fil)+selected_col;
		int index_kk = (columnas*selected_row)+selected_col;
		int index = (columnas*fil)+col;

		//Index for Matrix B
		if(fil > selected_row && col == selected_col)
		{
			
			int indexB = fil;
			B[indexB] =B [indexB]-B[col]*L[index_selected_col]/L[index_kk];
		}

		L[index] = L[index]-L[index_selected_row]*L[index_selected_col]/L[index_kk];
	}

}*/

__global__ void cuda_solve_Ux(float *U, float *B, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	for(int selected=filas-1; selected>=0; selected--) 
    {
    	selected_row = selected;
    	selected_col = selected;

		if((col < columnas && fil < filas) && (fil < selected_row && col < (selected_col+1)))
		{
			int index_selected_row = (columnas*selected_row)+col;
			int index_selected_col = (columnas*fil)+selected_col;
			int index_kk = (columnas*selected_row)+selected_col;
			
			int index = (columnas*fil)+col;

			if(fil<selected_row && col == selected_col)
			{
				int indexB = fil;
				//B[indexB] = B[col];
				//B[indexB] = B[col];
				B[indexB] = B[indexB]-B[col]*U[index_selected_col]/U[index_kk];;
			}

			U[index] = U[index]-U[index_selected_row]*U[index_selected_col]/U[index_kk];
			//L[index] = U[index]+L[index];		
		}
	}

}

/*__global__ void cuda_solve_Ux(float *U, float *B, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil < selected_row && col < (selected_col+1)))
	{
		int index_selected_row = (columnas*selected_row)+col;
		int index_selected_col = (columnas*fil)+selected_col;
		int index_kk = (columnas*selected_row)+selected_col;
		
		int index = (columnas*fil)+col;

		if(fil<selected_row && col == selected_col)
		{
			int indexB = fil;
			//B[indexB] = B[col];
			//B[indexB] = B[col];
			B[indexB] = B[indexB]-B[col]*U[index_selected_col]/U[index_kk];;
		}

		U[index] = U[index]-U[index_selected_row]*U[index_selected_col]/U[index_kk];
		//L[index] = U[index]+L[index];		
	}

}*/

__global__ void cuda_reduce_U(float *U, float *B, int filas, int columnas)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil == col))
	{
		int index = (columnas*fil)+col;
		B[fil] = B[fil]/U[index];

		U[index] = U[index]/U[index];
	}

}


int main()
{
	srand(time(NULL));
	int filas = 4;
	int columnas = 4;
	int N = filas;
	float *L = (float *)malloc(N * N * sizeof(float));
	float *xB = (float *)malloc(1 * N * sizeof(float));
	float *a = (float *)malloc(N * N * sizeof(float));
	float *P = (float *)malloc(N * N * sizeof(float));

	/*xB[0] = -1;
	xB[1] = 3;
	xB[2] = 2;
	xB[3] = 5;*/

	for(int i=0; i<N; i++)
	{
		xB[i] = rand()%5+1;
	}

//Create matrix L
	for (int i=0; i<N; i++) {
		L[i] = 0.0f;
		for (int j=0; j<N; j++) 
			if (i == j)
			{
				L[i * N + j] = 1.0f;
			} 
	}

	//srand(time(NULL));

//Create Matrix A
	float *dev_U;
	float *dev_L;
	float *dev_B;

	for(int i=0; i<filas; i++)
	{
		for(int j=0; j<columnas; j++)
		{
		  a[i*N+j] = rand()%10+1;
		  cout << a[i*N+j] << " - ";
		}
		cout << endl;
	}

	//a[0] = 0;

	/*for(int i=0; i<filas; i++)
	{
		for(int j=0; j<columnas; j++)
		{
			cout << a[i][j] << " - ";
		}
		cout << endl;
	}*/

	cudaMalloc((void**) &dev_U, filas*columnas*sizeof(float));
	cudaMalloc((void**) &dev_L, filas*columnas*sizeof(float));
	cudaMalloc((void**) &dev_B, filas*1*sizeof(float));

	cudaMemcpy(dev_U, a, filas*columnas*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_L, L, filas*columnas*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_B, xB, filas*1*sizeof(float), cudaMemcpyHostToDevice);

	dim3 dimThreadsBloque(Threads, Threads);

	float BFloat = (float) columnas / (float) Threads;
	int B = (int) ceil(BFloat);

	// El grid tendrá B número de bloques en x y y
	dim3 dimBloques(B, B);

	//int *dev_per;
	//cudaMalloc((void**) &dev_per, 1*sizeof(int));

	auto t11 = std::chrono::high_resolution_clock::now();
    
    //LU factorization
    for(int selected=0; selected<filas-1; selected++) 
    {
		cuda_LU_factorization<<<dimBloques, dimThreadsBloque>>>(dev_U, dev_L, filas, columnas, selected, selected);
    }

    /*for(int selected=0; selected<filas-1; selected++)
    {
    }*/
	cuda_solve_Lx<<<dimBloques, dimThreadsBloque>>>(dev_L, dev_B, filas, columnas, 0, 0);

    /*for(int selected=filas-1; selected>=0; selected--) 
    {
    }*/
	cuda_solve_Ux<<<dimBloques, dimThreadsBloque>>>(dev_U, dev_B, filas, columnas, 0, 0);

//Final reduce U
    cuda_reduce_U<<<dimBloques, dimThreadsBloque>>>(dev_U, dev_B, filas, columnas);
	
	auto t12 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
	cout << "Time to gauss elimination: " << duration << endl; 


	cudaMemcpy(a, dev_U, filas*columnas*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(L, dev_L, filas*columnas*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(xB, dev_B, filas*1*sizeof(float), cudaMemcpyDeviceToHost);

	//auto t12 = std::chrono::high_resolution_clock::now();

	cudaFree(dev_U);
	cudaFree(dev_L);

	cout << "print U: " << endl;
	for(int i=0; i<filas; i++)
	{
		for(int j=0; j<columnas; j++)
		{
			cout << a[i*N+j] << " - ";
		}
		cout << endl;
	}

	cout << "print L: " << endl;

	for(int i=0; i<filas; i++)
	{
		for(int j=0; j<columnas; j++)
		{
			cout << L[i*N+j] << " - ";
		}
		cout << endl;
	}

	cout << "print B: " << endl;
	for(int i=0; i<filas; i++)
	{
		cout << xB[i] << endl;
	}

	//cout << std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count() << endl;

  return 0;
}