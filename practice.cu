#include <iostream>
#include <time.h>
#include <string.h>
#include <chrono>

using namespace std;

#define Threads 32

//This method return the matrix L and U of a LU factorization
__global__ void LU_factorization(float *U, float *L, int filas, int columnas, int selected_row, int selected_col)
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

__global__ void solve_Lx(float *L, float *B, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil > selected_row && col > (selected_col-1)))
	{
		int index_selected_row = (columnas*selected_row)+col;
		int index_selected_col = (columnas*fil)+selected_col;
		int index_kk = (columnas*selected_row)+selected_col;
		
		int index = (columnas*fil)+col;

		L[index] = L[index]-L[index_selected_row]*L[index_selected_col]/L[index_kk];
		//L[index] = U[index]+L[index];		
	}

}

__global__ void solve_Ux(float *U, float *B, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil < selected_row && col < (selected_col+1)))
	{
		int index_selected_row = (columnas*selected_row)+col;
		int index_selected_col = (columnas*fil)+selected_col;
		int index_kk = (columnas*selected_row)+selected_col;
		
		int index = (columnas*fil)+col;

		U[index] = U[index]-U[index_selected_row]*U[index_selected_col]/U[index_kk];
		//L[index] = U[index]+L[index];		
	}

}

int main()
{
	int filas = 4;
	int columnas = 4;
	int N = filas;
	float *L = (float *)malloc(N * N * sizeof(float));
	float *a = (float *)malloc(N * N * sizeof(float));

	for (int i=0; i<N; i++) {
		L[i] = 0.0f;
		for (int j=0; j<N; j++) 
			if (i == j) L[i * N + j] = 1.0f;
	}

	//srand(time(NULL));

	float *dev_U;
	float *dev_L;

	for(int i=0; i<filas; i++)
	{
		for(int j=0; j<columnas; j++)
		{
		  a[i*N+j] = rand()%10+1;
		  cout << a[i*N+j] << " - ";
		}
		cout << endl;
	}

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

	cudaMemcpy(dev_U, a, filas*columnas*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_L, L, filas*columnas*sizeof(float), cudaMemcpyHostToDevice);

	dim3 dimThreadsBloque(Threads, Threads);

	float BFloat = (float) columnas / (float) Threads;
	int B = (int) ceil(BFloat);

	// El grid tendrá B número de bloques en x y y
	dim3 dimBloques(B, B);



	auto t11 = std::chrono::high_resolution_clock::now();
    
    //LU factorization
    for(int selected=0; selected<filas-1; selected++) 
    {
		LU_factorization<<<dimBloques, dimThreadsBloque>>>(dev_U, dev_L, filas, columnas, selected, selected);
    }

    for(int selected=filas-1; selected>=0; selected--) 
    {
		solve_Ux<<<dimBloques, dimThreadsBloque>>>(dev_U, dev_L, filas, columnas, selected, selected);
    }
	
	auto t12 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
	cout << "Time to gauss elimination: " << duration << endl; 


	cudaMemcpy(a, dev_U, filas*columnas*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(L, dev_L, filas*columnas*sizeof(float), cudaMemcpyDeviceToHost);

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

	//cout << std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count() << endl;

  return 0;
}