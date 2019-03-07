#include <iostream>
#include <time.h>
#include <string.h>
#include <chrono>

using namespace std;

#define Threads 32

__global__ void add(float *U, float *L, int filas, int columnas, int selected_row, int selected_col)
{
	int col = blockIdx.x*blockDim.x+threadIdx.x;
	int fil = blockIdx.y*blockDim.y+threadIdx.y;

	if((col < columnas && fil < filas) && (fil > selected_row && col > (selected_col-1)))
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

	for (int i=0; i<N; i++) {
		L[i] = 0.0f;
		for (int j=0; j<N; j++) 
			if (i == j) L[i * N + j] = 1.0f;
	}

	//srand(time(NULL));

	float *a = (float *)malloc(N * N * sizeof(float));;

	float *dev_a;
	float *dev_c;

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

	cudaMalloc((void**) &dev_a, filas*columnas*sizeof(float));
	cudaMalloc((void**) &dev_c, filas*columnas*sizeof(float));

	cudaMemcpy(dev_a, a, filas*columnas*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_c, L, filas*columnas*sizeof(float), cudaMemcpyHostToDevice);

	dim3 dimThreadsBloque(Threads, Threads);

	float BFloat = (float) columnas / (float) Threads;
	int B = (int) ceil(BFloat);

	// El grid tendrá B número de bloques en x y y
	dim3 dimBloques(B, B);



	auto t11 = std::chrono::high_resolution_clock::now();
    
    for(int selected=0; selected<3; selected++) 
    {
		add<<<dimBloques, dimThreadsBloque>>>(dev_a, dev_c, filas, columnas, selected, selected);
    }
	
	auto t12 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count();
	cout << "Time to gauss elimination: " << duration << endl; 


	cudaMemcpy(a, dev_a, filas*columnas*sizeof(float), cudaMemcpyDeviceToHost);

	//auto t12 = std::chrono::high_resolution_clock::now();

	cudaFree(dev_a);
	cudaFree(dev_c);

	for(int i=0; i<filas; i++)
	{
		for(int j=0; j<columnas; j++)
		{
			cout << a[i*N+j] << " - ";
		}
		cout << endl;
	}

	//cout << std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count() << endl;

  return 0;
}