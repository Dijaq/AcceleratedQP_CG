#ifndef CUDA_FUNCTIONS_H
#define CUDA_FUNCTIONS_H

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

__global__ void cuda_solve_Lx(float *L, float *B, int filas, int columnas, int selected_row, int selected_col)
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
			//int indexB = (columnas*(fil-1))+col;
			int indexB = fil;
			//int index_selected_row_B = (columnas*fil)+selected_col;
			//int selected_first = (columnas*col)+(fil-1);
			//B[indexB] = B[indexB]-B[0]*L[index_selected_col]/L[index_kk];
			B[indexB] =B [indexB]-B[col]*L[index_selected_col]/L[index_kk];
			//B[indexB] =-2;
			//B[col] = 1;
		}

		L[index] = L[index]-L[index_selected_row]*L[index_selected_col]/L[index_kk];

		//L[index] = U[index]+L[index];		
	}

}

__global__ void cuda_solve_Ux(float *U, float *B, int filas, int columnas, int selected_row, int selected_col)
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

}

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

#endif