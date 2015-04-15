#include<stdio.h>
#include <iostream>
#include"mat.h"
#include"csv.h"

using namespace std;

mat::mat(int M, int N)
{
	this->M = M;
	this->N = N;
	matrix = new long double*[M];
	int i,j;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new long double[N];
	}
	for (i = 0; i < M;i++)
	for (j = 0; j < N; j++)
		matrix[i][j] = 0;
}
mat::mat(char*name)
{

		int lM, lN;
		csvSize(name, lM, lN);
		this->M = lM;
		this->N = lN;
		matrix = new long double*[M];
		int i;
		for (i = 0; i < M; i++)
		{
			matrix[i] = new long double[N];
		}

		csvRead(name, matrix, M, N);	
	
}

mat::~mat()
{
	fprintf(stderr, "%d%d", M, N);
	fprintf(stderr, "\n", M, N);;
	;
	for (int i = 0; i < M; ++i)
		delete[] matrix[i];
	delete[] matrix;
}

bool mat::isMatrixSquared(void)
{
	if (this->N == this->M)
	{
		return true;
	}

	return false;
}

void mat::print_mat(void)
{

	int i, j;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(stderr, "%f\t", matrix[i][j]);
		}
		fprintf(stderr, "\n");
	}
}

void mat::product(mat& A, mat& B)
{
	
	int i, j, k;
	long double a = 0;
	if (A.M == B.N)
	{
	for (i = 0; i < A.M;i++)
	for (j = 0; j < B.N; j++)
	{
	for (k = 0; k < B.M; k++)
	a = a + A.matrix[i][k] * B.matrix[k][j];
	this->matrix[i][j] = a;
	a = 0;
	}
	}
}


void mat::transpuesta(void)
{
	// Assuming the matrix is already built, and has M and N defined.
	// We create an auxiliar Matrix where we will transpose it.
	long double ** newMatrix;

	newMatrix = new long double*[this->N];

	for (int i = 0; i < this->N; i++)
	{
		newMatrix[i] = new long double[this->M];
	}

	// We select each row and copy it on the newMatrix column.
	for (int indexOldRow = 0; indexOldRow < this->M; indexOldRow++)
	{
		long double * auxRow;
		auxRow = this->matrix[indexOldRow]; //We copy the entire row to an auxiliar.

		transpuestaCopyRowToColumn(newMatrix, auxRow, indexOldRow); //We transpose the row to a column on the new matrix.
	}

	transposeSetThisNewMatrix(newMatrix);
}

void mat::transpuestaCopyRowToColumn(long double** newMatrix, long double * rowToCopy, int rowNumber)
{
	int i;
	for (i = 0; i < this->N;i++)
	{
		// Copy each element of the oldRow to the new matrix fixed column rows.
		newMatrix[i][rowNumber] = rowToCopy[i];
	}
}

void mat::transposeSetThisNewMatrix(long double** newMatrix)
{

	//First we delete the old matrix, including its elements.
	for (int i = 0; i < M; ++i)
		delete[] matrix[i];
	delete[] matrix;

	//Now we swap the columns/rows
	int auxColumn = this->N;
	this->N = this->M;
	this->M = auxColumn;

	// We replace the old matrix pointer with the new one
	this->matrix = newMatrix;
}

void mat::swapRow(int row01, int row02)
{
	long double * auxRow;
	auxRow = this->matrix[row01];
	this->matrix[row01] = this->matrix[row02];
	this->matrix[row02] = auxRow;
}

mat::mat(char * specificMatrix,int M,int N)
{
	this->M = M;
	this->N = N;

	this->matrix = new long double*[M];
	int i;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new long double[N];
	}

	for (i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			matrix[i][j] = 0;

	if (specificMatrix == "identityMatrix")
	{
		if (M == N)
		{
			for (i = 0; i < this->N; i++)
			{
				this->matrix[i][i] = 1;
			}

			cout << "identity matrix created" << endl;
		}
	}

}