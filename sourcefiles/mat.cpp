#include<stdio.h>
#include <iostream>
#include"mat.h"
#include"csv.h"

using namespace std;

mat::mat(int M, int N)
{
	setVariablesToNull();

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
		setVariablesToNull();
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
		}
	}

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

long double ** mat::newAuxMatrixSamesize(void)
{
	long double ** auxMatrix;

	this->M = M;
	this->N = N;
	auxMatrix = new long double*[this->M];
	int i, j;
	for (i = 0; i < M; i++)
	{
		auxMatrix[i] = new long double[this->N];
	}
	for (i = 0; i < this->M; i++)
	for (j = 0; j < this->N; j++)
		auxMatrix[i][j] = 0;

	return auxMatrix;
}

long double ** mat::newAuxIdentityMatrix(void)
{
	if (this->M == this->N)
	{
		long double ** auxMatrix = newAuxMatrixSamesize();

		for (int i = 0; i < this->N; i++)
		{
			auxMatrix[i][i] = 1;
		}

		return auxMatrix;
	}

	return NULL;
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
	print_thisMatrix(this->matrix, this->M, this->N);
}

void mat::print_mat_L(void)
{
	if (this->matrixL != NULL)
		print_thisMatrix(this->matrixL, this->M, this->N);
}

void mat::print_mat_U(void)
{
	if (this->matrixU != NULL)
		print_thisMatrix(this->matrixU, this->M, this->N);
}

void mat::print_thisMatrix(long double** thisMatrix, int M, int N)
{
	int i, j;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(stderr, "%f\t", thisMatrix[i][j]);
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


/*
	Function that will create two matrix L and U necesary for the LUDecomposition algorithm.
*/
void mat::LUdecomposition(void)
{
	if (isMatrixSquared())
	{
		if ((this->matrixL == NULL) && (this->matrixU == NULL))
		{
			createLUDecompMatrix(); // Craetes too necesary matrix for LUDecomposition, same size as the original.

			int i, j,k;
			double sum;

			for (j = 0; j < this->N; j++) 
			{

				// Solving upper matrix
				for (i = 0; i <= j; i++) 
				{
					sum = 0.0;
					for (k = 0; k < i; k++) 
					{
						sum += matrixU[k][j] * matrixL[i][k];
					}

					matrixU[i][j] = matrix[i][j] - sum;
				}

				// Solving lower Matrix
				for (i = j + 1; i < this->N; i++)
				{
					sum = 0.0;
					for (k = 0; k < j; k++)
					{
						sum += matrixU[k][j] * matrixL[i][k];
					}

					matrixL[i][j] = (1.0 / matrixU[j][j]) * (matrix[i][j] - sum);
				}
			}
		}
	}
}

void mat::createLUDecompMatrix(void)
{
	this->matrixL = newAuxIdentityMatrix(); //Create lower Matrix with identity diagonal
	this->matrixU = newAuxMatrixSamesize(); //Create upperMatrix
}
