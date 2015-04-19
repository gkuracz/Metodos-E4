#include<stdio.h>
#include <iostream>
#include"mat.h"
#include"csv.h"
#include<math.h>

using namespace std;

mat::mat(int M, int N)
{
	setVariablesToNull();

	this->M = M;
	this->N = N;
	matrix = new long double*[M];
	int i, j;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new long double[N];
	}
	for (i = 0; i < M; i++)
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

mat::mat(char * specificMatrix, int M, int N)
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

long double** mat::product(long double** A, long double** B, int AM, int AN, int BM, int BN)
{

	int i, j, k;
	long double a = 0;
	long double **p = new_mat(AM, BM);
	if (AM == BN)
	{
		for (i = 0; i < AM; i++)
		for (j = 0; j < BN; j++)
		{
			for (k = 0; k < BM; k++)
				a = a + A[i][k] * B[k][j];
			p[i][j] = a;
			a = 0;
		}
	}
	return p;
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

long double mat::mat_det(void)
{
	int i;
	long double dl = 1, du = 1;
	if (M != N)
	{
		fprintf(stderr, "Error, the matrix must be square");
		return 0;
	}
	LUdecomposition();
	//print_mat_L();
	//print_mat_U();
	for (i = 0; i < N; i++)
	{
		dl = dl*matrixL[i][i];
		du = du*matrixU[i][i];
	}
	//fprintf(stderr, "DEBUG  det(A)=%f\n", dl*du);
	return dl*du;


}
void mat::mat_inv(void)
{
	int i, j, k;
	long double t;
	j = 0;
	if (mat_det())
	{
		long double **temp = newAuxMatrixSamesize();
		for (i = 0; i < N*M; i++)
			temp[i] = matrix[i];
		inv = newAuxIdentityMatrix();
		for (i = 0; i < N; i++)
		{
			if (i < N - 1)
			{
				while (temp[i][i] == 0)
				{
					if (i < N - 1)
					{
						swapRow(temp, i, i + 1);
						swapRow(temp, i + 1, N-1);
						swapRow(inv, i, i + 1);
						swapRow(inv, i + 1, N-1);
					}
					i = 0;
				}
			}
			for (j = 0, t = temp[i][i]; j < M; j++)
			{

				inv[i][j] = inv[i][j] / t;
				temp[i][j] = temp[i][j] / t;
			}
			if (i < N - 1)
			{
				for (k = i + 1; k < N; k++)
				{
					for (j = 0, t = temp[k][i]; j < M; j++)
					{

						inv[k][j] = inv[k][j] - inv[i][j] * t;
						temp[k][j] = temp[k][j] - temp[i][j] * t;
					}
				}
			}
		}
		for (k = N-1; k > 0; k--){
			for (i = N-1 - k, t = temp[N - 2 - i][k]; i < N - 1; i++)
			{
				for (j = 0; j < N; j++)
					inv[N - i - 2][j] = inv[N - 2 - i][j] - inv[k][j] * temp[N - 2 - i][k];
				temp[N - i - 2][k] = temp[N - 2 - i][k] - temp[k][k] * temp[N - 2 - i][k];

			}
		}

	}
	//print_thisMatrix(inv, N, N);
}

void mat::transpuestaCopyRowToColumn(long double** newMatrix, long double * rowToCopy, int rowNumber)
{
	int i;
	for (i = 0; i < this->N; i++)
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

void mat::swapRow(long double** matrix, int row01, int row02)
{
	long double * auxRow;
	auxRow = matrix[row01];
	matrix[row01] = matrix[row02];
	matrix[row02] = auxRow;
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

			int i, j, k;
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

void mat::qr(void)
{
	int j, i,k;
	Q = new_mat(M, N);
	long double** aux = newAuxMatrixSamesize();
	R = new_mat(N, N); 
	match_mats(matrix,M,N,Q,M,N);
	for (i = 0; i < N-1; i++)
	{
		R[i][i] = norm2ofvector(get_col(Q, i,M,N),M);

		for (j = 0; j < M; j++)Q[j][i] = Q[j][i] / R[i][i];
		match_mats(Q, M, N, aux, M, N);
		for (k = i+1; k < N; k++)
		{

			for (j = 0; j < M; j++)
			{

				Q[j][k] = aux[j][k] - aux[j][i] * product(transpose_of_col(get_col(aux, k, M, N), M), get_col(aux, i, M, N), 1, N, M, 1)[0][0];

			}

			R[i][k] = product(transpose_of_col(get_col(Q, i, M, N), M), get_col(matrix, k, M, N), 1, N, M, 1)[0][0];
		}
	}

	R[N-1][N-1]=norm2ofvector(get_col(Q, N-1, M, N), M);
	for (i = 0; i < M; i++)
		Q[i][N-1] = Q[i][N-1] / R[N-1][N-1];
	printf("Q:\n");
	print_thisMatrix(Q, M, N);
	printf("R:\n");
	print_thisMatrix(R,N,N);
	
}
long double mat::norm2ofvector(long double** v,int L)
//PASS AS COLUMN
{
	long double n = 0;
	for (int i = 0; i < L; i++)
	{
		n = n + v[i][0] * v[i][0];
	}
	return sqrt(n);
}

void mat::match_mats(long double ** A, int AM, int AN, long double ** B, int BM, int BN)
{
	if ((AN == BN) && (AM == BM))
		for (int i = 0; i < AM; i++)for (int j = 0; j < N;j++)
		B[i][j] = A[i][j];
}

long double** mat::get_col(long double **A,int c,int M, int N)
{
	long double** col = new_mat(M, 1);
	for (int i=0; i < M; i++)
		col[i][0] = A[i][c];
	return col;
}
long double ** mat::new_mat(int M, int N)
{
	long double ** matrix = new long double*[M];
	int i, j;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new long double[N];
	}
	for (i = 0; i < M; i++)
	for (j = 0; j < N; j++)
		matrix[i][j] = 0;
	return matrix;
}
long double ** mat::transpose_of_col(long double** col, int L)
{
	long double ** row = new_mat(1, L);
	for (int i=0; i < L; i++)
		row[0][i] = col[i][0];
	return row;

}

