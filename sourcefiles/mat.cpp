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
	matrix = new double*[M];
	int i, j;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new double[N];
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
	matrix = new double*[M];
	int i;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new double[N];
	}

	csvRead(name, matrix, M, N);

}

mat::mat(char * specificMatrix, int M, int N)
{
	this->M = M;
	this->N = N;

	this->matrix = new double*[M];
	int i;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new double[N];
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

void mat::cleanMatrix(double **matrix, int M, int N)
{
	int i,j;

	for (i = 0; i < M; i++)
	for (j = 0; j < N; j++)
		matrix[i][j] = 0;
}

double ** mat::newAuxMatrixSamesize(void)
{
	double ** auxMatrix;

	this->M = M;
	this->N = N;
	auxMatrix = new double*[this->M];
	int i, j;
	for (i = 0; i < M; i++)
	{
		auxMatrix[i] = new double[this->N];
	}
	for (i = 0; i < this->M; i++)
	for (j = 0; j < this->N; j++)
		auxMatrix[i][j] = 0;

	return auxMatrix;
}

double ** mat::newAuxIdentityMatrix(void)
{
	if (this->M == this->N)
	{
		double ** auxMatrix = newAuxMatrixSamesize();

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

void mat::print_mat_Q(void)
{
	if (this->Q != NULL)
		print_thisMatrix(this->Q, this->M, this->N);
}

void mat::print_mat_R(void)
{
	if (this->R != NULL)
		print_thisMatrix(this->R, this->M, this->N);
}

void mat::print_thisMatrix(double** thisMatrix, int M, int N)
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

double** mat::product(double** A, double** B, int AM, int AN, int BM, int BN)
{

	int i, j, k;
	double a = 0;
	double **p = new_mat(AM, BM);
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
	double ** newMatrix;

	newMatrix = new double*[this->N];

	for (int i = 0; i < this->N; i++)
	{
		newMatrix[i] = new double[this->M];
	}

	// We select each row and copy it on the newMatrix column.
	for (int indexOldRow = 0; indexOldRow < this->M; indexOldRow++)
	{
		double * auxRow;
		auxRow = this->matrix[indexOldRow]; //We copy the entire row to an auxiliar.

		transpuestaCopyRowToColumn(newMatrix, auxRow, indexOldRow); //We transpose the row to a column on the new matrix.
	}

	transposeSetThisNewMatrix(newMatrix);
}

double mat::mat_det(void)
{
	int i;
	double dl = 1, du = 1;
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
	double t;
	j = 0;
	if (mat_det())
	{
		double **temp = newAuxMatrixSamesize();
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

void mat::transpuestaCopyRowToColumn(double** newMatrix, double * rowToCopy, int rowNumber)
{
	int i;
	for (i = 0; i < this->N; i++)
	{
		// Copy each element of the oldRow to the new matrix fixed column rows.
		newMatrix[i][rowNumber] = rowToCopy[i];
	}
}

void mat::transposeSetThisNewMatrix(double** newMatrix)
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

void mat::swapRow(double** matrix, int row01, int row02)
{
	double * auxRow;
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
	double ** T = new_mat(1, M);
	double ** C = new_mat(M, N);
	double** aux = new_mat(M, N);
	double** P = new_mat(1,1);// AuxMatrixSamesize();
	R = new_mat(N, N); 
	match_mats(matrix,M,N,Q,M,N);

	for (i = 0; i < N-1; i++)
	{
		R[i][i] = norm2ofvector(get_col(Q, i,M,N),M);

		for (j = 0; j < M; j++)
		{
			Q[j][i] = (Q[j][i] / R[i][i]);
		}

		match_mats(Q, M, N, aux, M, N);

		for (k = i+1; k < N; k++)
		{
			f_get_col(aux, k, M, N, C);
			f_get_trans(C, M, N, T);
			f_get_col(aux, i, M, N, C);
			P = product(T, C, 1, N, M, 1);
			for (j = 0; j < M; j++)
			{				
				Q[j][k] = aux[j][k] - aux[j][i] * P[0][0];
			}
			kill(P, 1, 1);
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


void mat::f_get_col(double **A, int c, int M, int N, double **col)
{
	for (int i = 0; i < M; i++)
		col[i][0]=A[i][c];

}

void mat::f_get_trans(double **A, int M, int N, double **col)
{
	for (int i = 0; i < M; i++)
		col[0][i] = A[i][0];

}
void mat::kill(double** A, int M, int N)
{
	for (int i = 0; i < M; ++i)
		delete[] A[i];
	delete[] A;
}

double mat::norm2ofvector(double** v,int L)
//PASS AS COLUMN
{
	double n = 0;
	for (int i = 0; i < L; i++)
	{
		n = n + v[i][0] * v[i][0];
	}
	return sqrt(n);
}

void mat::match_mats(double ** A, int AM, int AN, double ** B, int BM, int BN)
{
	if ((AN == BN) && (AM == BM))
	{
		for (int i = 0; i < AM; i++)
		{
			for (int j = 0; j < AN; j++)
			{
				B[i][j] = A[i][j];
			}
		}
	}	
}

double** mat::get_col(double **A,int c,int M, int N)
{
	double** col = new_mat(M, 1);
	for (int i=0; i < M; i++)
		col[i][0] = A[i][c];
	return col;
}
double ** mat::new_mat(int M, int N)
{
	double ** matrix = new double*[M];
	int i, j;
	for (i = 0; i < M; i++)
	{
		matrix[i] = new double[N];
	}
	for (i = 0; i < M; i++)
	for (j = 0; j < N; j++)
		matrix[i][j] = 0;
	return matrix;
}
double ** mat::transpose_of_col(double** col, int L)
{
	double ** row = new_mat(1, L);
	for (int i=0; i < L; i++)
		row[0][i] = col[i][0];

	delete col;
	return row;

}

void solveLeastSquares(double **matrixQ, double **matrixR, double **matrixB, double ** matrixSolution, int M, int N)
{
	double *Y = new double [N];
	int i, k;

	//We will be solving the Rx = Y equation

	//First we set all the values to zero, to avoid floating and random number errors.
	for (k = 0; k<N; k++) 
	{
		Y[k] = 0;
		matrixSolution[0][k] = 0;
	}

	for (k = 0; k < N; k++)
	{
		for (i = 0; i < M; i++)
		{
			// By Multiplying  Q*B we will get (Q^T)*b
			Y[k] += matrixQ[i][k] * matrixB[i][0];
		}	
	}

	for (k = 1; k <= N; k++) //solve Rx=Y using back substitution
	{
		matrixSolution[0][N - k] += Y[N - k];
		for (i = N - 1; i >= 0; i--)
		{
			if (i > (N - k))
			{
				matrixSolution[0][N - k] -= matrixSolution[0][i] * matrixR[N - k][i];
			}		
		}
		matrixSolution[0][N - k] /= matrixR[N - k][N - k];
	}

	delete Y;
}

double** createTimeVectorMatrix(double ** matrixA, double ** timeVector, int totalRows)
{
	// We will be creating the ( 1 t t^2 ) matrix
	// in which it will be the A matrix in the Ax = Y equation
	int i;

	for (i = 0; i < totalRows; i++)
	{
		matrixA[i][0] = 1;
		matrixA[i][1] = timeVector[i][0];
		matrixA[i][2] = timeVector[i][0] * timeVector[i][0];
	}

	return matrixA;
}

double ** standarizationOfTimeValues(double ** tMatrix, int totalRows)
{
	double meanNumber = 0;
	double totalSum = 0;
	int i;

	for (i = 0; i < totalRows; i++)
		totalSum += tMatrix[i][0];

	meanNumber = (totalSum / totalRows);

	//We obtain the Standard Deviation
	for (i = 0, totalSum = 0; i < totalRows; i++)
	{
		totalSum += pow((tMatrix[i][0] - meanNumber), 2);
		
	}

	totalSum = (totalSum / totalRows);
	double standardDeviation = sqrt(totalSum);

	for (i = 0; i < totalRows; i++)
	{
		tMatrix[i][0] -= meanNumber;
		tMatrix[i][0] = (tMatrix[i][0] / standardDeviation);
	}

	return tMatrix;
}

double ** applyLogarithmToDataValues(double ** sMatrix, int totalRows)
{
	int i;

	for (i = 0; i < totalRows; i++)
	{
		sMatrix[i][0] = log(sMatrix[i][0]);
	}

	return sMatrix;
}


void mat::QRDecomposition(void) //using modified Gram-Schmidt
{
	int i, j, k;
	long double norm;

	//We create the QR matrixes where we will be saving all there data
	Q = new_mat(M, N);
	R = new_mat(N, N);

	cleanMatrix(R, N, N); //We clean R and set it to Zero.

	for (k = 0; k<N; k++) // for k= 1->n
	{
		//REPLACE THIS WITH THE NORMALIZE FUNCTION
		for (i = 0, norm = 0; i<M; i++) // R(k,k)= norm 2 of A(1:m,k)
			norm += pow(matrix[i][k], 2);
		R[k][k] = sqrt(norm);

		for (i = 0; i<M; i++)
			Q[i][k] = matrix[i][k] / R[k][k];
		//Q(1:m,k)=A(1:m,k)/R(k,k)

		for (j = k + 1; j<N; j++) //for j=k+1->n
		{
			for (i = 0; i<M; i++)
				R[k][j] += (Q[i][k] * matrix[i][j]);
			//R(k,j)=Q(1:m,k)R(k,j)
			for (i = 0; i<M; i++)
				matrix[i][j] -= (Q[i][k] * R[k][j]);
			//A(1:m,j)=A(1:m,j)-Q(1:m,k)R(k,j)
		}
	}
}
