
void solveLeastSquares(long double **matrixQ, long double **matrixR, long double **matrixB, long double ** matrixSolution, int M, int N);

long double** createTimeVectorMatrix(long double ** matrixA, long double ** timeVector, int totalRows);

long double ** standarizationOfTimeValues(long double ** tMatrix, int totalRows);

long double ** applyLogarithmToDataValues(long double ** sMatrix, int totalRows);


class mat
{
public:
	int N; // columna
	int M; // fila

	friend void solveLeastSquares(long double**matrixQ,long double **matrixR,long double **matrixB,long double * matrixSolution, int M, int N);

	long double** matrix; // Main Matrix varirable

	mat(int, int);//allocates memory for a m by n matrix
	mat(char*); // Constructor that recieved csv

	//if string inserted is "identityMatrix" it will create the identity matrix, must be squared.
	mat(char * specificMatrix, int M, int N);

	~mat();
	void print_mat(void); //Prints current Matrix
	void print_mat_L(void); //Print LUDecomp, L matrix only if it exists.
	void print_mat_U(void); //Print LUDecomp, U matrix only if it exists.
	void print_mat_Q(void); //Print QRDecomp, Q matrix only if it exists.
	void print_mat_R(void); //Print QRDecomp, R matrix only if it exists.
	void print_thisMatrix(long double** thisMatrix, int M, int N); //Prints selected matrix
	long double mat_det(void);
	

	long double** mat::product(long double** A,long double** B, int AM, int AN, int BM, int BN);

	mat& mat::operator=(mat B)
	{
		if (B.M == M)
		{
			if (B.N == N)
			{
				int i, j;
				for (i = 0; i < M; i++)
				{
					for (j = 0; j < N; j++)
					{
						matrix[i][j] = B.matrix[i][j];
					}
				}
			}
		}

		return *this;
	}

	//Verifies if matrix is squared
	bool isMatrixSquared(void);

	void swapRow(long double**, int row01, int row02); //Swaps the selected rows

	void transpuesta(void); //Transposes de current matrix and saves it.


	void LUdecomposition(void); // LUDecomposition, seperates them in matrixL and matrixU
	long double** matrixL;
	long double** matrixU;

	void mat_inv(void);
	long double** inv;

	void qr(void);
	long double** Q;
	long double** R;

	long double ** transpose_of_col(long double** col, int L);

	void setVariablesToNull(void)
	{
		matrixL = NULL;
		matrixU = NULL;
		matrix = NULL;
		inv = NULL;
	}
	void match_mats(long double **, int, int, long double **, int, int);

	long double** get_col(long double **, int col, int M, int N);
	long double ** new_mat(int M, int N);
	void f_get_col(long double **A, int c, int M, int N, long double **col);
	void f_get_trans(long double **A, int M, int N, long double **col);
	void kill(long double**, int, int);
	void cleanMatrix(long double **matrix, int M, int N);

	void QRDecomposition(void);
private:


	// Internal functions that are needed to transpose a matrix
	void transpuestaCopyRowToColumn(long double** newMatrix, long double * rowToCopy, int rowNumber);
	void transposeSetThisNewMatrix(long double** newMatrix);
	
	long double norm2ofvector(long double**, int L);
	//PASS AS COLUMN
	
	void createLUDecompMatrix(void);

	//Creates a matrix MxN, it will return the pointer to the matrix ** 
	long double ** newAuxMatrixSamesize(void);
	long double ** newAuxIdentityMatrix(void);


	void createIdentityMatrix(void);
	

};

