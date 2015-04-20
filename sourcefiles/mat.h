class mat
{
public:
	int N; // columna
	int M; // fila

	friend void solveLeastSquares(double**matrixQ, double **matrixR, double **matrixB, double * matrixSolution, int M, int N);

	double** matrix; // Main Matrix varirable

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
	void print_thisMatrix(double** thisMatrix, int M, int N); //Prints selected matrix
	double mat_det(void);

	double** mat::product(double** A, double** B, int AM, int AN, int BM, int BN);

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

	void swapRow(double**, int row01, int row02); //Swaps the selected rows

	void transpuesta(void); //Transposes de current matrix and saves it.


	void LUdecomposition(void); // LUDecomposition, seperates them in matrixL and matrixU
	double** matrixL;
	double** matrixU;

	void mat_inv(void);
	double** inv;

	void qr(void);
	double** Q;
	double** R;

	double ** transpose_of_col(double** col, int L);

	void setVariablesToNull(void)
	{
		matrixL = NULL;
		matrixU = NULL;
		matrix = NULL;
		inv = NULL;
	}
	void match_mats(double **, int, int, double **, int, int);

	double** get_col(double **,int col,int M,int N);
	double ** new_mat(int M, int N);


	void cleanMatrix(double **matrix, int M, int N);

	void QRDecomposition(void);
private:


	// Internal functions that are needed to transpose a matrix
	void transpuestaCopyRowToColumn(double** newMatrix, double * rowToCopy, int rowNumber);
	void transposeSetThisNewMatrix(double** newMatrix);
	
	double norm2ofvector(double**,int L);
	//PASS AS COLUMN
	
	void createLUDecompMatrix(void);

	//Creates a matrix MxN, it will return the pointer to the matrix ** 
	double ** newAuxMatrixSamesize(void);
	double ** newAuxIdentityMatrix(void);


	void createIdentityMatrix(void);
	

};

void solveLeastSquares(double **matrixQ, double **matrixR, double **matrixB, double ** matrixSolution, int M, int N);

double** createTimeVectorMatrix(double ** matrixA, double ** timeVector, int totalRows);

double ** standarizationOfTimeValues(double ** tMatrix, int totalRows);

double ** applyLogarithmToDataValues(double ** sMatrix, int totalRows);