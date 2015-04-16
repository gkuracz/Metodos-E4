class mat 
{
	public:
		int N;
		int M;

		long double** matrix; // Main Matrix varirable

		mat(int, int);//allocates memory for a m by n matrix
		mat(char*); // Constructor that recieved csv

		//if string inserted is "identityMatrix" it will create the identity matrix, must be squared.
		mat(char * specificMatrix, int M, int N);
		
		~mat();
		void print_mat(void); //Prints current Matrix
		void print_mat_L(void); //Print LUDecomp, L matrix only if it exists.
		void print_mat_U(void); //Print LUDecomp, U matrix only if it exists.
		void print_thisMatrix(long double** thisMatrix, int M, int N); //Prints selected matrix
		int mat_det(void);
		
		void product(mat& A, mat& B);
		
		

		//Verifies if matrix is squared
		bool isMatrixSquared(void);

		void swapRow(long double**, int row01, int row02); //Swaps the selected rows

		void transpuesta(void); //Transposes de current matrix and saves it.


		void LUdecomposition(void); // LUDecomposition, seperates them in matrixL and matrixU
			long double** matrixL;
			long double** matrixU;

		void mat_inv(void);
			long double** inv;

			void setVariablesToNull(void)
			{
				matrixL = NULL;
				matrixU = NULL;
				matrix  = NULL;
				inv		= NULL;
			}

private:
	

	// Internal functions that are needed to transpose a matrix
	void transpuestaCopyRowToColumn(long double** newMatrix, long double * rowToCopy, int rowNumber);
	void transposeSetThisNewMatrix(long double** newMatrix);

	void createLUDecompMatrix(void);

	//Creates a matrix MxN, it will return the pointer to the matrix ** 
	long double ** newAuxMatrixSamesize(void);
	long double ** newAuxIdentityMatrix(void);


	void createIdentityMatrix(void);
	
};
