class mat 
{
	public:
		mat(int, int);
		mat(char*);
		~mat();
		void print_mat(void);
		long double** matrix;
		void mat::product(mat& A, mat& B);
		int N;
		int M;

		void swapRow(int row01, int row02); //Swaps the selected rows

		void transpuesta(void); //Transposes de current matrix and saves it.

private:
	// Internal functions that are needed to transpose a matrix
	void transpuestaCopyRowToColumn(long double** newMatrix, long double * rowToCopy, int rowNumber);
	void transposeSetThisNewMatrix(long double** newMatrix);
};