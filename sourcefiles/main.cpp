#include<stdio.h>
#include <iostream>
#include <string>
#include"csv.h"
#include"mat.h"
#include "getFileInfo.h"

using namespace std;

int main(int argc, char** argv)
{

	if (argc<4)
	{
		cout << "############################################################" << endl << endl;
		cout << "[DEBUG] This program needs 3 arguments to run:" << endl;
		cout << "[DEBUG] nameA nameB nameX" << endl << endl << endl;
		cout << "[DEBUG] (nameA,nameB and nameX must be the .csv file name)." << endl << endl;
		cout << "############################################################" << endl << endl;
	}
	else
	{
		char* fileNameA = argv[1];
		char* fileNameB = argv[2];
		char* fileNameX = argv[3];

		cout << "[DEBUG] nameA: " << fileNameA << endl;
		cout << "[DEBUG] nameB: " << fileNameB << endl;
		cout << "[DEBUG] nameX: " << fileNameX << endl<<endl;

		/*mat timeRatio = mat(17610, 2); // we create the timeRatio, knowing the amount of lines beforehand
		fileSampleRead("spdc2693.txt", timeRatio.matrix, -1, -1); // we obtain all the data from the file

		mat matrixA = mat(timeRatio.M, 3); //Since we think our polinomialis Ao + A1x + A2x^2 = Y we have a vector 1x3

		long double ** tMatrix = NULL, ** sMatrix = NULL;

		tMatrix = timeRatio.get_col(timeRatio.matrix, 0, timeRatio.M, timeRatio.N); // We seperate from the data the time column
		sMatrix = timeRatio.get_col(timeRatio.matrix, 1, timeRatio.M, timeRatio.N);// We seperate from the data the data values column

		sMatrix = applyLogarithmToDataValues(sMatrix, timeRatio.M); // Since the original values correspond to S = exp{Ao + A1x + A2x^2 = Y}, we have to apply log to each value.
		tMatrix = standarizationOfTimeValues(tMatrix, timeRatio.M); // we standarized the time values so we can work with them without any unit

		matrixA.matrix = createTimeVectorMatrix(matrixA.matrix, tMatrix, timeRatio.M); // we create the time value A matrix ( 1  t  t^2 )

		delete tMatrix; // Since we used this matrix, we don't need it anymore.

		//matrixA.print_mat();
		matrixA.qr();
		//matrixA.QRDecomposition(); // we obtain the QR function to solveLeastSquares


		long double ** matrixX = new long double *[1]; // we create the solution matrix where all the values will go.
		matrixX[0] = new long double[matrixA.N];

		solveLeastSquares(matrixA.Q, matrixA.R, sMatrix, matrixX, timeRatio.M, matrixA.N);

		matrixA.print_mat_R();
		cout << endl;
		matrixA.print_thisMatrix(matrixX, 1, 3);

		//timeRatio.print_thisMatrix(tMatrix, timeRatio.M, 1);
		//timeRatio.print_thisMatrix(sMatrix, timeRatio.M, 1);

		delete matrixX;
		delete sMatrix;




		//.print_thisMatrix(sMatrix, timeRatio.M, 1);
		//timeRatio.print_mat();
		//timeRatio.qr();
		*/
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//	This is a working example of solving least squares with the QR decomposition algorithm method.
		//  - There is a set of examples that will work with this algorithm, such as Example01/02/03 on the csv files.
		//  - It will use the QRDecomposition algorithm.

		mat A= mat(fileNameA);

		mat b = mat(fileNameB);

		long double ** matrixX = new long double * [1]; // we create the solution matrix where all the values will go.
		matrixX[0] = new long double [A.N];
		printf("A:\n");
		A.print_mat();
		printf("b:\n");
		b.print_mat();
		A.qr();
		//A.QRDecomposition();

		solveLeastSquares(A.Q, A.R, b.matrix, matrixX, A.M, A.N);
		printf("X:\n");
		A.print_thisMatrix(matrixX,1,A.N);
		csvWrite(fileNameX, matrixX, 1, A.N);
		delete matrixX;
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%8/


		/*
		cout << endl << endl << "@@@@@@@@@@@@@@@@@@@" << endl << endl;
		A.LUdecomposition();
		cout << "[DEBUG]: MATRIX LOWER" << endl << endl;
		A.print_mat_L();
		cout << "[DEBUG]: MATRIX UPPER" << endl << endl;
		A.print_mat_U();*/

		/*cout << endl << endl << "@@@@@@@@@@@@@@@@@@@" << endl << endl;

		cout << "[DEBUG]: TRANSPOSE" << endl << endl;
		A.transpuesta();
		A.print_mat();

		cout << endl << endl << "@@@@@@@@@@@@@@@@@@@" << endl << endl;
		cout << "[DEBUG]: IDENTITY" << endl << endl;
		mat B = mat("identityMatrix",5,5);
		B.print_mat();*/

		/*mat B = mat("mat2.csv");
		B.print_mat();
		mat C = mat(A.M, B.N);
		C.product(A, B);
		C.print_mat();*/




	}

}

