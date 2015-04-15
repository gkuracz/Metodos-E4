#include<stdio.h>
#include <iostream>
#include"csv.h"
#include"mat.h"

using namespace std;

int main(void)
{
	mat A=mat("mat1.csv");
	A.print_mat();
	cout << endl << endl << "@@@@@@@@@@@@@@@@@@@" << endl << endl;
	A.transpuesta();
	A.print_mat();

	cout << endl << endl << "@@@@@@@@@@@@@@@@@@@" << endl << endl;

	mat B = mat("identityMatrix",5,5);
	B.print_mat();
	
	/*mat B = mat("mat2.csv");
	B.print_mat();
	mat C = mat(A.M, B.N);
	C.product(A, B);
	C.print_mat();*/

	printf("hello");
}

