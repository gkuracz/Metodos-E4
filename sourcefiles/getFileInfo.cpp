#include<stdio.h>
#include<string>
#include<stdlib.h>
#include <iostream>


#include"getFileInfo.h"
#include"mat.h"

#define _CRT_SECURE_NO_WARNINGS

using namespace std;

int fileSampleRead(char* filename, long double** mat, int M, int N)
{
	FILE* fp;
	
	fp = fopen(filename, "rt");
	if (fp == NULL)
		return 1;

	char c;
	string value;
	string date;

	bool typeOfData = DATE; //If false its value, if true its date
	bool wasDataComplete = true;

	int indexData = 0;
	int charsToSkip = 0;

	do 
	{
		c = fgetc(fp);

		if (charsToSkip < FILE_DATATITLESIZE)
		{
			charsToSkip++;
			cout << c;
		}
		else
		{

			if ((c != FILE_DATASEPERATOR) && (c != '\n')) //different from space, save data since its a value
			{
				if (typeOfData == DATE)
				{
					//It is the date of the current data
					date += c;
				}
				else if (typeOfData == VALUE)
				{
					value += c;
				}

				wasDataComplete = false;
			}
			else
			{
				if (wasDataComplete == false)
				{
					if (typeOfData == DATE)
					{
						//cout << date << "  ";
						mat[indexData][0] = atof(date.c_str()); //Saves  the current date
						date.clear(); //Cleans value stored
						typeOfData = VALUE; // Next information to be stored is the Value of this selected date

					}
					else
					if (typeOfData == VALUE)
					{
						//cout << value << endl;
						mat[indexData][1] = atof(value.c_str()); //Saves  the current date's data
						value.clear(); //Cleans value stored

						typeOfData = DATE; // Next information to be stored is the next time/data starting from next row.
						indexData++; //Next information recieved will be a new data with its DATE/VALUE
					}

					wasDataComplete = true; //The matrix slot was complete, so it will be set to true.
				}
			}
		}

	} while (c != EOF);

	return 0;
}

long double* createTimeVectorValues(long double * row, long double * auxiliarRow)
{
	//We will be creating the ( 1 t t^2 ) matrix, 1 row 3 columns
	auxiliarRow[0] = 1;
	auxiliarRow[1] = row[0];
	auxiliarRow[2] = row[0] * row[0];

	return auxiliarRow;
}