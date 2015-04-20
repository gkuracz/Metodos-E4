/****************************************************************
 Este archivo contiene rutinas para leer y escribir matrices de 
 números en punto flotante en archivos "cvs" - 
 "comma separated values".
 Estas rutinas deben ser tomadas como simples "ejemplos", y no
 código bug-free, bullet-proof ni bien comentado.
 Es código "sin garantías".
 El código fue probado en una máquina con Linux.
 24/02/2010
*****************************************************************/
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include"csv.h"
#define _CRT_SECURE_NO_WARNINGS
/* Devuelve el tamaño de una matriz.
ARGS:
 filename = nombre del archivo.
 M = número de filas (pasado por indirección).
 N = número de columnas (pasado por indirección)
RET:
 = 0 => No error.
 > 0 => Error.

Nota: en el caso en que el número de columnas de las filas difiera,
coloca -1 en N.
*/
int csvSize(char* filename, int& M, int& N)
{
	FILE* fp;
	char str[CSV_MAX_LINE_SIZE];
	int i, k, l;
	int err = 0;

	fp = fopen(filename,"rt");
	if( fp == NULL )
		return 1;
	M = 0;
	N = 0;

	while( fgets(str,CSV_MAX_LINE_SIZE,fp) != NULL )
	{
		M++;							// una línea más
		l = strlen(str);
		for( i = 0, k = 1; i < l; i++ ) // cuento número de columnas como número de comas+1
			if( str[i] == CSV_SEPARATOR )
				k++;
		if( M == 1 )
			N = k;
		if( N != k )					// número de columnas equivocado?
		{
			N = -1;
			break;
		}
	}

	fclose(fp);

	return (N<0?2:0);
}

int csvRead(char* filename, long double** mat, int M, int N)
{
	FILE* fp;
	char str[CSV_MAX_LINE_SIZE];
	char str2[CSV_MAX_LINE_SIZE];

	int i, j, k, l;
	int err = 0;

	fp = fopen(filename,"rt");
	if( fp == NULL )
		return 1;

	for( i = 0; i < M; i++ )
	{
		if( fgets(str,CSV_MAX_LINE_SIZE,fp) == NULL )
			break;
		j = 0;
		for( k = 0; k < N-1; k++, j++ )
		{
			for( l = 0; str[j] != CSV_SEPARATOR; l++, j++ )
				str2[l] = str[j];				
			str2[l] = 0x00;
			mat[i][k] = (long double) atof(str2);
		}
		mat[i][k] = (long double) atof(&str[j]);
	}

	fclose(fp);

	return 0;
}

int csvWrite(char* filename, long double** mat, int M, int N)
{
	FILE* fp;

	int i, k;
	int err = 0;

	fp = fopen(filename,"wt+");
	if( fp == NULL )
		return 1;

	for( i = 0; i < M; i++ )
	{
		for( k = 0; k < N-1; k++ )
			fprintf(fp,"%f,",(double)mat[i][k]);
		fprintf(fp,"%f\n",(double)mat[i][k]);
	}

	fclose(fp);

	return 0;
}




