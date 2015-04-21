#define CSV_MAX_LINE_SIZE 	1024
#define CSV_SEPARATOR		','

int csvSize(char* filename, int& M, int& N);
int csvRead(char* filename, long double** mat, int M, int N);
int csvWrite(char* filename, long double** mat, int M, int N);

