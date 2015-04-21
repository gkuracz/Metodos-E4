#define FILE_DATATITLE 	"Standard and Poor's 500 Index closing values from 1926 to 1993. First column contains the date(yymmdd), second column contains the value. These data are used in : E.Ley(1996) : \"On the Peculiar Distribution of the U.S.Stock Indices; \" forthcoming in The American Statistician.---  "
#define FILE_DATASEPERATOR ' '
#define FILE_DATATITLESIZE 280

enum typeOfData {DATE=0, VALUE};

int fileSampleSize(char* filename, int& M, int& N);
int fileSampleRead(char* filename, long double** mat, int M, int N);
int fileSampleWrite(char* 
