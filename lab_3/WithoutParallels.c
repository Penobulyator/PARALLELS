#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define N1 4096
#define N2 2048	
#define N3 1024

int rank = 0;
int size = 0;

double** createMatrix(int a, int b)
{
	//creates new matrix of size a*b
	double** out = (double**)malloc(a * sizeof(double*));
	for (int i = 0; i < a; i++)
		out[i] = (double*)calloc(b, sizeof(double));
	return out;
}

void deleteMatrix(double** const matrix, int a)
{
	//deletes matrix of size a*x
	for (int i = 0; i < a; i++)
		free(matrix[i]);
	free(matrix);
}

void printMatrix(double** matrix, int a, int b)
{
	for (int i = 0; i < a; i++)
	{
		for (int j = 0; j < b; j++)
			printf("%2f ", matrix[i][j]);
		printf("\n");
	}
	printf("\n");
}

void mult(double** C, double **A, double **B)
{
	//A is matrix of size N1*N2
	//B is matrix of size N2*N3
	//C is matrix of size N1*N3 (C = A*B)
	for (int i = 0; i < N1; ++i)
	{
		for (int j = 0; j < N3; ++j)
		{
			C[i][j] = 0;
			for (int k = 0; k < N2; ++k)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
}
int main()
{
	double** A = createMatrix(N1, N2);
	double** B = createMatrix(N2, N3);
	double** C = createMatrix(N1, N3);

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			A[i][j] = i * j + i + j;


	for (int i = 0; i < N2; i++)
		for (int j = 0; j < N3; j++)
			B[i][j] = i * j + i + j + 1;


	//printMatrix(A, N1, N2);
	//printMatrix(B, N2, N3);

	struct timeval tv1, tv2;
	gettimeofday(&tv1, NULL);

	mult(C, A, B);
	gettimeofday(&tv2, NULL);

	double dt_sec = (tv2.tv_sec - tv1.tv_sec);
	double dt_usec = (tv2.tv_usec - tv1.tv_usec);
	double dt = dt_sec + 1e-6*dt_usec;
	printf("time diff: %e \n", dt);

	//printMatrix(C, N1, N3);

	deleteMatrix(A, N1);
	deleteMatrix(B, N2);
	deleteMatrix(C, N1);

	return 0;
}