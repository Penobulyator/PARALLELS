#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h> 

#define N 1024
double** createMatrix()
{
	double** out = (double**)malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++)
		out[i] = (double*)calloc(N , sizeof(double));
	return out;
}
void deleteMatrix(double** const matrix)
{
	for (int i = 0; i < N; i++)
		free(matrix[i]);
	free(matrix);
}
double* createVector()
{
	double* out = (double*)calloc(N , sizeof(double));
	return out;
}
void printMatrix(double** matrix)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			printf("%f ", matrix[i][j]);
		printf("\r\n");
	}
	printf("\r\n");
}
void printVector(double* vector)
{
	for (int i = 0; i < N; i++)
	{
		printf("%f ", vector[i]);
		printf("\r\n");
	}
	printf("\r\n");
}
double metric(double* vector)
{
	double out = 0;
	for (int i = 0; i < N; i++)
		out += vector[i] * vector[i];
	return sqrt(out);
}
void countDiff(double* diff, double *const *const A, double *const b, double *const x)
{
	//count Ax
	for (int i = 0; i < N; i++)
	{
		diff[i] = 0;
		for (int j = 0; j < N; j++)
			diff[i] += A[i][j] * x[j];
	}
	//count Ax - b
	for (int i = 0; i < N; i++)
		diff[i] -= b[i];
}
void refreshRoot(double* root, double *const diff, const double t)
{
	//root = root - t * diff
	for (int i = 0; i < N; i++)
		root[i] = root[i] - t * diff[i];
}
double* root(double** const A, double* const b)
{
	double* cur_root = (double*)calloc(N, sizeof(double));
	double* diff = (double*)calloc(N, sizeof(double));
	const double t = 0.0001;
	const double epsilon = 0.001;
	while (1)
	{
		countDiff(diff, A, b, cur_root);
		if (metric(diff) / metric(b) < epsilon)
		{
			free(diff);
			return cur_root;
		}
		else
		{
			refreshRoot(cur_root, diff, t);
		}
	}
	
}
int main()
{
	double** A = createMatrix();
	for (int i = 0; i < N; i++)
		A[i][i] = 1;

	double* b = createVector();
	for (int i = 0; i < N; i++)
		b[i] = 1; 

	clock_t t = clock();
	double* x = root(A, b); 
	t = clock() - t;

	for (int i = 0; i < N; i++)
		printf("%2f\n", x[i]);

	printf("time diff: %e \n", ((double)t) / CLOCKS_PER_SEC);

	deleteMatrix(A);
	free(b);

	return 0;
}