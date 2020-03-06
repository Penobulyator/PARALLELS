#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define N  1000
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
	for (int i = 0; i < N; i++)
		out[i]= (double)(rand() % 10 + 1) / 10; //double from 0 to 1.0
	return out;
}
double* matrixMULTvectror(double** const matrix, double* const vector)
{
	double* out = (double*)calloc(N , sizeof(double));
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			out[i] += matrix[i][j] * vector[i];
	return out;
}
double* vectorMINUSvector(double* a, double* b)
{
	double* out = (double*)malloc(N * sizeof(double));
	for (int i = 0; i < N; i++)
		out[i] = a[i] - b[i];
	return out;
}
double metric(double* vector)
{
	double out = 0;
	for (int i = 0; i < N; i++)
		out += vector[i] * vector[i];
	return sqrt(out);
}
double* multScal(double* const vector, double num)
{
	double* out = (double*)malloc(N * sizeof(double));
	for (int i = 0; i < N; i++)
		out[i] = vector[i] * num;
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
double* root(double** const A, double* const b)
{
	double* cur_root = (double*)calloc(N, sizeof(double));
	double t = 0.0001;
	const double epsilon = 0.00001;
	while (1)
	{
		double* mult = matrixMULTvectror(A, cur_root);
		double* diff = vectorMINUSvector(mult, b);
		if (metric(diff) / metric(b) < epsilon)
		{
			free(mult);
			free(diff);
			return cur_root;
		}
		else
		{
			double* tdiff = multScal(diff, t);
			double* new_root = vectorMINUSvector(cur_root, tdiff);
			cur_root = new_root;
			free(tdiff);
			free(mult);
			free(diff);
		}
	}
	
}
int main()
{
	double** A = createMatrix();
	for (int i = 0; i < N; i++)
		A[i][i] = i + 1;
	double* b = createVector();
	for (int i = 0; i < N; i++)
		b[i] = i + 1;
	double* x = root(A, b);
	printVector(x);
	deleteMatrix(A);
	free(b);
	return 0;
}