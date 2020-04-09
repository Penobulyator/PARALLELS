#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int size = 0;    //count of processes
int rank = 0;    //number of process
#define N  100   //matrix size
#define M N/size //number of rows for one process
double** createMatrix()
{
	double** out = (double**)malloc(M * sizeof(double*));
	for (int i = 0; i < M; i++)
		out[i] = (double*)calloc(N, sizeof(double));
	return out;
}
void deleteMatrix(double** matrix)
{
	for (int i = 0; i < M; i++)
		free(matrix[i]);
	free(matrix);
}
double metric(double *const vector)
{
	double out = 0;
	for (int i = 0; i < M; i++)
		out += vector[i] * vector[i];
	return sqrt(out);
}
void printMatrix(double *const *const A)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
			printf("%lf ", A[i][j]);
		printf("\n");
	}
	printf("\n");
}
void printVector(double *const vector)
{
	for (int i = 0; i < M; i++)
	{
		printf("%lf ", vector[i]);
		printf("\n");
	}
	printf("\n");
}
void getFullRoot(double* full_root, double *const x)
{
	//send part of root to other processes
	for (int i = 0; i < size; i++)
		if (i != rank)
			MPI_Send(x, M, MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
	//receive missing parts of root from other processes and build a full root
	double* buf = (double*)calloc(M, sizeof(double));
	for (int i = 0; i < size; i++)
	{
		if (i == rank)
		{
			for (int j = 0; j < M; j++)
				full_root[i*M + j] = x[j];
		}
		else
		{
			MPI_Recv(buf, M, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int j = 0; j < M; j++)
				full_root[i*M + j] = buf[j];
		}
	}
	free(buf);
}
void countDiff(double* diff, double *const *const A, double *const b, double *const x, double *const full_root)
{
	//count Ax
	for (int i = 0; i < M; i++)
	{
		diff[i] = 0;
		for (int j = 0; j < N; j++)
			diff[i] += A[i][j] * full_root[j];
	}
	//count Ax - b
	for (int i = 0; i < M; i++)
		diff[i] -= b[i];
}
void refreshRoot(double* root, double *const diff, const double t)
{
	//root = root - t * diff
	for (int i = 0; i < M; i++)
		root[i] = root[i] - t * diff[i];
}
double* countRoot(double *const *const A, double *const b)
{
	double* x = (double*)calloc(M, sizeof(double));
	double* full_root = (double*)calloc(N, sizeof(double));
	double* diff = (double*)calloc(M, sizeof(double));
	const double t = 0.0001;
	const double epsilon = 0.001;
	while(1)
	{
		//get parts of root from other processes
		getFullRoot(full_root, x);
		//count Ax - b
		countDiff(diff, A, b, x, full_root);
		//count metric and exit if root is good
		double my_metric = metric(diff);
		double reduced_metric = 0;
		MPI_Allreduce(&my_metric, &reduced_metric, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (reduced_metric < epsilon)
		{
			//if root is good
			free(diff);
			free(full_root);
			return x;
		}
		else
		{
			//if root is bad
			refreshRoot(x, diff, t);
		}
	}
}
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);//start work of MPI
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get nuber of process
		MPI_Comm_size(MPI_COMM_WORLD, &size); // get count of processes
		//init part of matrix
		double** const A = createMatrix();
		for (int i = 0; i < M; i++)
			A[i][rank*M + i] = 1;
		//init part of b vector
		double* b = (double*)calloc(M, sizeof(double));
		for (int i = 0; i < M; i++)
		{
			b[i] = 1;
		}
		//get part of root
		double* root = countRoot(A, b);
		//print root and clear memory
		if (rank == 0)
		{
			printf("Root has finished\n");
			printVector(root);
		}
		deleteMatrix(A);
		free(b);
		free(root);
	}
	MPI_Finalize(); // end work of MPI
	return 0;
}