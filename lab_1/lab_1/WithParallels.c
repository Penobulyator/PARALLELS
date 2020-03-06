#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N  6
#define P 1
#define M N/P
double** createMatrix()
{
	double** out = (double**)malloc(M * sizeof(double*));
	for (int i = 0; i < M; i++)
		out[i] = (double*)calloc(N , sizeof(double));
	return out;
}
double* matrixMULTvectror(double** const matrix, double* const vector)
{
	double out[M] = {};// = createVector();
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			out[i] += matrix[i][j] * vector[i];
	return &(out[0]);
}
double* vectorMINUSvector(double* a, double* b)
{
	double out[M] = {};// = createVector();
	for (int i = 0; i < M; i++)
		out[i] = a[i] - b[i];
	return &(out[0]);
}
double metric(double* vector)
{
	double out = 0;
	for (int i = 0; i < M; i++)
		out += vector[i] * vector[i];
	return sqrt(out);
}
double* multScal(double* const vector, double num)
{
	double out[M] = {};// = createVector();
	for (int i = 0; i < M; i++)
		out[i] = vector[i] * num;
	return &(out[0]);
}
void printMatrix(double** matrix)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
			printf("%f ", matrix[i][j]);
		printf("\n");
	}
	printf("\n");
}
void printVector(double* vector)
{
	for (int i = 0; i < M; i++)
	{
		printf("%f ", vector[i]);
		printf("\n");
	}
	printf("\n");
}
double* get_full_root(double* my_root)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//send part of root to other processes
	for (int i = 0; i < P; i++)
		if (i != rank)
			MPI_Send(my_root, M, MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
	//receive missing parts from other processes of root and build a full root
	double output[N] = {};// = (double*)calloc(N, sizeof(double));
	double buf[M] = {};// = createVector();
	for (int i = 0; i < P; i++)
	{
		if (i == rank)
		{
			for (int j = 0; j < M; j++)
				output[rank*M + j] = my_root[j];
		}
		else
		{
			MPI_Recv(buf, M, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int j = 0; j < M; j++)
				output[i*M + j] = buf[j];
		}
	}
	//free(buf);
	return &(output[0]);
}

double* root(double** const A, double* const b)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получение номера процесса
	double my_root[M] = {};// = createVector();
	double t = 0.0001;
	const double epsilon = 0.001;
	for (int i = 0; i < 100; i++)
	{
		//count current difference between Ax and b
		double* full_root = get_full_root(my_root);
		if (rank == 0)
		{
			for (int j = 0; j < N; i++)
				printf("%d\n", full_root[j]);
			printf('\n');
		}
		double* mult = matrixMULTvectror(A, full_root);
		double* diff = vectorMINUSvector(mult, b);
		double my_metric = metric(diff) / metric(b);
		double reduced_metric = 0;
		MPI_Allreduce(&my_metric, &reduced_metric, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (reduced_metric < epsilon)
		{
			return my_root;
		}
		else
		{
			double* tdiff = multScal(diff, t);
			double* new_root = vectorMINUSvector(my_root, tdiff);
			my_root = new_root;
		}
	}

}
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	{
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получение номера процесса
		//init part of matrix
		double **my_matrix = createMatrix();
		for (int i = 0; i < M; i++)
			my_matrix[i][rank*M + i] = rank + i + 1;
		//init part of b vector
		double b[M];// = createVector();
		for (int i = 0; i < M; i++)
			b[i] = rank * M + i + 1;
		//get part of root
		if (rank == 0)
		{
			printf("Srarting root...\n");
		}
		double* x = root(my_matrix, b);
		if (rank == 0)
		{
			printVector(x);
		}
		free(b);
		free(x);
	}
	MPI_Finalize(); // Завершение работы MPI
	return 0;
}