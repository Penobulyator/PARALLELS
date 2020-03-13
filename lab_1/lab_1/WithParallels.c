#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N  1000
#define P 4
#define M N/P
double A[M][N] = { {} };
double b[M] = {};
double x[M] = {};
double full_root[N] = {};
double diff[M] = {}; //Ax - b
const double t = 0.0001;
const double epsilon = 0.001;
int rank = 0;
double metric()
{
	double out = 0;
	for (int i = 0; i < M; i++)
		out += diff[i] * diff[i];
	return sqrt(out);
}
void printMatrix()
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
			printf("%lf ", A[i][j]);
		printf("\n");
	}
	printf("\n");
}
void printVector(double* vector)
{
	for (int i = 0; i < M; i++)
	{
		printf("%lf ", vector[i]);
		printf("\n");
	}
	printf("\n");
}
void getFullRoot()
{
	//send part of root to other processes
	for (int i = 0; i < P; i++)
		if (i != rank)
			MPI_Send(x, M, MPI_DOUBLE, i, 123, MPI_COMM_WORLD);
	//receive missing parts from other processes of root and build a full root
	double buf[M] = {};// = createVector();
	for (int i = 0; i < P; i++)
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
	//free(buf);
}
void countDiff()
{
	for (int i = 0; i < M; i++)
	{
		diff[i] = 0;
		for (int j = 0; j < N; j++)
			diff[i] += A[i][j] * full_root[j];
	}
	for (int i = 0; i < M; i++)
		diff[i] -= b[i];
}
void refreshRoot()
{
	for (int i = 0; i < M; i++)
		x[i] = x[i] - t * diff[i];
}
void countRoot()
{
	while(1)
	{
		/*if (rank == 0)
		{
			printf("x:\n");
			printVector(x);
		}
		if (rank == 0)
		{
			printf("Full root:\n");
			for (int j = 0; j < N; j++)
				printf("%2f\n", full_root[j]);
			printf("\n");
		}*/
		//count current difference between Ax and b
		getFullRoot();
		/*if (rank == 0)
		{
			printf("Fool root:\n");
			for (int j = 0; j < N; j++)
				printf("%2f\n", full_root[j]);
			printf("\n");
		}*/
		countDiff();
		/*if (rank == 0)
		{
			printf("Ax - b:\n");
			printVector(diff);
		}*/
		double my_metric = metric();
		double reduced_metric = 0;
		MPI_Allreduce(&my_metric, &reduced_metric, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		//printf("%2f\n", reduced_metric);
		if (reduced_metric < epsilon)
		{
			return;
		}
		else
		{
			refreshRoot();
			/*if (rank == 0)
			{
				for (int j = 0; j < N; j++)
					printf("%2f\n", x[j]);
				printf("\n");
			}*/
		}
	}
	return;
}
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получение номера процесса
		//init part of matrix
		for (int i = 0; i < M; i++)
			A[i][rank*M + i] = 1;
		//init part of b vector and x
		for (int i = 0; i < M; i++)
		{
			b[i] = 1;
			x[i] = 0;
		}
		for (int i = 0; i < N; i++)
		{
			full_root[i] = 0;
		}
		//get part of root
		/*if (rank == 0)
		{
			printf("Srarting root...\n");
			printf("My matrix:\n");
			printMatrix();
			printf("My b:\n");
			printVector(b);
		}*/
		countRoot();
		if (rank == 0)
		{
			printf("Root has finished\n");
			printVector(x);
		}
	}
	MPI_Finalize(); // Завершение работы MPI
	return 0;
}