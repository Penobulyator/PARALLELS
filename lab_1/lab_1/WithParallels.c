#include <mpi.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#define N  6
#define P 3
#define M N/P
double** createMatrix()
{
	double** out = (double**)malloc(M * sizeof(double*));
	for (int i = 0; i < N; i++)
		out[i] = (double*)calloc(N, sizeof(double));
	return out;
}
void deleteMatrix(double** const matrix)
{
	for (int i = 0; i < M; i++)
		free(matrix[i]);
	free(matrix);
}
double* createVector()
{
	double* out = (double*)calloc(M, sizeof(double));
	return out;
}
double* matrixMULTvectror(double** const matrix, double* const vector)
{
	double* out = (double*)calloc(M, sizeof(double));
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			out[i] += matrix[i][j] * vector[i];
	return out;
}
double* vectorMINUSvector(double* a, double* b)
{
	double* out = (double*)malloc(M * sizeof(double));
	for (int i = 0; i < M; i++)
		out[i] = a[i] - b[i];
	return out;
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
	double* out = (double*)malloc(M * sizeof(double));
	for (int i = 0; i < M; i++)
		out[i] = vector[i] * num;
	return out;
}
void printMatrix(double** matrix)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
			printf("%f ", matrix[i][j]);
		printf("\r\n");
	}
	printf("\r\n");
}
void printVector(double* vector)
{
	FILE* fout = NULL;
	fopen_s(&fout, "out.txt", "w");
	for (int i = 0; i < M; i++)
	{
		fprintf(fout,"%f ", vector[i]);
		fprintf(fout, "\r\n");
	}
	fclose(fout);
	printf("\r\n");
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
	double* output = (double*)calloc(N, sizeof(double));
	double* buf = createVector();
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
	free(buf);
	return output;
}
double* root(double** const A, double* const b)
{
	double* my_root = (double*)calloc(M, sizeof(double));
	double t = 0.0001;
	const double epsilon = 0.00001;
	while (1)
	{
		//count current difference between Ax and b
		double* full_root = get_full_root(my_root);
		double* mult = matrixMULTvectror(A, full_root);
		double* diff = vectorMINUSvector(mult, b);
		double my_metric = metric(diff) / metric(b);
		double reduced_metric; 
		MPI_Reduce(&my_metric, &reduced_metric, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if (reduced_metric < epsilon)
		{
			free(full_root);
			free(mult);
			free(diff);
			return my_root;
		}
		else
		{
			double* tdiff = multScal(diff, t);
			double* new_root = vectorMINUSvector(my_root, tdiff);
			my_root = new_root;
			free(full_root);
			free(tdiff);
			free(mult);
			free(diff);
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
		double** my_matrix = createMatrix();
		for (int i = 0; i < M; i++)
			my_matrix[rank*M + i] = rank + i + 1;
		//init part of b vector
		double* b = createVector();
		for (int i = 0; i < M; i++)
			b[i] = rank * M + i + 1;
		//get part of root
		double* x = root(my_matrix, b);
		if (rank == 0)
		{
			printVector(x);

		}
		//deleteMatrix(my_matrix);
		//free(b);
		//free(x);
	}
	MPI_Finalize(); // Завершение работы MPI
	return 0;
}