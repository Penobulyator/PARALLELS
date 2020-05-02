#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#define N1 8
#define N2 4
#define N3 4

int rank = 0;
int size = 0;

int p1 = 0;
int p2 = 0;

#define A_STATUS 0
#define B_STATUS 1
#define C_STATUS 2

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

void setProcessesNumber()
{
	//find p1 and p2 such that:
	//p1 * p2 = rank
	//N1 % p1 = 0
	//N3 % p2 = 0
	for (int i = 2; i < size - 1; i++)
		if ((size - 1) % i == 0)
		{
			int prob_p1 = i;
			int prob_p2 = size / i;

			if (N1 % prob_p1 == 0 && N3 % prob_p2 == 0)
			{
				p1 = prob_p1;
				p2 = prob_p2;
				return;
			}
		}
	//if there are no suitable values
	printf("Bad number of processes\n");
	exit(0);
}

/*ROOT PROCESS FUNCTIONS*/
void sendChunks(double **A, double**B)	
{
	//A is matrix of size N1*N2
	//B is matrix of size N2*N3

	//send parts of matrix A
	for (int lineNumber = 0; lineNumber < N1; lineNumber++)
	{

		MPI_Send(A[lineNumber], N2, MPI_DOUBLE, lineNumber / p2 + 1, A_STATUS, MPI_COMM_WORLD);
		/*printf("Sending this to process number %d: ", lineNumber / p2 + 1);
		for (int i = 0; i < N2; i++)
			printf("%2f ", A[lineNumber][i]);
		printf("\n");*/
	}
	//send parts of matrix B
	for (int columnNumber = 0; columnNumber < p2; columnNumber++)
	{
		for (int lineNumber = 0; lineNumber < N2; lineNumber++)
		{
				MPI_Send(&B[lineNumber][(N3 / p2) * columnNumber], N3 / p2, MPI_DOUBLE, 1 + columnNumber*p1, B_STATUS, MPI_COMM_WORLD);
			/*printf("Sending this to process number %d: ", p1*columnNumber + 1);
			for (int i = 0; i < N3/p2; i++)
				printf("%2f ", B[lineNumber][(N3 / p2) * columnNumber + i]);
			printf("\n");*/
		}
	}

}

void receiveChunks(double **C)
{
	//C is matrix of size N1*N3
	//we are receiving N1/p1 chunks of size N3/p2 (from one process)
	for (int processNumber = 1; processNumber <= p1*p2; processNumber++)
	{
		for (int lineNumber = 0; lineNumber < N1 / p1; lineNumber++)
		{
			MPI_Recv(&C[((processNumber - 1) % p1)*p2 +lineNumber][(N3 / p2)*((processNumber - 1) / p1)], N3 / p2, MPI_DOUBLE, processNumber, C_STATUS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

}
/*BUISNESS PROCESSES FUNCTIONS*/
void mult(double** C, double **A, double **B)
{
	//A is matrix of size (N1/p1)*N2
	//B is matrix of size N2*(N3/p2)
	//C is matrix of size (N1/p1)*(N3/p2) (C = A*B)
	
	for (int i = 0; i < N1 / p1; ++i)
	{
		for (int j = 0; j < N3 / p2; ++j)
		{
			C[i][j] = 0;
			for (int k = 0; k < N2; ++k)
				C[i][j] += A[i][k] * B[k][j];
		}
	}
	//	printf("Process number %d has finished mult\n", rank);
}
void receiveOneChunk(double **A, double **B)
{
	//A is matrix of size (N1/p1)*N2
	//B is matrix of size N2*(N3/p2)

	//receive part of A matrix (from left)
	int providerProcess = 0;
	if (rank - p1 >= 0)
	{
		providerProcess = rank - p1;
	}
	for (int i = 0; i < N1 / p1; i++)
	{
		MPI_Recv(A[i], N2, MPI_DOUBLE, providerProcess, A_STATUS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	//send part of A matrix to next process (to the right)
	if (rank + p1 < size)
		for (int i = 0; i < N1 / p1; i++)
		{
			MPI_Send(A[i], N2, MPI_DOUBLE, rank + p1, A_STATUS, MPI_COMM_WORLD);
		}
	//printf("Process number %d received A from process number %d and sended A to process number %d\n", rank, providerProcess, rank + p1);

	//receive part of B matrix
	providerProcess = 0;
	if ((rank - 1) % p1 != 0)
		providerProcess = rank - 1;
	for (int i = 0; i < N2; i++)
	{
		MPI_Recv(B[i], N3 / p2, MPI_DOUBLE, providerProcess, B_STATUS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	//send part of B matrix to next process (downword)
	if (rank % p1 != 0)
	{
		for (int i = 0; i < N2; i++)
		{
			MPI_Send(B[i], N3 / p2, MPI_DOUBLE, rank + 1, B_STATUS, MPI_COMM_WORLD);
		}
		//printf("Process number %d received B from process number %d and sended B to process number %d\n", rank, providerProcess, rank + 1);
	}
}
void sendOneChunk(double** C)
{
	//A is matrix of size (N1/p1)*(N3/p2)
	for (int i = 0; i < N1 / p1; i++)
	{
		MPI_Send(C[i], N3 / p2, MPI_DOUBLE, 0, C_STATUS, MPI_COMM_WORLD);
	}

}
/*MAIN*/
int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);//start work of MPI
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get number of process
		MPI_Comm_size(MPI_COMM_WORLD, &size); // get count of processes

		setProcessesNumber();

		if (rank == 0)
		{

			double** A = createMatrix(N1, N2);
			double** B = createMatrix(N2, N3);
			double** C = createMatrix(N1, N3);

			
			for (int i = 0; i < N1; i++)
				for (int j = 0; j < N2; j++)
					A[i][j] = i * j + i + j ;

			for (int i = 0; i < N2; i++)
				for (int j = 0; j < N3; j++)
					B[i][j] = i * j + i + j + 1;
			printf("p1 = %d, p2 = %d\n", p1, p2);

			printf("A:\n");
			printMatrix(A, N1, N2);

			printf("B:\n");
			printMatrix(B, N2, N3);

			setProcessesNumber();

			
			sendChunks(A, B);

			receiveChunks(C);


			printf("C:\n");
			printMatrix(C, N1, N3);

			deleteMatrix(A, N1);
			deleteMatrix(B, N2);
			deleteMatrix(C, N1);
		}
		else
		{
			double** A = createMatrix(N1/p1, N2);
			double** B = createMatrix(N2, N3/p2);

			receiveOneChunk(A, B);

			double** C = createMatrix(N1 / p1, N3 / p2);
			mult(C, A, B);

			sendOneChunk(C);

			deleteMatrix(A, N1 / p1);
			deleteMatrix(B, N2);
			deleteMatrix(C, N1/p1);
			
		}
	}
	MPI_Finalize(); // end work of MPI

	return 0;
}