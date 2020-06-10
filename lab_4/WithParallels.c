#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

int size = -1;
int rank = -1;

#define N 32

#define Nx N 
#define Ny N
#define Nz N / size

#define D 2.0

#define hx (D / (Nx - 1))
#define hy (D / (Ny - 1))
#define hz (D / (Nz - 1))

int min_z;
int max_z;

#define eps 1e-8
#define a 0


double func(double i, double j, double k)
{
	double x = -1.0 + i * hx;
	double y = -1.0 + j * hy;
	double z = -1.0 + (min_z + k) * hz;

	return x * x + y * y + z * z;
}

double*** allocGrid() {
	double*** output = (double***)malloc(Nx * sizeof(double**));
	for (int i = 0; i < Nx; i++)
	{
		output[i] = (double**)malloc(Ny * sizeof(double*));
		for (int j = 0; j < Ny; j++)
			output[i][j] = (double*)calloc(Nz, sizeof(double));

	}
	return output;
}

void deleteGrid(double*** grid)
{
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
			free(grid[i][j]);
		free(grid[i]);
	}
	free(grid);
}

void setBourderyConditions(double*** f)
{
	//bourdery conditions for Z
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
		{
			f[i][j][0] = func(i, j, 0);
			f[i][j][Nz - 1] = func(i, j, Nz - 1);
		}
	//bourdery conditions for Y
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Nz; j++)
		{
			f[i][0][j] = func(i, 0, j);
			f[i][Ny - 1][j] = func(i, Ny - 1, j);
		}
	//bourdery conditions for X
	for (int i = 0; i < Ny; i++)
		for (int j = 0; j < Nz; j++)
		{
			f[0][i][j] = func(0, i, j);
			f[Nx - 1][i][j] = func(Nx - 1, i, j);
		}

}

double countResidual(double*** grid, double*** prevItarationGrid)
{
	double residual = 0;
	for (int i = 1; i < Nx - 2; i++)
		for (int j = 1; j < Ny - 2; j++)
			for (int k = 1; k < Nz - 2; k++)
			{
				double t = fabs(grid[i][j][k] - prevItarationGrid[i][j][k]);
				if (t > residual)
					residual = t;
			}
	return residual;
}

void exchangeValues(double*** grid, double*** prevItarationGrid)
{
	if (rank == size - 1)
	{
		for (int i = 1; i < Nx - 2; i++)
			for (int j = 1; j < Ny - 2; j++)
			{
				double buf = 0;
				//send upper layer
				MPI_Send(&prevItarationGrid[i][j][0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

				//receive upper layer
				MPI_Recv(&buf, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				grid[i][j][0] = 1 / (2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a) * (
					(prevItarationGrid[i + 1][j][0] + prevItarationGrid[i - 1][j][0]) / (hx * hx) +
					(prevItarationGrid[i][j + 1][0] + prevItarationGrid[i][j - 1][0]) / (hy * hy) +
					(prevItarationGrid[i][j][1] + buf) / (hz * hz)
					- 6);
			}
	}
	else if (rank == 0)
	{
		for (int i = 1; i < Nx - 2; i++)
			for (int j = 1; j < Ny - 2; j++)
			{
				double buf = 0;
				//send bottom layer
				MPI_Send(&prevItarationGrid[i][j][Ny - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
				//receive bottom layer
				MPI_Recv(&buf, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				grid[i][j][Nz - 1] = 1 / (2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a) * (
					(prevItarationGrid[i + 1][j][Nz - 1] + prevItarationGrid[i - 1][j][Nz - 1]) / (hx * hx) +
					(prevItarationGrid[i][j + 1][Nz - 1] + prevItarationGrid[i][j - 1][Nz - 1]) / (hy * hy) +
					(prevItarationGrid[i][j][Nz - 2] + buf) / (hz * hz)
					- 6);

			}
	}
	else
	{
		for (int i = 1; i < Nx - 2; i++)
			for (int j = 1; j < Ny - 2; j++)
			{
				double buf;
				//send upper layer
				MPI_Send(&prevItarationGrid[i][j][0], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

				//send bottom layer
				MPI_Send(&prevItarationGrid[i][j][Ny - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

				//receive upper layer
				MPI_Recv(&buf, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				grid[i][j][0] = 1 / (2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a) * (
					(prevItarationGrid[i + 1][j][0] + prevItarationGrid[i - 1][j][0]) / (hx * hx) +
					(prevItarationGrid[i][j + 1][0] + prevItarationGrid[i][j - 1][0]) / (hy * hy) +
					(prevItarationGrid[i][j][1] + buf) / (hz * hz)
					- 6);

				//receive bottom layer
				MPI_Recv(&buf, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				grid[i][j][Nz - 1] = 1 / (2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a) * (
					(prevItarationGrid[i + 1][j][Nz - 1] + prevItarationGrid[i - 1][j][Nz - 1]) / (hx * hx) +
					(prevItarationGrid[i][j + 1][Nz - 1] + prevItarationGrid[i][j - 1][Nz - 1]) / (hy * hy) +
					(prevItarationGrid[i][j][Nz - 2] + buf) / (hz * hz)
					- 6);
			}
	}
}
void refreshGrid(double*** grid, double*** prevItarationGrid)
{
	//struct timeval tv1, tv2;

	//gettimeofday(&tv1, NULL);
	exchangeValues(grid, prevItarationGrid);
	//gettimeofday(&tv2, NULL);

	/*double dt_sec = (tv2.tv_sec - tv1.tv_sec);
	double dt_usec = (tv2.tv_usec - tv1.tv_usec);
	double dt = dt_sec + 1e-6*dt_usec;
	if (rank == 0)
		printf("time diff: %e \n", dt);*/

	for (int i = 1; i < Nx - 2; i++)
		for (int j = 1; j < Ny - 2; j++)
			for (int k = 2; k < Nz - 2; k++)
			{
				grid[i][j][k] = 1 / (2 / (hx * hx) + 2 / (hy * hy) + 2 / (hz * hz) + a) * (
					(prevItarationGrid[i + 1][j][k] + prevItarationGrid[i - 1][j][k]) / (hx * hx) +
					(prevItarationGrid[i][j + 1][k] + prevItarationGrid[i][j - 1][k]) / (hy * hy) +
					(prevItarationGrid[i][j][k + 1] + prevItarationGrid[i][j][k - 1]) / (hz * hz)
					- 6);
			}
}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);//start work of MPI
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get number of process
		MPI_Comm_size(MPI_COMM_WORLD, &size); // get count of processes
		min_z =  rank * Nz;
		max_z = (rank + 1) * Nz;

		double*** f = allocGrid();
		setBourderyConditions(f);

		double residual = 1e8;
		double*** f1 = allocGrid();

		MPI_Barrier(MPI_COMM_WORLD);

		struct timeval tv1, tv2;
		gettimeofday(&tv1, NULL);
		while (residual > eps)
		{
			//fefresh f1
			refreshGrid(f1, f);
			//count residual
			double my_residual = countResidual(f1, f);
			MPI_Allreduce(&my_residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			//if (rank == 0)
				//printf("r: %2f\n", residual);
			//refresh f
			for (int i = 1; i < Nx - 2; i++)
				for (int j = 1; j < Ny - 2; j++)
					for (int k = 1; k < Nz - 2; k++)
						f[i][j][k] = f1[i][j][k];

		}
		gettimeofday(&tv2, NULL);
		double dt_sec = (tv2.tv_sec - tv1.tv_sec);
		double dt_usec = (tv2.tv_usec - tv1.tv_usec);
		double dt = dt_sec + 1e-6*dt_usec;
		if (rank == 0)
			printf("time diff: %e \n", dt);

		double err = 0;
		int x, y, z;
		for (int i = 1; i < Nx - 2; i++)
			for (int j = 1; j < Ny - 2; j++)
				for (int k = 1; k < Nz - 2; k++)
				{
					double t = fabs(func(i, j, k) - f[i][j][k]);
					if (t > err)
					{
						err = t;
					}
				}

		printf("err: %2f\n", err);
		deleteGrid(f1);
		deleteGrid(f);
	}
	MPI_Finalize(); // end work of MPI
	return 0;
}