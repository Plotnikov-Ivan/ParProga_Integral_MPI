#include <stdio.h>
#include "mpi.h"
#include <math.h>

/* Doc

func: y = 1/ln(x)

*/

double rectangleIntegration(int integration_method, int n, double a, double b);

int main(int argc, char* argv[])
{
	double a = 1, b = 400;
	int n = 80000000;
	double part_of_Square, part_of_Square1, part_of_Square2, total_Square1, total_Square2, total_Square = 0;
	double timeStart, timeEnd;

	int procNum, procRank;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	double procStep = (b - a) / procNum;
	double a1 = a + procStep * (double)procRank;
	double b1 = a + procStep * ((double)procRank + 1);


	if (procRank == procNum - 1)
		b1 = b;

	timeStart = MPI_Wtime(); //get time

	part_of_Square = rectangleIntegration(1, n / procNum, a1, b1);
	MPI_Reduce(&part_of_Square, &total_Square, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	part_of_Square1 = rectangleIntegration(2, n / procNum, a1, b1);
	MPI_Reduce(&part_of_Square1, &total_Square1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);


	part_of_Square2 = rectangleIntegration(3, n / procNum, a1, b1);
	MPI_Reduce(&part_of_Square2, &total_Square2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	timeEnd = MPI_Wtime();


	if (procRank == 0)
	{
		printf("Square right: %lf \n", total_Square);
		printf("Square middle: %lf \n", total_Square1);
		printf("Square left: %lf \n", total_Square2);
		printf("Time (in seconds): %lf \n", timeEnd - timeStart);
	}

	MPI_Finalize();

	return 0;
}

double rectangleIntegration(int integration_method, int n, double a, double b)
{
	double x, X, step = (b - a) / n;
	double result_Square = 0;
	//Rectangle methods
	if (integration_method == 1)
		x = a + step;
	else if (integration_method == 2)
		x = a + 3 * step / 2;
	else
		x = a + 2 * step;

	//calculating integral
	for (int i = 0; i < n; i++) {
		X = x + step * i;
		result_Square += 1 / log(X * step);
	}

	return result_Square;
}

