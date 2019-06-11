#pragma comment(lib,"mpi.lib") 
#include"mpi.h"
#include<iostream>
#include"windows.h"
#include<cmath>
#include <time.h>
using namespace std;
double f1(double a)
{
	return (4.0 / (1.0 + a*a));
}
double f2(double a)
{
	if (int(a) % 2 == 1) {
		return (-4) / (2 * a + 1);
	}
	else {
		return 4 / (2 * a + 1);
	}
}
double f3(double a)
{
	if (int(a) % 2 == 1) {
		return (-4) / ((2 * a + 1)*pow(5,int(2*a+1)))-(-1)/((2*a+1)*pow(239,int(2*a+1)));
	}
	else {
		return 4 / ((2 * a + 1)*pow(5, int(2 * a + 1))) - 1 / ((2 * a + 1)*pow(239, int(2 * a + 1)));
	}
}
int main(int argc, char* argv[])
{
	int done = 0, n, myid, numprocs, i, rc;
	double PI25DT = 3.141592653589793238462643;
	double mypi, pi, h, sum, x, a;
	double starttime=0.0, endtime = 0.0;

	MPI_Status status;
	MPI_Request request;
	MPI_Init(&argc, &argv); /*MPI的初始化函数*/
	MPI_Comm_rank(MPI_COMM_WORLD, &myid); /*该进程的编号*/
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs); /*总的进程数目*/


	if (myid == 0) {
		n = 10000;
		starttime = MPI_Wtime();
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/* f1(x)面积积分法
	h = 1.0 / (double)n;   // width of each sub
	sum = 0.0;
	for (i = myid + 1; i <= n; i += numprocs)  
	{
		x = h * ((double)i - 0.5);
		sum += f1(x);
	}
	mypi = h * sum;
	*/

	//幂级数法与改进幂级数法
	/*
	sum = 0.0;
	for (i = myid; i <= n; i += numprocs)
	{
		sum += f3(i);//f(2),f(3)
	}
	mypi =4*sum/n;
	*/
	
	//蒙特卡洛法
	double count = 0;
	double num = 0;
	srand((unsigned)time(NULL));
	for (i = myid; i <= n; i += numprocs) {
		double pointx = double(rand() / double(RAND_MAX));
		double pointy = double(rand() / double(RAND_MAX));
		double distance = pointx*pointx + pointy*pointy;
		if (distance <= 1.0) {
			count += 1;
		}
		num += 1;
	}
	double pi_one = double(count / n);
	mypi = double(4.0*pi_one);
	
	//随机积分法
	/*
	double one_f = 0.0;
	srand((unsigned)time(NULL));
	for (i = myid; i <= n; i += numprocs)
	{
		double x = double(rand() / double(RAND_MAX));
		double f = f1(x);
		double x2 = 1.0;
		double x1 = 0.0;
		one_f += f*(x2 - x1);
	}
	mypi = double(one_f / n);
	*/

	MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (myid == 0) {
		endtime = MPI_Wtime();//
		printf("pi is approximately %.24f, Error is %.24f\n", pi, fabs(pi - PI25DT));
		cout << "time is " << endtime - starttime << endl;
		getchar();
	}
	MPI_Finalize();

	return 0;
}