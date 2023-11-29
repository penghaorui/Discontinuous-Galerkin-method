/*this is the main program for the MPI version of the DG code, the example model here is a simple two layer geologic model
 with a slope. For this example, each core is assigned with one partition. Here, two cores are used in total.
 Each core has a copy of all the Vx, Vz,and stress variables and the broadcast of MPI transfers computation of the other cores*/
#include "dg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <malloc.h>
#include <omp.h>
#include "fun_save_to_file.c"
#include "fun_DG.c"
#include "MPIFunction.c"
#include "memory.h"
#include <mpi.h>

#include <sys/time.h>
double dtime()
{
	double tseconds = 0.0;
	struct timeval mytime;
	gettimeofday(&mytime, (struct timezone *)0);
	tseconds = (double)(mytime.tv_sec + mytime.tv_usec * 1.0e-6);
	return tseconds;
}

int main(int argc, char *argv[])
{	
	/*initializing MPI functions*/
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &size);            

	/*claim time variables and counters*/
	double start, finish;
	double duration = 0.0;
	long int i, j, k, n;
	long int num;
	double sum;

	// read mesh data
	char *path_triangle;
	char *path_point;
	char *path_neighbor;
	char *path;
	path_triangle = "input_files/mesh/slope_model_converted.1.ele";
	path_neighbor = "input_files/mesh/slope_model_converted.1.neigh";
	path_point = "input_files/mesh/slope_model_converted.1.node";

	// read precomputed DG matrice
	printf("*** read DG matrice\n");
	path = "input_files/precomputed_matrix/FeketeNew.txt";
	static double FeketeNew[Np][2] = {{0}};
	ReadDouble(path, *FeketeNew, Np, 2);
	path = "input_files/precomputed_matrix/S1.txt";
	static double S1[Np][N + 1] = {{0}};
	ReadDouble(path, *S1, Np, N + 1);
	path = "input_files/precomputed_matrix/S2.txt";
	static double S2[Np][N + 1] = {{0}};
	ReadDouble(path, *S2, Np, N + 1);
	path = "input_files/precomputed_matrix/S3.txt";
	static double S3[Np][N + 1] = {{0}};
	ReadDouble(path, *S3, Np, N + 1);
	path = "input_files/precomputed_matrix/V.txt";
	static double V[Np][Np] = {{0}};
	ReadDouble(path, *V, Np, Np);
	path = "input_files/precomputed_matrix/Vb.txt";
	static double Vb[N + 1][N + 1] = {{0}};
	ReadDouble(path, *Vb, N + 1, N + 1);
	path = "input_files/precomputed_matrix/Vr.txt";
	static double Vr[Np][Np] = {{0}};
	ReadDouble(path, *Vr, Np, Np);
	path = "input_files/precomputed_matrix/Vs.txt";
	static double Vs[Np][Np] = {{0}};
	ReadDouble(path, *Vs, Np, Np);
	path = "input_files/precomputed_matrix/inverV.txt";
	static double inverV[Np][Np] = {{0}};
	ReadDouble(path, *inverV, Np, Np);
	path = "input_files/precomputed_matrix/inverVb.txt";
	static double inverVb[N + 1][N + 1] = {{0}};
	ReadDouble(path, *inverVb, N + 1, N + 1);

	// Read mesh size
	long Tri = Calculate_Tri(path_triangle);
	long Point_Num = Calculate_Point_Num(path_point);

	// claim material variable: Vp, Density and modulus
	double *Vp;
	Vp = (double *)malloc(sizeof(double) * Tri);
	double *Den;
	Den = (double *)malloc(sizeof(double) * Tri);
	double *Miu;
	Miu = (double *)malloc(sizeof(double) * Tri);
	double *Lambda;
	Lambda = (double *)malloc(sizeof(double) * Tri);

	/*Node and neighboring edge information of the triangular mesh*/
	long *Node_Neighbor;
	Node_Neighbor = (long *)malloc(sizeof(long) * Tri * 3 * (N + 1));
	long *Neighbor;
	Neighbor = (long *)malloc(sizeof(long) * Tri * 3);
	long *Edge;
	Edge = (long *)malloc(sizeof(long) * Tri * 3);

	/*flag of the boundary edge for the whole model*/
	long *Edge_flag;
	Edge_flag = (long *)malloc(sizeof(long) * Tri * 3);

	/*absorbing function*/
	double *absorb_fun;
	absorb_fun = (double *)malloc(sizeof(double) * Tri * Np);

	/*neighbor mark to check if the neighboring triangle is on the edge of the partitioned area*/
	long *Neighbor_Mark;
	Neighbor_Mark = (long *)malloc(sizeof(long) * Tri * 3);

	/*number of triangles in each partitioned area*/
	long *Series;
	Series = (long *)malloc(sizeof(long) * Process_Num);

	/*reordering from zero for the number of triangles in each partitioned area*/
	long *Series_order;
	Series_order = (long *)malloc(sizeof(long) * Process_Num);

	/*number of triangles on the border of th partitioned area */
	long *Bound_Series;
	Bound_Series = (long *)malloc(sizeof(long) * Process_Num);

	/*triangle partition information*/
	long *Triangle_part;
	Triangle_part = (long *)malloc(sizeof(long) * Tri);

	/*Jacobi matrix of the triangles and the three edges*/
	double *Jacobi;
	Jacobi = (double *)malloc(sizeof(double) * Tri);
	double *Jacobi_Edge;
	Jacobi_Edge = (double *)malloc(sizeof(double) * Tri * 3);

	/*matrix for the central flux*/
	double *Q;
	Q = (double *)malloc(sizeof(double) * Tri * 3 * 5 * 2);

	/*matrix used in the DG theory*/
	double *A;
	A = (double *)malloc(sizeof(double) * Tri * 5);
	double *B;
	B = (double *)malloc(sizeof(double) * Tri * 5);
	double *Dx;
	Dx = (double *)malloc(sizeof(double) * Tri * Np * Np);
	double *Dz;
	Dz = (double *)malloc(sizeof(double) * Tri * Np * Np);
	double *S;
	S = (double *)malloc(sizeof(double) * 3 * Np * (N + 1));
	long *ibe;
	ibe = (long *)malloc(sizeof(long) * 3 * (N + 1));

	/*normal vectors of the three edges of each triangle*/
	double *normal1;
	normal1 = (double *)malloc(sizeof(double) * Tri * 2);
	double *normal2;
	normal2 = (double *)malloc(sizeof(double) * Tri * 2);
	double *normal3;
	normal3 = (double *)malloc(sizeof(double) * Tri * 2);

	/*variable to stored recorded two-component data at given location*/
	double *record_x;
	record_x = (double *)malloc(sizeof(double) * TIME * Channel);
	memset(record_x, 0, sizeof(double) * TIME * Channel);
	double *record_z;
	record_z = (double *)malloc(sizeof(double) * TIME * Channel);
	memset(record_z, 0, sizeof(double) * TIME * Channel);

	/*variable for the source configuration*/
	double *fun;
	fun = (double *)malloc(sizeof(double) * Source_tri * Np);

	// Input mesh: nodes, triangles and neighboring triangles
	double *Point;
	Point = (double *)malloc(sizeof(double) * Point_Num * 2);
	long *Triangle;
	Triangle = (long *)malloc(sizeof(long) * Tri * 3);
	long *neighbor;
	neighbor = (long *)malloc(sizeof(long) * Tri * 3);

	/*numbered physical attributes of each triangle*/
	long *attribute;
	attribute = (long *)malloc(sizeof(long) * Tri);

	printf("*** read mesh\n");
	Mesh_Input(Tri, Point_Num, path_triangle, path_neighbor, path_point, Point, Triangle, neighbor, attribute);

	// set neighbor nodes
	NeighborNodes(Tri, Node_Neighbor, Triangle, Edge, Neighbor, ibe);

	/*read the partition information*/
	path = "input_files/mesh/slope_partition.mesh.epart.2";
	long *epart;
	epart = (long *)malloc(sizeof(long) * Tri);
	epartinput(path, epart, Tri);

	// calculate normal vector and Jacobi matrix
	Normals(Tri, Point, Triangle, normal1, normal2, normal3);
	Jacobian(Tri, Jacobi_Edge, Jacobi, Point, Triangle);

	// set source and receiver
	printf("\n*** get  source\n\n");
	double source_position[2] = {1200, 1200};
	long *source;
	source = (long *)malloc(sizeof(long) * Source_tri * 2);
	Get_Source(Tri, Point_Num, Point, Triangle, source, source_position);

	double receiver_position[2] = {1600, 1600};
	long *receiver;
	receiver = (long *)malloc(sizeof(long) * Channel * 2);
	printf("\n*** get receiver\n\n");

	Get_Receiver(Tri, Point_Num, Point, Triangle, receiver, receiver_position);


	// set source by DG theory
	for (num = 0; num < Source_tri; num++)
	{
		for (i = 0; i < Np; i++)
		{
			n = source[num * 2 + 0] - 1;
			k = source[num * 2 + 1] - 1;
			sum = 0;
			for (j = 0; j < Np; j++)
			{
				sum = sum + V[i][j] * V[k][j];
			}
			fun[num * Np + i] = sum / Jacobi[n];
		}
	}

	// set model property
	double Cp1 = 1500;
	double Cs1 = 0;
	double den1 = 1000;
	double miu1 = den1 * Cs1 * Cs1;
	double lambda1 = den1 * Cp1 * Cp1 - 2 * miu1;

	double Cp2 = 3000;
	double Cs2 = 1800;
	double den2 = 2200;
	double miu2 = den2 * Cs2 * Cs2;
	double lambda2 = den2 * Cp2 * Cp2 - 2 * miu2;

	for (i = 0; i < Tri; i++)
	{
		if (attribute[i] == 1 || attribute[i] == 2)
		{
			Vp[i] = Cp1;
			Lambda[i] = lambda1;
			Miu[i] = miu1;
			Den[i] = den1;
		}

		if (attribute[i] == 3 || attribute[i] == 4)
		{
			Vp[i] = Cp2;
			Lambda[i] = lambda2;
			Miu[i] = miu2;
			Den[i] = den2;
		}
	}

	// set time step dt and ricker source delay time t0 by the peak frequency Fp
	double dt = Calculate_Time_Step(Tri, Point, Triangle, Vp);
	double t0;
	t0 = Ricker_Delay_Time(dt);
	printf("step dt  %f delay t0 %f\n", dt, t0);

	/*get the total number of triangle inside and outside the absorbing boundary*/
	long Absorb_Num;
	Absorb_Num = Get_Absorb_Square_Num(Tri, Point, Triangle);
	long Inside_Num;
	Inside_Num = Get_Inside_Square_Num(Tri, Point, Triangle);

	/*set the triangles in the absorbing boundary*/
	long *absorbing_triangle;
	absorbing_triangle = (long *)malloc(sizeof(long) * Absorb_Num);
	long *inside_triangle;
	inside_triangle = (long *)malloc(sizeof(long) * Inside_Num);
	Get_absorb_square(Tri, Point, Triangle, absorbing_triangle, inside_triangle);
	
	/*set the absorbing function*/
	Absorbing_Fun(Tri, Absorb_Num, 0.1, absorb_fun, absorbing_triangle, Point, Triangle, FeketeNew);

	/*set matrix A and B that represent the material properties by DG theory*/
	Elastic_A_B(Tri, A, B, Lambda, Miu, Den);

	/*set matrix S by DG theory*/
	for (j = 0; j < Np; j++)
	{
		for (n = 0; n < N + 1; n++)
		{
			S[0 * Np * (N + 1) + j * (N + 1) + n] = S1[j][n];
			S[1 * Np * (N + 1) + j * (N + 1) + n] = S2[j][n];
			S[2 * Np * (N + 1) + j * (N + 1) + n] = S3[j][n];
		}
	}

	/*set neighboring triangles and edge odering*/
	Neighbor_Edge(Tri, Neighbor, neighbor, Edge, Triangle, Edge_flag);
	IBE(ibe);
	
	/*set neighboring nodes*/
	NeighborNodes(Tri, Node_Neighbor, Triangle, Edge, Neighbor, ibe);

	//compute Dx and Dz matrix by DG theory
	Dx_Dz(Tri, Dx, Dz, Vr, Vs, inverV, Triangle, Point);

	/*set the neighbor partition mark that is used to determine if the neighbor triangle is in a different partition*/
	Neighbor_Partition_Mark(Tri, Neighbor, epart, Neighbor_Mark);

	// calculate boudnary triangle number and set boundary triangles that separate different partitions 
	long int Boundary_Num;
	Boundary_Num = Calculate_Boundary_Num(Tri, Neighbor_Mark);

	long int *Boundary_Triangle;
	Boundary_Triangle = (long *)malloc(sizeof(long) * Boundary_Num);

	// find boudnary triangle and neighbor mark and reordered  triangles
	Partition_Boundary_Neighbor(Tri, Neighbor, epart, Neighbor_Mark, Series, Series_order, Bound_Series, Triangle_part, Boundary_Triangle);

	//compute the flux by DG theory
	Central_flux(Tri, Jacobi, Jacobi_Edge, Q, A, B, normal1, normal2, normal3);

	//reoder the boundary triangles for each partition from zero
	long *Bound_Series_order;
	Bound_Series_order = (long *)malloc(sizeof(long) * Process_Num);
	j = 0;
	for (i = 0; i < Process_Num; i++)
	{
		if (rank == 0)
		{
			printf("Series[%ld] = %ld\t ", i, Series[i]);
			printf("Series_order[%ld] = %ld\n ", i, Series_order[i]);
		}

		Bound_Series_order[i] = j;
		j += Bound_Series[i];
		if (rank == 0)
		{
			printf("Bound_Series[%ld] = %ld\t", i, Bound_Series[i]);
			printf("Bound_Series_order[%ld] = %ld\n", i, Bound_Series_order[i]);
		}
	}
	printf("Tri=%ld\tbound=%ld\n", Tri, Boundary_Num);

	/*save the node information for plotting the wavefield snapshot using matlab later*/
	double *All_node;
	All_node = (double *)malloc(sizeof(double) * Tri * Np * 2);
	set_All_node(Tri, All_node, Point, Triangle, FeketeNew);
	path = "output_files/All_node.txt";
	Save_Data(path, All_node, Tri * Np, 2);

	/*set MPI barrier for synchronization*/
	fflush(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	/*free variables that are not used later*/
	free(neighbor);
	free(Point);
	free(Triangle);
	free(normal1);
	free(normal2);
	free(normal3);
	free(Lambda);
	free(Miu);
	free(Den);
	free(epart);

	/*variables of x and z component of particle velocities of all triangles that are gathered from every core when
	trying to save to files. They are not part of the DG computation process*/
	double *Vx1_recv;
	Vx1_recv = (double *)malloc(sizeof(double) * Tri * Np);

	double *Vz1_recv;
	Vz1_recv = (double *)malloc(sizeof(double) * Tri * Np);

	/*variables of x and z component of particle velocities that are used for computing*/
	double *Vx1;
	Vx1 = (double *)malloc(sizeof(double) * Tri * Np);
	double *Vz1;
	Vz1 = (double *)malloc(sizeof(double) * Tri * Np);

	/*variables of three components of stresses that are used for computing*/
	double *Txx1;
	Txx1 = (double *)malloc(sizeof(double) * Tri * Np);
	double *Tzz1;
	Tzz1 = (double *)malloc(sizeof(double) * Tri * Np);
	double *Txz1;
	Txz1 = (double *)malloc(sizeof(double) * Tri * Np);

	/*variable for computing Ricker wavelet as the source time function*/
	double w;

	/*variables Vx and Vz of boundary triangles (between partitions) that are sent from every core before
	Vx and Vz are used to compute stresses, they are not part of the DG computation process*/
	double *Vx1_b;
	Vx1_b = (double *)malloc(sizeof(double) * Boundary_Num * Np);
	double *Vz1_b;
	Vz1_b = (double *)malloc(sizeof(double) * Boundary_Num * Np);

	/*variables Vx and Vz of boundary triangles (between partitions) that receive broadcast from every core before
	Vx and Vz are used to compute stresses, they are not part of the DG computation process*/
	double *Vx1_b_r;
	Vx1_b_r = (double *)malloc(sizeof(double) * Boundary_Num * Np);
	double *Vz1_b_r;
	Vz1_b_r = (double *)malloc(sizeof(double) * Boundary_Num * Np);

	/*variables stress of boundary triangles (between partitions) that are sent from every core before
	stresses are used to compute Vx and Vz, they are not part of the DG computation process*/
	double *Txx1_b;
	Txx1_b = (double *)malloc(sizeof(double) * Boundary_Num * Np);
	double *Tzz1_b;
	Tzz1_b = (double *)malloc(sizeof(double) * Boundary_Num * Np);
	double *Txz1_b;
	Txz1_b = (double *)malloc(sizeof(double) * Boundary_Num * Np);

	/*variables stress of boundary triangles (between partitions) that receive broadcast from every core before
	stresses are used to compute Vx and Vz, they are not part of the DG computation process*/
	double *Txx1_b_r;
	Txx1_b_r = (double *)malloc(sizeof(double) * Boundary_Num * Np);
	double *Tzz1_b_r;
	Tzz1_b_r = (double *)malloc(sizeof(double) * Boundary_Num * Np);
	double *Txz1_b_r;
	Txz1_b_r = (double *)malloc(sizeof(double) * Boundary_Num * Np);

	/*counting and ordering of variables (Vx,Vz and stresses) for the boundary triangles (between partitions)*/
	int *Bound_Series_sr;
	Bound_Series_sr = (int *)malloc(sizeof(int) * Process_Num);
	int *Bound_Series_order_sr;
	Bound_Series_order_sr = (int *)malloc(sizeof(int) * Process_Num);

	for (i = 0; i < Process_Num; i++)
	{
		Bound_Series_order_sr[i] = (int)Bound_Series_order[i] * Np;
		Bound_Series_sr[i] = (int)Bound_Series[i] * Np;
	}

	
	int m;
	long time;

	/*start computation*/
	start = dtime();
	for (m = 0; m < 1; m++)
	{
		memset(Vx1, 0, sizeof(double) * Np * Tri);
		memset(Vz1, 0, sizeof(double) * Np * Tri);
		memset(Txx1, 0, sizeof(double) * Np * Tri);
		memset(Tzz1, 0, sizeof(double) * Np * Tri);
		memset(Txz1, 0, sizeof(double) * Np * Tri);

		if (rank == 0)
		{
			printf("source %ld\n", source[m * 2]);
		}

		for (time = 0; time < TIME; time++) // time stepping
		{

			int send_size = Bound_Series[rank] * Np;

			/*get source time function*/
			w = Ricker_Source(time, dt, t0);

			/*apply absorbing boundary conditions*/
			Absorb_Boundary(rank, Tri, Vx1, Vz1, Txx1, Tzz1, Txz1, absorb_fun, Triangle_part, Series, Series_order);

			/*apply source to the designated triangle*/
			source_kernel(Vz1, source, fun, w, m);

			/*gather stress of boundary triangles (between partitions)*/
			for (i = Bound_Series_order[rank]; i < Bound_Series_order[rank] + Bound_Series[rank]; i++)
			{
				int num = Boundary_Triangle[i] - 1;
				for (j = 0; j < Np; j++)
				{
					int id_b = i * Np + j;
					int id_num = num * Np + j;
					Txx1_b[id_b] = Txx1[id_num];
					Tzz1_b[id_b] = Tzz1[id_num];
					Txz1_b[id_b] = Txz1[id_num];
				}
			}

			/*send stress of boundary triangles (between partitions) to every core*/
			double *txx_send = &Txx1_b[Bound_Series_order_sr[rank]];
			double *tzz_send = &Tzz1_b[Bound_Series_order_sr[rank]];
			double *txz_send = &Txz1_b[Bound_Series_order_sr[rank]];
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgatherv(txx_send, send_size, MPI_DOUBLE, Txx1_b_r, Bound_Series_sr, Bound_Series_order_sr, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Allgatherv(tzz_send, send_size, MPI_DOUBLE, Tzz1_b_r, Bound_Series_sr, Bound_Series_order_sr, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Allgatherv(txz_send, send_size, MPI_DOUBLE, Txz1_b_r, Bound_Series_sr, Bound_Series_order_sr, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			
			/*receive stress of boundary triangles (between partitions) from every core*/
			for (i = 0; i < Boundary_Num; i++)
			{
				int num = Boundary_Triangle[i] - 1;
				for (j = 0; j < Np; j++)
				{
					int id_b = i * Np + j;
					int id_num = num * Np + j;
					Txx1[id_num] = Txx1_b_r[id_b];
					Tzz1[id_num] = Tzz1_b_r[id_b];
					Txz1[id_num] = Txz1_b_r[id_b];
				}
			}

			/*use stress to update Vx and Vz*/
			LU_Velocity_kernel(rank, Tri, dt, Vx1, Vz1, Txx1, Tzz1, Txz1, Triangle_part, Neighbor_Mark, Series, Series_order, A, B, Dx, Dz, Q, S, ibe, Node_Neighbor, Neighbor, Edge_flag);

			/*gather Vx and Vz of boundary triangles (between partitions)*/
			for (i = Bound_Series_order[rank]; i < Bound_Series_order[rank] + Bound_Series[rank]; i++)
			{
				int num = Boundary_Triangle[i] - 1;
				for (j = 0; j < Np; j++)
				{
					int id_b = i * Np + j;
					int id_num = num * Np + j;
					Vx1_b[id_b] = Vx1[id_num];
					Vz1_b[id_b] = Vz1[id_num];
				}
			}

			/*send Vx and Vz of boundary triangles (between partitions) to every core*/
			double *vx_send = &Vx1_b[Bound_Series_order_sr[rank]];
			double *vz_send = &Vz1_b[Bound_Series_order_sr[rank]];
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgatherv(vx_send, send_size, MPI_DOUBLE, Vx1_b_r, Bound_Series_sr, Bound_Series_order_sr, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Allgatherv(vz_send, send_size, MPI_DOUBLE, Vz1_b_r, Bound_Series_sr, Bound_Series_order_sr, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

			/*receive Vx and Vz of boundary triangles (between partitions) from every core*/
			for (i = 0; i < Boundary_Num; i++)
			{
				int num = Boundary_Triangle[i] - 1;
				for (j = 0; j < Np; j++)
				{
					int id_b = i * Np + j;
					int id_num = num * Np + j;
					Vx1[id_num] = Vx1_b_r[id_b];
					Vz1[id_num] = Vz1_b_r[id_b];
				}
			}

			/*use Vx and Vz to update stress*/
			LU_Stress_Elastic_kernel(rank, Tri, dt, Txx1, Tzz1, Txz1, Vx1, Vz1, Triangle_part, Neighbor_Mark, Series, Series_order, A, B, Dx, Dz, Q, S,
									 ibe, Node_Neighbor, Neighbor, Edge_flag);

			
			/*save point record of the wavefield*/
			MPI_Barrier(MPI_COMM_WORLD);
			Record(Vx1,Vz1,time,record_x,record_z,receiver);	

			finish = dtime();
			duration = (double)(finish - start);
			if (rank == 0 && (time + 1) % 10 == 0)
			{
				printf("time:%ld,\tduration is %f s\n", time, duration);
			}

			if ((time + 1) % 100 == 0)
			{	
				/*gather Vx and Vz from every core to the core of rank 0*/
				MPI_Barrier(MPI_COMM_WORLD);
				if (rank != 0)
				{
					MPI_Send(Vx1, Tri * Np, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
					MPI_Send(Vz1, Tri * Np, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
				}
				else
				{
					for (i = 1; i < size; i++)
					{
						MPI_Recv(Vx1_recv, Tri * Np, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
						MPI_Recv(Vz1_recv, Tri * Np, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
						for (k = Series_order[i]; k < (Series_order[i] + Series[i]); k++)
						{
							num = Triangle_part[k] - 1;
							for (j = 0; j < Np; j++)
							{
								Vx1[num * Np + j] = Vx1_recv[num * Np + j];
								Vz1[num * Np + j] = Vz1_recv[num * Np + j];
							}
						}
					}
				}
				MPI_Barrier(MPI_COMM_WORLD);

				/*After the assembly, save Vx and Vz to files */
				if (rank == 0)
				{
					char filename[50];
					sprintf(filename, "output_files/Vx%ld", (time + 1));
					printf("%s\n", filename);
					Save_Data(filename, Vx1, Np * Tri, 1);

					sprintf(filename, "output_files/Vz%ld", (time + 1));
					printf("%s\n", filename);
					Save_Data(filename, Vz1, Np * Tri, 1);
				}
			}

		} // time end
	}

	free(record_x);
	free(record_z);

	MPI_Finalize();
	return (0);
}
