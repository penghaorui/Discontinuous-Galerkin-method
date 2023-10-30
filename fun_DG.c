/*this file contains functions that are both used in DG_leap_Frog_MPI.c and DG_leap_Frog_OpenMP.c,
including reading mesh files, configuring meshes (seting neighboring nodes, edges and so on), 
computing time step, absoribing boundary, get source and receiver, and computing DG matrix, Ricker source time function.*/

/*prototyping functions*/
long Calculate_Tri(char path_triangle[]);

long Calculate_Point_Num(char path_point[]);

void Mesh_Input(long Tri, long Point_Num, char path_triangle[], char path_neighbor[], char path_point[], double *Point,
				long *Triangle, long *neighbor, long *attribute);

void EdgeLength(long Tri, double *Point, long *Triangle, double *Length);

double Ricker_Delay_Time(double dt);

double InDiameter_Min(long Tri, double *Point, long *Triangle, double *Length);

double Calculate_Time_Step(long Tri, double *Point, long *Triangle, double *Vp);

long Get_Absorb_Square_Num(long Tri, double *Point, long *Triangle);

long Get_Inside_Square_Num(long Tri, double *Point, long *Triangle);

void Get_absorb_square(long Tri, double *Point, long *Triangle, long *absorbing_triangle, long *inside_triangle);

void Get_Source(long Tri, long Point_Num, double *Point, long *Triangle, long *source, double *source_position);

void Get_Receiver(long Tri, long Point_Num, double *Point, long *Triangle, long *receiver, double *receiver_position);

void Jacobian(long Tri, double *Jacobi_Edge, double *Jacobi, double *Point, long int *Triangle);

void Angle(long Tri, double *Point, long int *Triangle, double *angle);

void Normals(long Tri, double *Point, long int *Triangle, double *normal1, double *normal2, double *normal3);

void Neighbor_Edge(long Tri, long int *Neighbor, long int *neighbor, long int *Edge, long int *Triangle, long int *Edge_flag);

void Neighbor_Partition_Mark(long Tri, long int *Neighbor, long int *epart, long int *Neighbor_Mark);

long Calculate_Boundary_Num(long Tri, long int *Neighbor_Mark);

void Partition_Boundary_Neighbor(long Tri, long int *Neighbor, long int *epart, long int *Neighbor_Mark, long *Series,
								 long *Series_order, long *Bound_Series, long *Triangle_part, long *Boundary_Triangle);

void set_All_node(long Tri, double *All_node, double *Point, long *Triangle, double FeketeNew[][2]);

void Rstoxy(double position[][2], double vertice[][2], double FeketeNew[][2]);

void Absorb_Tri_Node(long Absorb_Num, double *Absorb_tri_nodes, long int *Absorb_tri, double *Point, long int *Triangle, double FeketeNew[][2]);

void Absorbing_Fun(long Tri, long Absorb_Num, double beta, double *absorb_fun, long int *Absorb_tri, double *point,
				   long int *triangle, double FeketeNew[][2]);

void IBE(long int *ibe);

void Elastic_A_B(long Tri, double *A, double *B, double *Lambda, double *Miu, double *Den);

void Central_flux(long Tri, double *Jacobi, double *Jacobi_Edge, double *Q, double *A, double *B, double *normal1, double *normal2, double *normal3);

void Dx_Dz(long Tri, double *Dx, double *Dz, double Vr[][Np], double Vs[][Np], double inverV[][Np], long int *Triangle, double *Point);

void NeighborNodes(long Tri, long int *Node_Neighbor, long int *Triangle, long int *Edge, long int *Neighbor, long int *ibe);

double Ricker_Source(long int time, double dt, double t0);

/*****************************************introducing functions**********************************/

/*get total triangle number*/
long Calculate_Tri(char path_triangle[])
{
	long *Tri_Num;
	Tri_Num = (long *)malloc(sizeof(long) * 1);
	ReadInt(path_triangle, Tri_Num, 1, 1);
	printf("Triangle number %ld\n", *Tri_Num);
	return (*Tri_Num);
}

/*get total node number*/
long Calculate_Point_Num(char path_point[])
{
	long *P_Num;
	P_Num = (long *)malloc(sizeof(long) * 1);
	ReadInt(path_point, P_Num, 1, 1);
	printf("Point number %ld\n", *P_Num);
	return (*P_Num);
}

/*read mesh files*/
void Mesh_Input(long Tri, long Point_Num, char path_triangle[], char path_neighbor[], char path_point[], double *Point, long *Triangle, long *neighbor, long *attribute)
{
	long i, j;

	long *TRIANGLE;
	TRIANGLE = (long *)malloc(sizeof(long) * (Tri + 1) * 5);
	ReadInt(path_triangle, TRIANGLE, (Tri + 1), 5);
	long *NEIGHBOR;
	NEIGHBOR = (long *)malloc(sizeof(long) * (Tri + 1) * 4);
	ReadInt(path_neighbor, NEIGHBOR, (Tri + 1), 4);
	double *POINT;
	POINT = (double *)malloc(sizeof(double) * (Point_Num + 1) * 4);
	ReadDouble(path_point, POINT, (Point_Num + 1), 4);

	for (i = 0; i < Tri; i++)
	{
		for (j = 1; j <= 3; j++)
		{
			Triangle[i * 3 + j - 1] = TRIANGLE[(i + 1) * 5 + j];
			neighbor[i * 3 + j - 1] = NEIGHBOR[(i + 1) * 4 + j];
		}
		attribute[i] = TRIANGLE[(i + 1) * 5 + 4];
	}

	for (i = 0; i < Point_Num; i++)
	{
		Point[i * 2 + 0] = POINT[(i + 1) * 4 + 1];
		Point[i * 2 + 1] = POINT[(i + 1) * 4 + 2];
	}
}

/*compute edge length of each triangle*/
void EdgeLength(long Tri, double *Point, long *Triangle, double *Length)
{
	long int i, num1, num2, num3;
	double P1[2] = {0};
	double P2[2] = {0};
	double P3[2] = {0};
	for (i = 0; i < Tri; i++)
	{
		num1 = Triangle[i * 3 + 0] - 1;
		num2 = Triangle[i * 3 + 1] - 1;
		num3 = Triangle[i * 3 + 2] - 1;
		P1[0] = Point[num1 * 2 + 0];
		P1[1] = Point[num1 * 2 + 1];
		P2[0] = Point[num2 * 2 + 0];
		P2[1] = Point[num2 * 2 + 1];
		P3[0] = Point[num3 * 2 + 0];
		P3[1] = Point[num3 * 2 + 1];

		Length[i * 3 + 0] = sqrt(pow(P1[0] - P2[0], 2) + pow(P1[1] - P2[1], 2));
		Length[i * 3 + 1] = sqrt(pow(P2[0] - P3[0], 2) + pow(P2[1] - P3[1], 2));
		Length[i * 3 + 2] = sqrt(pow(P1[0] - P3[0], 2) + pow(P1[1] - P3[1], 2));
	}
}

/*compute minimum diameter of a circle that can be put inside a triangle, this is used to compute time step dt*/
double InDiameter_Min(long Tri, double *Point, long *Triangle, double *Length)
{
	EdgeLength(Tri, Point, Triangle, Length);
	double r, a, b, c, s, S, min;
	long i;
	double *radius;
	radius = (double *)malloc(sizeof(double) * Tri);

	for (i = 0; i < Tri; i++)
	{
		a = Length[i * 3 + 0];
		b = Length[i * 3 + 1];
		c = Length[i * 3 + 2];
		s = 0.5 * (a + b + c);
		S = sqrt(s * (s - a) * (s - b) * (s - c));

		radius[i] = 2 * S / s;
	}

	min = radius[0];
	for (i = 1; i < Tri; i++)
	{
		if (radius[i] < min)
		{
			min = radius[i];
		}
	}

	return (min);
}

/*compute time step dt*/
double Calculate_Time_Step(long Tri, double *Point, long *Triangle, double *Vp)
{
	long i, j, num;
	double *Length;
	Length = (double *)malloc(sizeof(double) * Tri * 3);
	double min, dt;
	min = 0;
	double Vmax = Vp[0];
	for (num = 1; num < Tri; num++)
	{
		if (Vp[num] > Vmax)
		{
			Vmax = Vp[num];
		}
	}
	EdgeLength(Tri, Point, Triangle, Length);
	min = InDiameter_Min(Tri, Point, Triangle, Length);
	dt = min / Vmax / (2 * N + 1) * 0.99;
	printf("***warning: dt can be adapted if problem occurs!***\n");
	return (dt);
}

/*get number of triangles inside absorbing boundary*/
long Get_Absorb_Square_Num(long Tri, double *Point, long *Triangle)
{
	long num, i, count, k, p;
	count = 0;
	double vertice[3][2];
	for (num = 0; num < Tri; num++)
	{
		for (i = 0; i < 3; i++)
		{
			k = Triangle[num * 3 + i] - 1;
			vertice[i][0] = Point[k * 2 + 0];
			vertice[i][1] = Point[k * 2 + 1];

			if (vertice[i][0] > (bound_x + model_w) || vertice[i][0] < bound_x || vertice[i][1] > (bound_z + model_h) || vertice[i][1] < bound_z)
			{
				count++;
				break;
			}
		}
	}
	printf("Absorbing Triangle Number %ld \n", count);
	return (count);
}

/*get number of triangles outside absorbing boundary, i.e., inside actual model*/
long Get_Inside_Square_Num(long Tri, double *Point, long *Triangle)
{
	long num, i, count, k, p;
	count = 0;
	double vertice[3][2];
	for (num = 0; num < Tri; num++)
	{
		for (i = 0; i < 3; i++)
		{
			k = Triangle[num * 3 + i] - 1;
			vertice[i][0] = Point[k * 2 + 0];
			vertice[i][1] = Point[k * 2 + 1];

			if (vertice[i][0] <= (bound_x + model_w) || vertice[i][0] >= bound_x || vertice[i][1] <= (bound_z + model_h) || vertice[i][1] >= bound_z)
			{
				count++;
				break;
			}
		}
	}

	return (count);
}

/*set triangles inside absorbing boundary*/
void Get_absorb_square(long Tri, double *Point, long *Triangle, long *absorbing_triangle, long *inside_triangle)
{
	long num, i, count1, count2, k, p;
	count1 = 0;
	count2 = 0;
	double vertice[3][2];
	for (num = 0; num < Tri; num++)
	{
		for (i = 0; i < 3; i++)
		{
			k = Triangle[num * 3 + i] - 1;
			vertice[i][0] = Point[k * 2 + 0];
			vertice[i][1] = Point[k * 2 + 1];

			if (vertice[i][0] > (bound_x + model_w) || vertice[i][0] < bound_x || vertice[i][1] > (bound_z + model_h) || vertice[i][1] < bound_z)
			{
				absorbing_triangle[count1] = num + 1;
				count1++;
				break;
			}

			if (vertice[i][0] <= (bound_x + model_w) || vertice[i][0] >= bound_x || vertice[i][1] <= (bound_z + model_h) || vertice[i][1] >= bound_z)
			{
				inside_triangle[count2] = num + 1;
				count2++;
				break;
			}
		}
	}
}

/*get the position for applying source */
void Get_Source(long Tri, long Point_Num, double *Point, long *Triangle, long *source, double *source_position)
{
	long num, i, count, k, source_node;
	count = 0;
	long flag = 0;
	for (i = 0; i < Point_Num; i++)
	{
		double error_x = Point[i * 2 + 0] - source_position[0];
		double error_z = Point[i * 2 + 1] - source_position[1];

		if (error_x < 0)
		{
			error_x = -error_x;
		}

		if (error_z < 0)
		{
			error_z = -error_z;
		}

		if (error_x < 1e-3 && error_z < 1e-3)
		{
			source_node = i + 1;
			flag = 1;
			printf("source_node %ld\n", source_node);
			break;
		}
	}

	if (flag == 0)
	{
		printf("\nsource not found\n");
	}

	if (flag == 1)
	{
		for (num = 0; num < Tri; num++)
		{
			for (i = 0; i < 3; i++)
			{
				k = Triangle[num * 3 + i];
				if (k == source_node)
				{
					source[count * 2 + 0] = num + 1;
					printf("count %ld  ", count);
					switch (i)
					{
					case 0:
						source[count * 2 + 1] = 6;
						break;
					case 1:
						source[count * 2 + 1] = 1;
						break;
					case 2:
						source[count * 2 + 1] = 3;
						break;
					}
					printf("source[%ld]  %ld  %ld\n", count, source[count * 2 + 0], source[count * 2 + 1]);
					count++;
					break;
				}
			}

			if (count == Source_tri)
			{
				break;
			}
		}
	}
}

/*get the position for applying receiver */
void Get_Receiver(long Tri, long Point_Num, double *Point, long *Triangle, long *receiver, double *receiver_position)
{
	long num, i, j, count, k, p, n;
	count = 0;
	long *receiver_node;
	receiver_node = (long *)malloc(sizeof(long) * Channel);
	long Num[Channel];
	for (i = 0; i < Point_Num; i++)
	{
		for (j = 0; j < Channel; j++)
		{
			double error_x = Point[i * 2 + 0] - receiver_position[j * 2 + 0];
			double error_z = Point[i * 2 + 1] - receiver_position[j * 2 + 1];

			if (error_x < 0)
			{
				error_x = -error_x;
			}

			if (error_z < 0)
			{
				error_z = -error_z;
			}

			if (error_x < 1e-2 && error_z < 1e-2)
			{
				receiver_node[j] = i + 1;
				break;
			}
		}
	}

	for (n = 0; n < Channel; n++)
	{
		for (num = 0; num < Tri; num++)
		{
			for (i = 0; i < 3; i++)
			{
				k = Triangle[num * 3 + i];
				if (k == receiver_node[n])
				{
					receiver[n * 2 + 0] = num + 1;
					switch (i)
					{
					case 0:
						receiver[n * 2 + 1] = 6;
						Num[n] = i;
						break;
					case 1:
						receiver[n * 2 + 1] = 1;
						Num[n] = i;
						break;
					case 2:
						receiver[n * 2 + 1] = 3;
						Num[n] = i;
						break;
					}
				}
			}
		}
		printf("receiver[%ld] %ld  receiver_position_x[%ld] %f receiver_position_z[%ld] %f\n", n, receiver[n * 2 + 0], n, receiver_position[n * 2 + 0], n, receiver_position[n * 2 + 1]);
	}
}

/*compute Jacobian of each triangle and edge*/
void Jacobian(long Tri, double *Jacobi_Edge, double *Jacobi, double *Point, long int *Triangle)
{
	double *xr;
	xr = (double *)malloc(sizeof(double) * Tri);
	double *yr;
	yr = (double *)malloc(sizeof(double) * Tri);
	double *xs;
	xs = (double *)malloc(sizeof(double) * Tri);
	double *ys;
	ys = (double *)malloc(sizeof(double) * Tri);
	double *Length;
	Length = (double *)malloc(sizeof(double) * Tri * 3);

	EdgeLength(Tri, Point, Triangle, Length);
	long int i, j, num1, num2, num3;
	for (i = 0; i < Tri; i++)
	{
		num1 = Triangle[i * 3 + 0] - 1;
		num2 = Triangle[i * 3 + 1] - 1;
		num3 = Triangle[i * 3 + 2] - 1;
		xr[i] = (Point[num2 * 2 + 0] - Point[num1 * 2 + 0]) * 0.5;
		yr[i] = (Point[num2 * 2 + 1] - Point[num1 * 2 + 1]) * 0.5;

		xs[i] = (Point[num3 * 2 + 0] - Point[num1 * 2 + 0]) * 0.5;
		ys[i] = (Point[num3 * 2 + 1] - Point[num1 * 2 + 1]) * 0.5;

		Jacobi[i] = xr[i] * ys[i] - xs[i] * yr[i];
		for (j = 0; j < 3; j++)
		{
			Jacobi_Edge[i * 3 + j] = Length[i * 3 + j] / 2;
		}
	}
}

/*compute angles between two edges*/
void Angle(long Tri, double *Point, long int *Triangle, double *angle)
{
	long int i;
	double *Length;
	Length = (double *)malloc(sizeof(double) * Tri * 3);
	double L1, L2, L3;
	EdgeLength(Tri, Point, Triangle, Length);

	for (i = 0; i < Tri; i++)
	{
		L1 = Length[i * 3 + 0];
		L2 = Length[i * 3 + 1];
		L3 = Length[i * 3 + 2];

		angle[i * 3 + 0] = acos((L1 * L1 + L3 * L3 - L2 * L2) / (2 * L1 * L3));
		angle[i * 3 + 1] = acos((L1 * L1 + L2 * L2 - L3 * L3) / (2 * L1 * L2));
		angle[i * 3 + 2] = acos((L2 * L2 + L3 * L3 - L1 * L1) / (2 * L2 * L3));
	}
}

/*compute normal vectors of each edge*/
void Normals(long Tri, double *Point, long int *Triangle, double *normal1, double *normal2, double *normal3)
{
	long int i, j;
	long int num1, num2, num3;
	double A[2];
	double P1B1;
	double B1[2];
	double AB1[2];
	double P2B2;
	double B2[2];
	double AB2[2];
	double P3B3;
	double B3[2];
	double AB3[2];

	double unitx[2] = {1, 0};
	double unity[2] = {0, 1};
	double P1[2] = {0};
	double P2[2] = {0};
	double P3[2] = {0};
	double *Length;
	Length = (double *)malloc(sizeof(double) * Tri * 3);
	double *angle;
	angle = (double *)malloc(sizeof(double) * Tri * 3);
	EdgeLength(Tri, Point, Triangle, Length);
	Angle(Tri, Point, Triangle, angle);

	for (i = 0; i < Tri; i++)
	{
		num1 = Triangle[i * 3 + 0] - 1;
		num2 = Triangle[i * 3 + 1] - 1;
		num3 = Triangle[i * 3 + 2] - 1;
		P1[0] = Point[num1 * 2 + 0];
		P1[1] = Point[num1 * 2 + 1];
		P2[0] = Point[num2 * 2 + 0];
		P2[1] = Point[num2 * 2 + 1];
		P3[0] = Point[num3 * 2 + 0];
		P3[1] = Point[num3 * 2 + 1];

		for (j = 0; j < 2; j++)
		{
			A[j] = (Length[i * 3 + 1] * P1[j] + Length[i * 3 + 2] * P2[j] + Length[i * 3 + 0] * P3[j]) / (Length[i * 3 + 1] + Length[i * 3 + 2] + Length[i * 3 + 0]);
		}

		P1B1 = sqrt(pow(A[0] - P1[0], 2) + pow(A[1] - P1[1], 2)) * cos(0.5 * angle[i * 3 + 0]);

		for (j = 0; j < 2; j++)
		{
			B1[j] = P1[j] + (P1B1 / Length[i * 3 + 0]) * (P2[j] - P1[j]);
			AB1[j] = B1[j] - A[j];
		}

		normal1[i * 2 + 0] = (AB1[0] * unitx[0] + AB1[1] * unitx[1]) / sqrt(pow(AB1[0], 2) + pow(AB1[1], 2));
		normal1[i * 2 + 1] = (AB1[0] * unity[0] + AB1[1] * unity[1]) / sqrt(pow(AB1[0], 2) + pow(AB1[1], 2));

		P2B2 = sqrt(pow(A[0] - P2[0], 2) + pow(A[1] - P2[1], 2)) * cos(0.5 * angle[i * 3 + 1]);

		for (j = 0; j < 2; j++)
		{
			B2[j] = P2[j] + (P2B2 / Length[i * 3 + 1]) * (P3[j] - P2[j]);
			AB2[j] = B2[j] - A[j];
		}

		normal2[i * 2 + 0] = (AB2[0] * unitx[0] + AB2[1] * unitx[1]) / sqrt(pow(AB2[0], 2) + pow(AB2[1], 2));
		normal2[i * 2 + 1] = (AB2[0] * unity[0] + AB2[1] * unity[1]) / sqrt(pow(AB2[0], 2) + pow(AB2[1], 2));

		P3B3 = sqrt(pow(A[0] - P3[0], 2) + pow(A[1] - P3[1], 2)) * cos(0.5 * angle[i * 3 + 2]);

		for (j = 0; j < 2; j++)
		{
			B3[j] = P3[j] + (P3B3 / Length[i * 3 + 2]) * (P1[j] - P3[j]);
			AB3[j] = B3[j] - A[j];
		}

		normal3[i * 2 + 0] = (AB3[0] * unitx[0] + AB3[1] * unitx[1]) / sqrt(pow(AB3[0], 2) + pow(AB3[1], 2));
		normal3[i * 2 + 1] = (AB3[0] * unity[0] + AB3[1] * unity[1]) / sqrt(pow(AB3[0], 2) + pow(AB3[1], 2));
	}
}

/*set neighboring triangles and their edge*/
void Neighbor_Edge(long Tri, long int *Neighbor, long int *neighbor, long int *Edge, long int *Triangle, long int *Edge_flag)
{
	long int p1, p2, p3;
	long int P1, P2, P3;
	long int i, j;
	long int tri;
	for (i = 0; i < Tri; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (neighbor[i * 3 + j] == -1)
			{
				neighbor[i * 3 + j] = 0;
			}
		}
	}

	for (i = 0; i < Tri; i++)
	{
		Neighbor[i * 3 + 1] = neighbor[i * 3 + 0];
		Neighbor[i * 3 + 2] = neighbor[i * 3 + 1];
		Neighbor[i * 3 + 0] = neighbor[i * 3 + 2];

		P1 = Triangle[i * 3 + 0];
		P2 = Triangle[i * 3 + 1];
		P3 = Triangle[i * 3 + 2];

		tri = Neighbor[i * 3 + 0] - 1;
		p1 = Triangle[tri * 3 + 0];
		p2 = Triangle[tri * 3 + 1];
		p3 = Triangle[tri * 3 + 2];

		if (Neighbor[i * 3 + 0])
		{
			tri = Neighbor[i * 3 + 0] - 1;
			p1 = Triangle[tri * 3 + 0];
			p2 = Triangle[tri * 3 + 1];
			p3 = Triangle[tri * 3 + 2];

			if (P1 == p1)
			{
				if (P2 == p2)
				{
					Edge[i * 3 + 0] = 1;
				}
				if (P2 == p3)
				{
					Edge[i * 3 + 0] = 3;
				}
			}

			if (P1 == p2)
			{
				if (P2 == p1)
				{
					Edge[i * 3 + 0] = 1;
				}
				if (P2 == p3)
				{
					Edge[i * 3 + 0] = 2;
				}
			}

			if (P1 == p3)
			{
				if (P2 == p1)
				{
					Edge[i * 3 + 0] = 3;
				}
				if (P2 == p2)
				{
					Edge[i * 3 + 0] = 2;
				}
			}
		}

		tri = Neighbor[i * 3 + 1] - 1;
		p1 = Triangle[tri * 3 + 0];
		p2 = Triangle[tri * 3 + 1];
		p3 = Triangle[tri * 3 + 2];
		if (Neighbor[i * 3 + 1])
		{
			if (P2 == p1)
			{
				if (P3 == p2)
				{
					Edge[i * 3 + 1] = 1;
				}
				if (P3 == p3)
				{
					Edge[i * 3 + 1] = 3;
				}
			}

			if (P2 == p2)
			{
				if (P3 == p1)
				{
					Edge[i * 3 + 1] = 1;
				}
				if (P3 == p3)
				{
					Edge[i * 3 + 1] = 2;
				}
			}

			if (P2 == p3)
			{
				if (P3 == p1)
				{
					Edge[i * 3 + 1] = 3;
				}
				if (P3 == p2)
				{
					Edge[i * 3 + 1] = 2;
				}
			}
		}

		tri = Neighbor[i * 3 + 2] - 1;
		p1 = Triangle[tri * 3 + 0];
		p2 = Triangle[tri * 3 + 1];
		p3 = Triangle[tri * 3 + 2];
		if (Neighbor[i * 3 + 2])
		{
			if (P3 == p1)
			{
				if (P1 == p2)
				{
					Edge[i * 3 + 2] = 1;
				}
				if (P1 == p3)
				{
					Edge[i * 3 + 2] = 3;
				}
			}

			if (P3 == p2)
			{
				if (P1 == p1)
				{
					Edge[i * 3 + 2] = 1;
				}
				if (P1 == p3)
				{
					Edge[i * 3 + 2] = 2;
				}
			}

			if (P3 == p3)
			{
				if (P1 == p1)
				{
					Edge[i * 3 + 2] = 3;
				}
				if (P1 == p2)
				{
					Edge[i * 3 + 2] = 2;
				}
			}
		}
	}

	for (i = 0; i < Tri; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (Neighbor[i * 3 + j] == 0)
			{
				Neighbor[i * 3 + j] = 1;
			}

			if (Edge[i * 3 + j])
			{
				Edge_flag[i * 3 + j] = Edge[i * 3 + j] / Edge[i * 3 + j];
			}
			else
			{
				Edge_flag[i * 3 + j] = 0;
			}
		}
	}
}

/*mark if neighboring triangle is in the other partition area*/
void Neighbor_Partition_Mark(long Tri, long int *Neighbor, long int *epart, long int *Neighbor_Mark)
{
	// calculate neighboring partition mark
	long j, n, num;
	for (num = 0; num < Tri; num++)
	{
		for (j = 0; j < 3; j++)
		{
			n = Neighbor[num * 3 + j] - 1;

			if (n == -1)
				Neighbor_Mark[num * 3 + j] = -1;
			else
				Neighbor_Mark[num * 3 + j] = epart[n];
		}
	}
}

/*count boundary triangle (between partitions) numbers*/
long Calculate_Boundary_Num(long Tri, long int *Neighbor_Mark)
{
	long j, n, num, mark, n0, n1, n2;

	long Boundary_Num = 0;
	for (num = 0; num < Tri; num++)
	{
		n0 = Neighbor_Mark[num * 3 + 0];
		n1 = Neighbor_Mark[num * 3 + 1];
		n2 = Neighbor_Mark[num * 3 + 2];

		if ((n0 != n1 || n0 != n2 || n1 != n2) && (n0 != -1 && n1 != -1 && n2 != -1))
		{
			Boundary_Num = Boundary_Num + 1;
		}

		if (n1 != n2 && n0 == -1 && n1 != -1 && n2 != -1)
		{
			Boundary_Num = Boundary_Num + 1;
		}

		if (n0 != n2 && n1 == -1 && n0 != -1 && n2 != -1)
		{
			Boundary_Num = Boundary_Num + 1;
		}

		if (n0 != n1 && n2 == -1 && n0 != -1 && n1 != -1)
		{
			Boundary_Num = Boundary_Num + 1;
		}
	}

	return (Boundary_Num);
}

/*set neighboring triangles in the boundary between partitions*/
void Partition_Boundary_Neighbor(long Tri, long int *Neighbor, long int *epart, long int *Neighbor_Mark, long *Series, long *Series_order, long *Bound_Series, long *Triangle_part, long *Boundary_Triangle)
{
	long j, m, n, num, i, k, order, size, sum;
	m = 0;
	long order_sum = 0;
	long size_sum = 0;
	for (num = 0; num < Process_Num; num++)
	{
		order = 0;
		size = 0;
		for (j = 0; j < Tri; j++)
		{
			if (epart[j] == num)
			{
				Triangle_part[order_sum + order] = j + 1;

				order++;

				for (i = 0; i < 3; i++)
				{
					if (Neighbor_Mark[j * 3 + i] != num && Neighbor_Mark[j * 3 + i] != -1)
					{
						Boundary_Triangle[m] = j + 1;
						m++;
						size++;
						break;
					}
				}
			}
		}

		Bound_Series[num] = size; // 进程中边界的三角形数量 number of triangles inside  boundary between different partitions in one process
		Series[num] = order;	  // 进程中三角形的数量 number of triangles inside partitions in one process

		Series_order[num] = order_sum;
		order_sum += order;
	}
	printf("m is %ld\n", m);
}


/*set all node for saving to a file, which can be used for visualizing wavefield snapshots*/
void set_All_node(long Tri, double *All_node, double *Point, long *Triangle, double FeketeNew[][2])
{
	long int i, num, m, p;
	double r, s;
	double vertice[3][2];
	for (num = 0; num < Tri; num++)
	{

		for (m = 0; m < 3; m++)
		{
			p = Triangle[num * 3 + m] - 1;
			vertice[m][0] = Point[p * 2 + 0];
			vertice[m][1] = Point[p * 2 + 1];
		}
		for (i = 0; i < Np; i++)
		{
			r = FeketeNew[i][0];
			s = FeketeNew[i][1];
			All_node[num * Np * 2 + i * 2 + 0] = -(r + s) / 2 * vertice[0][0] + (r + 1) / 2 * vertice[1][0] + (s + 1) / 2 * vertice[2][0];
			All_node[num * Np * 2 + i * 2 + 1] = -(r + s) / 2 * vertice[0][1] + (r + 1) / 2 * vertice[1][1] + (s + 1) / 2 * vertice[2][1];
		}
	}
}

/*coordinates transform, from local triangle to the model coordinates*/
void Rstoxy(double position[][2], double vertice[][2], double FeketeNew[][2])
{
	long int i;
	double r, s;

	for (i = 0; i < Np; i++)
	{
		r = FeketeNew[i][0];
		s = FeketeNew[i][1];
		position[i][0] = -(r + s) / 2 * vertice[0][0] + (r + 1) / 2 * vertice[1][0] + (s + 1) / 2 * vertice[2][0];
		position[i][1] = -(r + s) / 2 * vertice[0][1] + (r + 1) / 2 * vertice[1][1] + (s + 1) / 2 * vertice[2][1];
	}
}

/*set nodes inside absorbing bounary*/
void Absorb_Tri_Node(long Absorb_Num, double *Absorb_tri_nodes, long int *Absorb_tri, double *Point, long int *Triangle, double FeketeNew[][2])
{
	long int i, j, k, num;
	double position[Np][2];
	double vertice[3][2];

	for (num = 0; num < Absorb_Num; num++)
	{
		for (i = 0; i < 3; i++)
		{
			k = Absorb_tri[num] - 1;
			j = Triangle[k * 3 + i] - 1;
			vertice[i][0] = Point[j * 2 + 0];
			vertice[i][1] = Point[j * 2 + 1];
		}

		Rstoxy(position, vertice, FeketeNew);
		for (i = 0; i < Np; i++)
		{
			Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] = position[i][0];
			Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] = position[i][1];
		}
	}
}


/*compute absorbing functions*/
void Absorbing_Fun(long Tri, long Absorb_Num, double beta, double *absorb_fun, long int *Absorb_tri, double *point, long int *triangle, double FeketeNew[][2])
{
	long int i, j, num, k, order;

	double *Absorb_tri_nodes;
	Absorb_tri_nodes = (double *)malloc(sizeof(double) * Absorb_Num * Np * 2);
	Absorb_Tri_Node(Absorb_Num, Absorb_tri_nodes, Absorb_tri, point, triangle, FeketeNew);
	double position[Np][2];

	double distance, dx, dy;
	double Lx = bound_x;
	double Ly = bound_z;

	for (num = 0; num < Tri * Np; num++)
	{
		absorb_fun[num] = 1;
	}

	for (num = 0; num < Absorb_Num; num++)
	{
		for (i = 0; i < Np; i++)
		{
			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] < bound_x && (Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] > bound_z && Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] < bound_z + model_h))
			{
				k = Absorb_tri[num] - 1;
				distance = bound_x - Absorb_tri_nodes[num * Np * 2 + i * 2 + 0];
				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Lx);
			}

			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] > bound_x + model_w && (Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] > bound_z && Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] < bound_z + model_h))
			{
				k = Absorb_tri[num] - 1;
				distance = Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] - bound_x - model_w;
				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Lx);
			}

			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] < bound_z && (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] > bound_x && Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] < bound_x + model_w))
			{
				k = Absorb_tri[num] - 1;
				distance = bound_z - Absorb_tri_nodes[num * Np * 2 + i * 2 + 1];
				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Ly);
			}

			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] > bound_z + model_h && (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] > bound_x && Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] < bound_x + model_w))
			{
				k = Absorb_tri[num] - 1;
				distance = Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] - bound_z - model_h;
				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Ly);
			}

			dx = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] - bound_x);
			dy = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] - bound_z);
			distance = sqrt(dx * dx + dy * dy);
			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] < bound_x && Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] < bound_z && distance < (bound_x))
			{
				k = Absorb_tri[num] - 1;

				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Lx);
			}

			dx = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] - (bound_x + model_w));
			dy = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] - bound_z);
			distance = sqrt(dx * dx + dy * dy);
			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] > bound_x + model_w && Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] < bound_z && distance < (bound_x))
			{
				k = Absorb_tri[num] - 1;

				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Lx);
			}

			dx = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] - (bound_x + model_w));
			dy = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] - (bound_z + model_h));
			distance = sqrt(dx * dx + dy * dy);

			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] > bound_x + model_w && Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] > bound_z + model_h && distance < (bound_x))
			{
				k = Absorb_tri[num] - 1;

				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Lx);
			}
			dx = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] - bound_x);
			dy = fabs(Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] - (bound_z + model_h));
			distance = sqrt(dx * dx + dy * dy);

			if (Absorb_tri_nodes[num * Np * 2 + i * 2 + 0] < bound_x && Absorb_tri_nodes[num * Np * 2 + i * 2 + 1] > bound_z + model_h && distance < (bound_x))
			{
				k = Absorb_tri[num] - 1;

				order = k * Np + i;
				absorb_fun[order] = 0.5 + 0.5 * cos(beta * distance * Pi / Lx);
			}
		}
	}
}

/*set IBE node ordering by DG theory*/
void IBE(long int *ibe)
{
	if (N == 4)
	{
		ibe[0 * (N + 1) + 0] = 8;
		ibe[0 * (N + 1) + 1] = 7;
		ibe[0 * (N + 1) + 2] = 6;
		ibe[0 * (N + 1) + 3] = 5;
		ibe[0 * (N + 1) + 4] = 4;
		ibe[1 * (N + 1) + 0] = 4;
		ibe[1 * (N + 1) + 1] = 15;
		ibe[1 * (N + 1) + 2] = 14;
		ibe[1 * (N + 1) + 3] = 13;
		ibe[1 * (N + 1) + 4] = 9;
		ibe[2 * (N + 1) + 0] = 8;
		ibe[2 * (N + 1) + 1] = 12;
		ibe[2 * (N + 1) + 2] = 11;
		ibe[2 * (N + 1) + 3] = 10;
		ibe[2 * (N + 1) + 4] = 9;
	}

	if (N == 3)
	{
		ibe[0 * (N + 1) + 0] = 10;
		ibe[0 * (N + 1) + 1] = 8;
		ibe[0 * (N + 1) + 2] = 5;
		ibe[0 * (N + 1) + 3] = 1;
		ibe[1 * (N + 1) + 0] = 1;
		ibe[1 * (N + 1) + 1] = 2;
		ibe[1 * (N + 1) + 2] = 3;
		ibe[1 * (N + 1) + 3] = 4;
		ibe[2 * (N + 1) + 0] = 10;
		ibe[2 * (N + 1) + 1] = 9;
		ibe[2 * (N + 1) + 2] = 7;
		ibe[2 * (N + 1) + 3] = 4;
	}

	if (N == 2)
	{
		ibe[0 * (N + 1) + 0] = 6;
		ibe[0 * (N + 1) + 1] = 4;
		ibe[0 * (N + 1) + 2] = 1;
		ibe[1 * (N + 1) + 0] = 1;
		ibe[1 * (N + 1) + 1] = 2;
		ibe[1 * (N + 1) + 2] = 3;
		ibe[2 * (N + 1) + 0] = 6;
		ibe[2 * (N + 1) + 1] = 5;
		ibe[2 * (N + 1) + 2] = 3;
	}

	if (N == 1)
	{
		ibe[0 * (N + 1) + 0] = 3;
		ibe[0 * (N + 1) + 1] = 1;
		ibe[1 * (N + 1) + 0] = 1;
		ibe[1 * (N + 1) + 1] = 2;
		ibe[2 * (N + 1) + 0] = 3;
		ibe[2 * (N + 1) + 1] = 2;
	}
}

/*compute A and B matrix by DG theory*/
void Elastic_A_B(long Tri, double *A, double *B, double *Lambda, double *Miu, double *Den)
{
	long int num;
	for (num = 0; num < Tri; num++)
	{
		A[num * 5 + 0] = 1 / Den[num];
		A[num * 5 + 1] = 1 / Den[num];
		A[num * 5 + 2] = Lambda[num] + 2 * Miu[num];
		A[num * 5 + 3] = Lambda[num];
		A[num * 5 + 4] = Miu[num];

		B[num * 5 + 0] = 1 / Den[num];
		B[num * 5 + 1] = 1 / Den[num];
		B[num * 5 + 2] = Lambda[num];
		B[num * 5 + 3] = Lambda[num] + 2 * Miu[num];
		B[num * 5 + 4] = Miu[num];
	}
}

/*compute central flux by DG theory*/
void Central_flux(long Tri, double *Jacobi, double *Jacobi_Edge, double *Q, double *A, double *B, double *normal1, double *normal2, double *normal3)
{
	long int num, j, m, n;
	double nx, nz;

	for (num = 0; num < Tri; num++)
	{
		for (j = 0; j < 3; j++)
		{
			if (j == 0)
			{
				nx = normal1[num * 2 + 0];
				nz = normal1[num * 2 + 1];
			}

			if (j == 1)
			{
				nx = normal2[num * 2 + 0];
				nz = normal2[num * 2 + 1];
			}

			if (j == 2)
			{
				nx = normal3[num * 2 + 0];
				nz = normal3[num * 2 + 1];
			}

			Q[num * 3 * 5 * 2 + j * 5 * 2 + 0 * 2 + 0] = 0.5 * nx * A[num * 5 + 0] * Jacobi_Edge[num * 3 + j] / Jacobi[num];
			Q[num * 3 * 5 * 2 + j * 5 * 2 + 0 * 2 + 1] = 0.5 * nz * B[num * 5 + 0] * Jacobi_Edge[num * 3 + j] / Jacobi[num];

			Q[num * 3 * 5 * 2 + j * 5 * 2 + 1 * 2 + 0] = 0.5 * nz * A[num * 5 + 1] * Jacobi_Edge[num * 3 + j] / Jacobi[num];
			Q[num * 3 * 5 * 2 + j * 5 * 2 + 1 * 2 + 1] = 0.5 * nx * B[num * 5 + 1] * Jacobi_Edge[num * 3 + j] / Jacobi[num];

			Q[num * 3 * 5 * 2 + j * 5 * 2 + 2 * 2 + 0] = 0.5 * nx * A[num * 5 + 2] * Jacobi_Edge[num * 3 + j] / Jacobi[num];
			Q[num * 3 * 5 * 2 + j * 5 * 2 + 2 * 2 + 1] = 0.5 * nz * B[num * 5 + 2] * Jacobi_Edge[num * 3 + j] / Jacobi[num];

			Q[num * 3 * 5 * 2 + j * 5 * 2 + 3 * 2 + 0] = 0.5 * nx * A[num * 5 + 3] * Jacobi_Edge[num * 3 + j] / Jacobi[num];
			Q[num * 3 * 5 * 2 + j * 5 * 2 + 3 * 2 + 1] = 0.5 * nz * B[num * 5 + 3] * Jacobi_Edge[num * 3 + j] / Jacobi[num];

			Q[num * 3 * 5 * 2 + j * 5 * 2 + 4 * 2 + 0] = 0.5 * nz * A[num * 5 + 4] * Jacobi_Edge[num * 3 + j] / Jacobi[num];
			Q[num * 3 * 5 * 2 + j * 5 * 2 + 4 * 2 + 1] = 0.5 * nx * B[num * 5 + 4] * Jacobi_Edge[num * 3 + j] / Jacobi[num];
		}
	}
}

/*compute Dx and Dz matrix by DG theory*/
void Dx_Dz(long Tri, double *Dx, double *Dz, double Vr[][Np], double Vs[][Np], double inverV[][Np], long int *Triangle, double *Point)
{
	long int num, i, j, m, n, p, k, q;
	double X[3];
	double Z[3];
	double Xr, Zr, Xs, Zs;
	double Rx, Rz, Sx, Sz;
	double J, f1, f2;
	double Dxd[Np][Np];
	double Dzd[Np][Np];

	for (num = 0; num < Tri; num++)
	{
		for (j = 0; j < 3; j++)
		{
			m = Triangle[num * 3 + j] - 1;
			X[j] = Point[m * 2 + 0];
			Z[j] = Point[m * 2 + 1];
		}

		Xr = (X[1] - X[0]) * 0.5;
		Zr = (Z[1] - Z[0]) * 0.5;
		Xs = (X[2] - X[0]) * 0.5;
		Zs = (Z[2] - Z[0]) * 0.5;
		J = Xr * Zs - Xs * Zr;
		Rx = Zs / J;
		Rz = -Xs / J;
		Sx = -Zr / J;
		Sz = Xr / J;

		for (i = 0; i < Np; i++)
		{
			for (n = 0; n < Np; n++)
			{
				Dxd[i][n] = Rx * Vr[i][n] + Sx * Vs[i][n];
				Dzd[i][n] = Rz * Vr[i][n] + Sz * Vs[i][n];
			}
		}

		for (p = 0; p < Np; p++)
		{
			for (q = 0; q < Np; q++)
			{
				f1 = 0;
				f2 = 0;
				for (k = 0; k < Np; k++)
				{
					f1 = f1 + Dxd[p][k] * inverV[k][q];
					f2 = f2 + Dzd[p][k] * inverV[k][q];
				}
				Dx[num * Np * Np + p * Np + q] = f1;
				Dz[num * Np * Np + p * Np + q] = f2;
			}
		}
	}
}

/*set nodes shared by neighbors*/
void NeighborNodes(long Tri, long int *Node_Neighbor, long int *Triangle, long int *Edge, long int *Neighbor, long int *ibe)
{
	long int P[3];
	long int order[N + 1];
	long int num, i, j, p, edge, neighbor;

	for (num = 0; num < Tri; num++)
	{
		P[0] = Triangle[num * 3 + 0];
		P[1] = Triangle[num * 3 + 1];
		P[2] = Triangle[num * 3 + 2];

		for (j = 0; j < 3; j++)
		{
			edge = Edge[num * 3 + j] - 1;
			neighbor = Neighbor[num * 3 + j] - 1;
			if (edge != -1)
			{
				if (j == 2)
				{
					if (edge == 2)
						p = Triangle[neighbor * 3 + 0];
					else
						p = Triangle[neighbor * 3 + edge];

					for (i = 0; i < N + 1; i++)
					{
						order[i] = ibe[edge * (N + 1) + i];
					}

					if (P[0] == p)
					{
						for (i = 0; i < N + 1; i++)
						{
							Node_Neighbor[num * 3 * (N + 1) + j * (N + 1) + i] = order[i];
						}
					}
					else
					{
						for (i = 0; i < N + 1; i++)
						{
							Node_Neighbor[num * 3 * (N + 1) + j * (N + 1) + i] = order[N - i];
						}
					}
				}
				else
				{
					if (edge == 2)
						p = Triangle[neighbor * 3 + 0];
					else
						p = Triangle[neighbor * 3 + edge];

					for (i = 0; i < N + 1; i++)
					{
						order[i] = ibe[edge * (N + 1) + i];
					}
					if (P[j] == p)
					{
						for (i = 0; i < N + 1; i++)
						{
							Node_Neighbor[num * 3 * (N + 1) + j * (N + 1) + i] = order[i];
						}
					}
					else
					{
						for (i = 0; i < N + 1; i++)
						{
							Node_Neighbor[num * 3 * (N + 1) + j * (N + 1) + i] = order[N - i];
						}
					}
				}
			}
			else
			{
				for (i = 0; i < N + 1; i++)
				{
					Node_Neighbor[num * 3 * (N + 1) + j * (N + 1) + i] = 1;
				}
			}
		}
	}
}

/*compute Ricker wavelet*/
double Ricker_Source(long int time, double dt, double t0)
{
	double s;
	double t;
	t = dt * time;
	s = (1 - 2 * Pi * Pi * Fp * Fp * (t - t0) * (t - t0)) * exp(-Pi * Pi * Fp * Fp * (t - t0) * (t - t0));
	return (s);
}

/*compute time delay of Ricker wavelet*/
double Ricker_Delay_Time(double dt)
{
	double t1, t2, s1, s2;
	long T = 0;
	int flag = 0;
	while (flag < 1)
	{
		t1 = T * dt;
		s1 = (1 - 2 * Pi * Pi * Fp * Fp * (t1) * (t1)) * exp(-Pi * Pi * Fp * Fp * (t1) * (t1));
		t2 = (T + 1) * dt;
		s2 = (1 - 2 * Pi * Pi * Fp * Fp * (t2) * (t2)) * exp(-Pi * Pi * Fp * Fp * (t2) * (t2));
		if (s2 > s1 && s2 < 1e-8)
		{
			flag++;
		}
		else
		{
			T++;
		}
	}
	return (t2 * 2.6);
}