/*this file contains functions used for the MPI code*/

/*prototyping functions*/

/*apply absorbing boundary*/
void Absorb_Boundary(int rank, long Tri, double *Vx, double *Vz, double *Txx, double *Tzz, double *Txz,
					 double *absorb_fun_d, long *Triangle_part, long *Series, long *Series_order);

void source_kernel(double *Vz, long *source, double *fun, double w, long m);
void Record(double *Vx, double *Vz, long time, double *record_x_d, double *record_z_d, long *receiver_d);

void LU_Velocity_kernel(int rank, long Tri, double Dt, double *Vx1, double *Vz1, double *Txx1, double *Tzz1, double *Txz1,
						long *Triangle_part, long *Neighbor_Mark, long *Series, long *Series_order,
						double *A, double *B, double *Dx, double *Dz, double *Q, double *S, long int *ibe,
						long int *nodes_neighbor, long int *Neighbor, long int *Edge_flag);

void LU_Stress_Elastic_kernel(int rank, long Tri, double Dt, double *Txx1, double *Tzz1, double *Txz1, double *Vx1, double *Vz1,
							  long *Triangle_part, long *Neighbor_Mark, long *Series, long *Series_order,
							  double *A, double *B, double *Dx, double *Dz, double *Q, double *S, long int *ibe,
							  long int *nodes_neighbor, long int *Neighbor, long int *Edge_flag);

/**************************************************************************************************/

/*apply absorbing boundary*/
void Absorb_Boundary(int rank, long Tri, double *Vx, double *Vz, double *Txx, double *Tzz, double *Txz, double *absorb_fun_d, long *Triangle_part, long *Series, long *Series_order)
{
	long m, k, i, id, num;
	m = rank;
	{
#pragma omp parallel for private(num, id, i)
		for (k = Series_order[m]; k < (Series_order[m] + Series[m]); k++)
		{
			num = Triangle_part[k] - 1;

			for (i = 0; i < Np; i++)
			{
				id = num * Np + i;
				Vx[id] *= absorb_fun_d[id];
				Vz[id] *= absorb_fun_d[id];
				Txx[id] *= absorb_fun_d[id];
				Tzz[id] *= absorb_fun_d[id];
				Txz[id] *= absorb_fun_d[id];
			}
		}
	}
}

/*apply souce to the designated triangle*/
void source_kernel(double *Vz, long *source, double *fun, double w, long m)
{
	int i;
	long order;
	for (i = 0; i < Np; i++)
	{
		order = Np * (source[m * 2 + 0] - 1) + i;
		Vz[order] += fun[m * Np + i] * w;
	}
}

/*save Vx and Vz record at a given position*/
void Record(double *Vx, double *Vz, long time, double *record_x_d, double *record_z_d, long *receiver_d)
{
	long i, k, receiver_tri, n;

	for (k = 0; k < Channel; k++)
	{
		receiver_tri = receiver_d[k * 3 + 0];
		n = receiver_d[k * 3 + 1];
		record_x_d[time * Channel + k] = Vx[Np * (receiver_tri - 1) + n - 1];
		record_z_d[time * Channel + k] = Vz[Np * (receiver_tri - 1) + n - 1];
	}
}

/*use stress to update Vx and Vz*/
void LU_Velocity_kernel(int rank, long Tri, double Dt, double *Vx1, double *Vz1, double *Txx1, double *Tzz1, double *Txz1,
						long *Triangle_part, long *Neighbor_Mark, long *Series, long *Series_order,
						double *A, double *B, double *Dx, double *Dz, double *Q, double *S, long int *ibe, long int *nodes_neighbor, long int *Neighbor, long int *Edge_flag)
{
	long i, j, k, m, n, neighbor, num;
	long int order;
	long int order1;
	long int order2;

	double sum1, sum2, sum3, sum4;

	m = rank;
	{
#pragma omp parallel for private(num, i, j, n, sum1, sum2, sum3, sum4, order, order1, order2, neighbor)
		for (k = Series_order[m]; k < (Series_order[m] + Series[m]); k++)
		{
			long ibe_out[N + 1] = {0};
			double Dx_Txx[Np] = {0};
			double Dx_Txz[Np] = {0};
			double Dz_Tzz[Np] = {0};
			double Dz_Txz[Np] = {0};
			double A_B_Vx[Np] = {0};
			double A_B_Vz[Np] = {0};
			double dTxx[N + 1] = {0};
			double dTzz[N + 1] = {0};
			double dTxz[N + 1] = {0};
			double S_dTxx[Np] = {0};
			double S_dTzz[Np] = {0};
			double S_dTxz[Np] = {0};
			double Q_S_dVx[Np][3] = {{0}};
			double Q_S_dVz[Np][3] = {{0}};
			num = Triangle_part[k] - 1;

			for (i = 0; i < Np; i++)
			{
				sum1 = 0;
				sum2 = 0;
				sum3 = 0;
				sum4 = 0;
				for (n = 0; n < Np; n++)
				{
					order = Np * num + n;

					sum1 += Dx[num * Np * Np + i * Np + n] * Txx1[order];
					sum2 += Dx[num * Np * Np + i * Np + n] * Txz1[order];
					sum3 += Dz[num * Np * Np + i * Np + n] * Tzz1[order];
					sum4 += Dz[num * Np * Np + i * Np + n] * Txz1[order];
				}

				Dx_Txx[i] = sum1;
				Dx_Txz[i] = sum2;
				Dz_Tzz[i] = sum3;
				Dz_Txz[i] = sum4;
			}

			for (i = 0; i < Np; i++)
			{
				A_B_Vx[i] = A[num * 5 + 0] * Dx_Txx[i] + B[num * 5 + 0] * Dz_Txz[i];
				A_B_Vz[i] = A[num * 5 + 1] * Dx_Txz[i] + B[num * 5 + 1] * Dz_Tzz[i];
			}

			for (j = 0; j < 3; j++)
			{
				neighbor = Neighbor[num * 3 + j] - 1;
				for (i = 0; i < N + 1; i++)
				{
					ibe_out[i] = nodes_neighbor[num * 3 * (N + 1) + j * (N + 1) + i] - 1;
					order1 = Np * neighbor + ibe_out[i];
					order2 = Np * num + ibe[j * (N + 1) + i] - 1;

					dTxx[i] = Txx1[order1] * Edge_flag[num * 3 + j] - Txx1[order2];
					dTzz[i] = Tzz1[order1] * Edge_flag[num * 3 + j] - Tzz1[order2];
					dTxz[i] = Txz1[order1] * Edge_flag[num * 3 + j] - Txz1[order2];
				}

				// S*(Uke -Ue)
				for (i = 0; i < Np; i++)
				{
					sum1 = 0;
					sum2 = 0;
					sum3 = 0;

					for (n = 0; n < N + 1; n++)
					{
						sum1 += S[j * (N + 1) * Np + i * (N + 1) + n] * dTxx[n];
						sum2 += S[j * (N + 1) * Np + i * (N + 1) + n] * dTzz[n];
						sum3 += S[j * (N + 1) * Np + i * (N + 1) + n] * dTxz[n];
					}

					S_dTxx[i] = sum1;
					S_dTzz[i] = sum2;
					S_dTxz[i] = sum3;
				}

				// Q*S*(Uke-Ue)
				for (i = 0; i < Np; i++)
				{
					Q_S_dVx[i][j] = Q[num * 3 * 5 * 2 + j * 5 * 2 + 0 * 2 + 0] * S_dTxx[i] + Q[num * 3 * 5 * 2 + j * 5 * 2 + 0 * 2 + 1] * S_dTxz[i];
					Q_S_dVz[i][j] = Q[num * 3 * 5 * 2 + j * 5 * 2 + 1 * 2 + 0] * S_dTzz[i] + Q[num * 3 * 5 * 2 + j * 5 * 2 + 1 * 2 + 1] * S_dTxz[i];
				}
			}

			for (i = 0; i < Np; i++)
			{
				order = Np * num + i;
				Vx1[order] = (A_B_Vx[i] + Q_S_dVx[i][0] + Q_S_dVx[i][1] + Q_S_dVx[i][2]) * Dt + Vx1[order];
				Vz1[order] = (A_B_Vz[i] + Q_S_dVz[i][0] + Q_S_dVz[i][1] + Q_S_dVz[i][2]) * Dt + Vz1[order];
			}
		}
	}
}

/*use Vx and Vz to update stress*/
void LU_Stress_Elastic_kernel(int rank, long Tri, double Dt, double *Txx1, double *Tzz1, double *Txz1, double *Vx1, double *Vz1,
							  long *Triangle_part, long *Neighbor_Mark, long *Series, long *Series_order,
							  double *A, double *B, double *Dx, double *Dz, double *Q, double *S, long int *ibe, long int *nodes_neighbor, long int *Neighbor, long int *Edge_flag)
{
	long int i, j, k, m, n, neighbor;
	long int order;
	long int order1;
	long int order2;

	double sum1, sum2, sum3, sum4;
	double flag = 0;
	double flag1 = 0;
	double flag2 = 0;
	double flag3 = 0;

	long num;
	//	for(m = 0; m < Process_Num; m++)
	m = rank;
	{
#pragma omp parallel for private(num, i, j, n, sum1, sum2, sum3, sum4, order, order1, order2, neighbor)
		for (k = Series_order[m]; k < (Series_order[m] + Series[m]); k++)
		{
			double Dx_Vx[Np] = {0};
			double Dx_Vz[Np] = {0};
			double Dz_Vx[Np] = {0};
			double Dz_Vz[Np] = {0};
			double A_B_Txx[Np] = {0};
			double A_B_Tzz[Np] = {0};
			double A_B_Txz[Np] = {0};
			double dVx[N + 1] = {0};
			double dVz[N + 1] = {0};
			double S_dVx[Np] = {0};
			double S_dVz[Np] = {0};
			double Q_S_dTxx[Np][3] = {{0}};
			double Q_S_dTzz[Np][3] = {{0}};
			double Q_S_dTxz[Np][3] = {{0}};
			double Q_S_dE[Np][3] = {{0}};
			long int ibe_out[N + 1] = {0};

			num = Triangle_part[k] - 1;
			for (i = 0; i < Np; i++)
			{
				sum1 = 0;
				sum2 = 0;
				sum3 = 0;
				sum4 = 0;

				for (n = 0; n < Np; n++)
				{
					order = Np * num + n;
					sum1 += Dx[num * Np * Np + i * Np + n] * Vx1[order];
					sum2 += Dx[num * Np * Np + i * Np + n] * Vz1[order];
					sum3 += Dz[num * Np * Np + i * Np + n] * Vx1[order];
					sum4 += Dz[num * Np * Np + i * Np + n] * Vz1[order];
				}

				Dx_Vx[i] = sum1;
				Dz_Vx[i] = sum3;
				Dx_Vz[i] = sum2;
				Dz_Vz[i] = sum4;
			}

			for (i = 0; i < Np; i++)
			{
				A_B_Txx[i] = A[num * 5 + 2] * Dx_Vx[i] + B[num * 5 + 2] * Dz_Vz[i];
				A_B_Tzz[i] = A[num * 5 + 3] * Dx_Vx[i] + B[num * 5 + 3] * Dz_Vz[i];
				A_B_Txz[i] = A[num * 5 + 4] * Dx_Vz[i] + B[num * 5 + 4] * Dz_Vx[i];
			}

			for (j = 0; j < 3; j++)
			{
				neighbor = Neighbor[num * 3 + j] - 1;
				for (i = 0; i < N + 1; i++)
				{
					ibe_out[i] = nodes_neighbor[num * 3 * (N + 1) + j * (N + 1) + i] - 1;
					order1 = Np * neighbor + ibe_out[i];
					order2 = Np * num + ibe[j * (N + 1) + i] - 1;
					dVx[i] = Vx1[order1] * Edge_flag[num * 3 + j] - Vx1[order2];
					dVz[i] = Vz1[order1] * Edge_flag[num * 3 + j] - Vz1[order2];
				}

				// S*(Uke -Ue)
				for (i = 0; i < Np; i++)
				{
					sum1 = 0;
					sum2 = 0;

					for (n = 0; n < N + 1; n++)
					{
						sum1 += S[j * (N + 1) * Np + i * (N + 1) + n] * dVx[n];
						sum2 += S[j * (N + 1) * Np + i * (N + 1) + n] * dVz[n];
					}

					S_dVx[i] = sum1;
					S_dVz[i] = sum2;
				}

				// Q*S*(Uke-Ue)
				for (i = 0; i < Np; i++)
				{
					Q_S_dTxx[i][j] = Q[num * 3 * 5 * 2 + j * 5 * 2 + 2 * 2 + 0] * S_dVx[i] + Q[num * 3 * 5 * 2 + j * 5 * 2 + 2 * 2 + 1] * S_dVz[i];
					Q_S_dTzz[i][j] = Q[num * 3 * 5 * 2 + j * 5 * 2 + 3 * 2 + 0] * S_dVx[i] + Q[num * 3 * 5 * 2 + j * 5 * 2 + 3 * 2 + 1] * S_dVz[i];
					Q_S_dTxz[i][j] = Q[num * 3 * 5 * 2 + j * 5 * 2 + 4 * 2 + 0] * S_dVx[i] + Q[num * 3 * 5 * 2 + j * 5 * 2 + 4 * 2 + 1] * S_dVz[i];
				}
			}

			for (i = 0; i < Np; i++)
			{
				order = Np * num + i;

				Txx1[order] = (A_B_Txx[i] + Q_S_dTxx[i][0] + Q_S_dTxx[i][1] + Q_S_dTxx[i][2]) * Dt + Txx1[order];

				Tzz1[order] = (A_B_Tzz[i] + Q_S_dTzz[i][0] + Q_S_dTzz[i][1] + Q_S_dTzz[i][2]) * Dt + Tzz1[order];

				Txz1[order] = (A_B_Txz[i] + Q_S_dTxz[i][0] + Q_S_dTxz[i][1] + Q_S_dTxz[i][2]) * Dt + Txz1[order];
			}
		}
	}
}
