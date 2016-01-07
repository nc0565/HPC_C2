/*
** code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the bhatnagar-gross-krook collision step.
**
** the 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** a 2d grid:
**
**           cols
**       --- --- ---
**      | d | e | f |
** rows  --- --- ---
**      | a | b | c |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1d array:
**
**  --- --- --- --- --- ---
** | a | b | c | d | e | f |
**  --- --- --- --- --- ---
**
** grid indicies are:
**
**          ny
**          ^       cols(jj)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(ii) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   d2q9-bgk.exe input.params obstacles.dat
**
** be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
// #include <stddef.h>

#include "mpi.h"
#include "lbm.h"
#define MASTER 0

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int com_size;               // The number of ranks in the comunicator
	MPI_Status status;          // Struct for recieve

	param_t params /*= {0,0,0,0,0.0,0.0,0.0}*/;            /* struct to hold parameter values */
	int    ii;                  /*  generic counter */
	struct timeval timstr;      /* structure to hold elapsed time */
	struct rusage ru;           /* structure to hold CPU time--system and user */
	double tic,toc;             /* doubleing point numbers to calculate elapsed wallclock time */
	double usrtim;              /* doubleing point number to record elapsed user CPU time */
	double systim;              /* doubleing point number to record elapsed system CPU time */
	char * final_state_file = NULL;
	char * av_vels_file = NULL;
	char * param_file = NULL;

	accel_area_t accel_area;

	speed_t* cells     = NULL;    /* grid containing fluid densities */
	speed_t* tmp_cells = NULL;    /* scratch space */
	int*     obstacles = NULL;    /* grid indicating which cells are blocked */
	double*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */

	// Check world
	MPI_Comm_size(MPI_COMM_WORLD, &com_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &params.my_rank);

	// Initialise
	if (params.my_rank == MASTER)
	{
		parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);
		initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels);         
	}

	// Create mpi type for param
	int blocklengths[2] = {8,3};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Aint displacements[2];
	MPI_Get_address(&params.nx, &displacements[0]);
	MPI_Get_address(&params.density, &displacements[1]);
	displacements[1] -= displacements[0];
	displacements[0] = 0;

	MPI_Datatype mpi_param;
	MPI_Type_create_struct(2, blocklengths, displacements, types, &mpi_param);
	MPI_Type_commit(&mpi_param);

	// Broadcast params
	MPI_Bcast(&params, 1, mpi_param, MASTER, MPI_COMM_WORLD);

	// Restore rank numbers
	MPI_Comm_rank(MPI_COMM_WORLD, &params.my_rank);

	MPI_Datatype mpi_accel_area;
	MPI_Type_create_struct(1, (int[1]){2}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_INT}, &mpi_accel_area);
	MPI_Type_commit(&mpi_accel_area);

	// Broadcast accel_area
	MPI_Bcast(&accel_area, 1, mpi_accel_area, MASTER, MPI_COMM_WORLD);


	// Set stripe and buffer direction and size
	params.grid_fat = ((params.nx - params.ny) >= -200)? 1:0; /* 1 if the grid is square or fat, 0 if it's tall */
	//int params.local_nrows;    // Number of rows in the current rank         Done in params
	//int params.local_ncols;    // Number of cols in the current rank
	// int last_nrows = -1;     // Number of rows in the last rank
	// int last_ncols = -1;     // Number of cols in the last rank
	int prev = (params.my_rank == 0)? (com_size-1) : params.my_rank-1;
	int next = (params.my_rank+1) % com_size;
	// printf("Ramk=%d, prev=%d, next=%d\n", params.my_rank, prev, next);

	// Declare buffers
	speed_t* local_work_space = NULL;
	speed_t* local_temp_space = NULL;
	speed_t* send_buff = NULL;
	speed_t* read_buff = NULL;
	int*     local_obstacles = NULL;


	calculate_local_stripes(&params, com_size, &send_buff, &read_buff,
	 &local_work_space, &local_temp_space, &local_obstacles);
	// printf("Rank:%d params.local_nrows=%d params.local_ncols=%d\n", params.my_rank, params.local_nrows, params.local_ncols);

	// printf("Rank:%d\nnx=%d ny=%d maxIt=%d reyDim=%d lrows=%d lcols=%d density=%f accel=%f omega=%f\n\n"
	// 	, params.my_rank, params.nx, params.ny, params.max_iters, params.reynolds_dim, params.local_nrows, params.local_ncols
	// 	, params.density, params.accel, params.omega);
	// MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

	// Create mpi type for spee_t
	MPI_Datatype mpi_speed_t;
	MPI_Type_create_struct(1, (int[1]){9}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_DOUBLE}, &mpi_speed_t);
	MPI_Type_commit(&mpi_speed_t);

	MPI_Datatype mpi_row;
	MPI_Datatype mpi_vels_North;
	MPI_Datatype mpi_vels_South;
	// Could maybe take the exchange out of the if and just use the row type generally.

	// scatter, shift and exchange based on striping
	if (params.grid_fat!=0)
	{               // Row stripes are contigous in C memory
		// create data type for
		MPI_Type_contiguous(params.local_ncols, mpi_speed_t, &mpi_row);
		MPI_Type_commit(&mpi_row);

		MPI_Scatter(cells, params.local_nrows, mpi_row,
		 &local_work_space[params.local_ncols], params.local_nrows, mpi_row,
		  MASTER, MPI_COMM_WORLD);

		// Shift cells out of the halo space
		// for (int k = (params.local_ncols*params.local_nrows)-1; k >= params.local_ncols; k--)
		// {
		//     local_work_space[params.local_ncols+k]= local_work_space[k];
		//     // local_work_space[k] = (speed_t){0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		// }    

		MPI_Scatter(obstacles, params.local_ncols*(params.local_nrows), MPI_INT,
		 &local_obstacles[params.local_ncols], params.local_ncols*(params.local_nrows), MPI_INT,
		  MASTER, MPI_COMM_WORLD);

		// Shift obstacles out of the halo space
		// for (int k = (params.local_ncols*params.local_nrows)-1; k >= params.local_ncols; k--)
		// {
		//     local_obstacles[params.local_ncols+k]= local_obstacles[k];
		//     // local_work_space[k] = (speed_t){0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		// }

		// Exchange back
		// MPI_Sendrecv(&local_work_space[params.local_ncols], 1, mpi_row, prev, INITIALISE,
		//  &local_work_space[params.local_ncols*(params.local_nrows+1)], 1, mpi_row, next, INITIALISE,
		//   MPI_COMM_WORLD, &status);

		// // Exchange forward
		// MPI_Sendrecv(&local_work_space[params.local_ncols*params.local_nrows], 1, mpi_row, next, INITIALISE,
		//  local_work_space, 1, mpi_row, prev, INITIALISE,
		//    MPI_COMM_WORLD, &status);

		int* blok_len = malloc(sizeof(int)*params.local_ncols*2);
		int* vel_disps = malloc(sizeof(MPI_Aint)*params.local_ncols*2);
		// MPI_Aint base;           
		blok_len[0] = 1;
		vel_disps[0] = 2;
		blok_len[1] = 2;
		vel_disps[1] = 5;
		for (int t = 2; t < params.local_ncols*2; ++t)
		{
			blok_len[t]   = 1;
			vel_disps[t]  = vel_disps[t-2]+9;

			blok_len[++t] = 2;
			vel_disps[t]  = vel_disps[t-2]+9;

		}
		/*MPI_Get_address(&(local_temp_space[params.local_ncols*(params.local_nrows+1)].speeds[2]), &vel_disps[0]);
		blok_len[0] = 1;
		MPI_Get_address(&(local_temp_space[params.local_ncols*(params.local_nrows+1)].speeds[5]), &vel_disps[1]);
		vel_disps[1] -= vel_disps[0];
		blok_len[1] = 1;
		MPI_Get_address(&(local_temp_space[params.local_ncols*(params.local_nrows+1)].speeds[6]), &vel_disps[2]);
		vel_disps[2] -= vel_disps[0];
		blok_len[2] = 1;
		// if (params.my_rank==0) {printf("blen[0]=1\ndips[%d]=%d\ndips[%d]=%d\ndips[%d]=%d\n"
			// , 0,vel_disps[0],1,vel_disps[1],2,vel_disps[2]);}
		for (int t = 1; t < params.local_ncols; ++t)
		{
			blok_len[t*3] = 1;
			MPI_Get_address(&(local_temp_space[t+(params.local_ncols*(params.local_nrows+1))].speeds[2]), &vel_disps[t*3]);
			vel_disps[t*3] -= vel_disps[0];
			blok_len[(t*3)+1] = 1;
			MPI_Get_address(&(local_temp_space[t+(params.local_ncols*(params.local_nrows+1))].speeds[5]), &vel_disps[(t*3)+1]);
			vel_disps[(t*3)+1] -= vel_disps[0];
			blok_len[(t*3)+2] = 1;
			MPI_Get_address(&(local_temp_space[t+(params.local_ncols*(params.local_nrows+1))].speeds[6]), &vel_disps[(t*3)+2]);
			vel_disps[(t*3)+2] -= vel_disps[0];
			// if (params.my_rank==0) {printf("dips[%d]=%d\ndips[%d]=%d\ndips[%d]=%d\n"
			// ,t*3,vel_disps[t*3],(t*3)+1,vel_disps[(t*3)+1],(t*3)+2,vel_disps[(t*3)+2]);}
		}
		MPI_Get_address(&local_temp_space[params.local_ncols*(params.local_nrows+1)], &base);
		vel_disps[0] -= base;*/

		MPI_Type_indexed(params.local_ncols*2, blok_len,
			vel_disps, MPI_DOUBLE,
			&mpi_vels_North);

		MPI_Type_commit(&mpi_vels_North);

		blok_len[0] = 1;
		vel_disps[0] = 4;
		blok_len[1] = 2;
		vel_disps[1] = 7;
		for (int t = 2; t < params.local_ncols*2; ++t)
		{
			blok_len[t]   = 1;
			vel_disps[t]  = vel_disps[t-2]+9;

			blok_len[++t] = 2;
			vel_disps[t]  = vel_disps[t-2]+9;

		}

		/*MPI_Get_address(&(local_temp_space[0].speeds[4]), &vel_disps[0]);
		blok_len[0] = 1;
		MPI_Get_address(&(local_temp_space[0].speeds[7]), &vel_disps[1]);
		vel_disps[1] -= vel_disps[0];
		blok_len[1] = 1;
		MPI_Get_address(&(local_temp_space[0].speeds[8]), &vel_disps[2]);
		vel_disps[2] -= vel_disps[0];
		blok_len[2] = 1;
		// if (params.my_rank==0) {printf("blen[0]=1\ndips[%d]=%d\ndips[%d]=%d\ndips[%d]=%d\n"
			// , 0,vel_disps[0],1,vel_disps[1],2,vel_disps[2]);}
		for (int t = 1; t < params.local_ncols; ++t)
		{
			blok_len[t*3] = 1;
			MPI_Get_address(&(local_temp_space[t+(0)].speeds[4]), &vel_disps[t*3]);
			vel_disps[t*3] -= vel_disps[0];
			blok_len[(t*3)+1] = 1;
			MPI_Get_address(&(local_temp_space[t+(0)].speeds[7]), &vel_disps[(t*3)+1]);
			vel_disps[(t*3)+1] -= vel_disps[0];
			blok_len[(t*3)+2] = 1;
			MPI_Get_address(&(local_temp_space[t+(0)].speeds[8]), &vel_disps[(t*3)+2]);
			vel_disps[(t*3)+2] -= vel_disps[0];
			// if (params.my_rank==0) {printf("dips[%d]=%d\ndips[%d]=%d\ndips[%d]=%d\n"
			// ,t*3,vel_disps[t*3],(t*3)+1,vel_disps[(t*3)+1],(t*3)+2,vel_disps[(t*3)+2]);}
		}
		MPI_Get_address(&local_temp_space[params.local_ncols*(params.local_nrows+1)], &base);
		vel_disps[0] -= base;*/

		MPI_Type_indexed(params.local_ncols*2, blok_len,
			vel_disps, MPI_DOUBLE,
			&mpi_vels_South);

		MPI_Type_commit(&mpi_vels_South);

		// for (int i = 0; i < 6; ++i)
		// if (params.my_rank==1)
		// printf("blok_len[%d]=%d\nvel_disps[%d]=%d\n",
		//   i, blok_len[i],  i, vel_disps[i]);

		free(blok_len);
		free(vel_disps);
	}
	else
	{               // Columb stripes need a strided data type
		// Vector
		// MPI_Type_commit(&mpi_last_row);

	}


	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

	int acel_row_rank;
	if (params.grid_fat!=0)
	{
		// Calculate the rank and nex index that the acel row exists in
		if (accel_area.col_or_row==ACCEL_ROW){
			int lrows = params.ny/com_size;
			acel_row_rank = accel_area.idx /lrows;
			if ((accel_area.idx %( lrows)) !=0)  // If there's a remainder, check if the rank should be imcremented.
			{ 
				acel_row_rank = ((accel_area.idx /(lrows))== com_size)? com_size : acel_row_rank; 
				// Adjust the index for the local rank
				if (acel_row_rank>0)
				{
					accel_area.idx -= (acel_row_rank!=com_size)? (acel_row_rank)*lrows:(acel_row_rank-1)*lrows;;
					if (accel_area.idx<1) accel_area.idx = 1;
				}
			}
			else
			{
				if (acel_row_rank!=0)
				{
					acel_row_rank--;
					accel_area.idx -= (acel_row_rank)*2;
				}
			}

		}
	}
	else
	{
		// cacl row rank for col-wise
	}

	if (params.my_rank==MASTER)
	{
		/* iterate for max_iters timesteps */
		gettimeofday(&timstr,NULL);
		tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	}

	if (params.grid_fat!=0)
	{       // Default striping is row-wise
		for (ii = 0; ii < params.max_iters; ii++)
		{
			// timestep(params, accel_area, cells, tmp_cells, obstacles);
			if(accel_area.col_or_row==ACCEL_COLUMN)    // Can send through all ranks
			{ 
				accelerate_flow_Colum_RW(params,accel_area,local_work_space,local_obstacles);
				// don't need to acel or exchane halo here, as prop only pushes non halo to temp
			}
			else                                // The row is in one rank only, uses the calculated rank and index
			{if (params.my_rank==acel_row_rank) { accelerate_flow_Row_RW(params,accel_area,local_work_space,local_obstacles); }}

			propagate_row_wise(params,local_work_space,local_temp_space);

			// Individual velocities moved betweenn cells, so have to recieve temp cells and manually replace velocities
			
			if (params.my_rank==1)
			{
				for (int i = 4; i < params.local_ncols-96; ++i)
				{
					// local_temp_space[i].speeds[2] = 5.1f;
					printf("South_Source:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
					   , i, local_temp_space[i].speeds[6]
					   , local_temp_space[i].speeds[2]
					   , local_temp_space[i].speeds[5]
					   , local_temp_space[i].speeds[3]
					   , local_temp_space[i].speeds[0]
					   , local_temp_space[i].speeds[1]
					   , local_temp_space[i].speeds[7]
					   , local_temp_space[i].speeds[4]
					   , local_temp_space[i].speeds[8]);
				}
			}
			if (params.my_rank==0)
			{
				for (int i = 4; i < params.local_ncols-96; ++i)
				{
					printf("South_Before:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
					   , i, local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[6]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[2]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[5]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[3]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[0]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[1]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[7]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[4]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[8]);
				}
			}

			// Exchange temp halo back
			MPI_Sendrecv(&local_temp_space[0], 1, mpi_vels_South, prev, HALO_VELS,
			 &local_temp_space[params.local_ncols*params.local_nrows], 1, mpi_vels_South, next, HALO_VELS,
			  MPI_COMM_WORLD, &status);

			if (params.my_rank==0)
			{
				for (int i = 4; i < params.local_ncols-96; ++i)
				{
					printf("South_After:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
					, i, local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[6]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[2]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[5]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[3]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[0]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[1]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[7]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[4]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[8]);
				}
			}

// MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

			/*if (params.my_rank==0)
			{
			    for (int i = 49; i < params.local_ncols-51; ++i)
			    {
			        // local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[7] = 5.1f;
			        printf("North_Source:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
			           , i, local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[6]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[2]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[5]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[3]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[0]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[1]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[7]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[4]
			           , local_temp_space[i+params.local_ncols*(params.local_nrows+1)].speeds[8]);
			    }
			}
			if (params.my_rank==1)
			{
			    for (int i = 49; i < params.local_ncols-51; ++i)
			    {
			        printf("North_Before:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
			           , i, local_temp_space[i+params.local_ncols].speeds[6]
			           , local_temp_space[i+params.local_ncols].speeds[2]
			           , local_temp_space[i+params.local_ncols].speeds[5]
			           , local_temp_space[i+params.local_ncols].speeds[3]
			           , local_temp_space[i+params.local_ncols].speeds[0]
			           , local_temp_space[i+params.local_ncols].speeds[1]
			           , local_temp_space[i+params.local_ncols].speeds[7]
			           , local_temp_space[i+params.local_ncols].speeds[4]
			           , local_temp_space[i+params.local_ncols].speeds[8]);
			    }
			}*/

			// Exchange temp halo forward
			MPI_Sendrecv(&local_temp_space[params.local_ncols*(params.local_nrows+1)], 1, mpi_vels_North, next, HALO_VELS,
			 &local_temp_space[params.local_ncols], 1, mpi_vels_North, prev, HALO_VELS,
			   MPI_COMM_WORLD, &status);
			

			/*if (params.my_rank==1)
			{
			    for (int i = 49; i < params.local_ncols-51; ++i)
			    {
			        printf("North_After:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
			        , i, local_temp_space[i+params.local_ncols].speeds[6]
			           , local_temp_space[i+params.local_ncols].speeds[2]
			           , local_temp_space[i+params.local_ncols].speeds[5]
			           , local_temp_space[i+params.local_ncols].speeds[3]
			           , local_temp_space[i+params.local_ncols].speeds[0]
			           , local_temp_space[i+params.local_ncols].speeds[1]
			           , local_temp_space[i+params.local_ncols].speeds[7]
			           , local_temp_space[i+params.local_ncols].speeds[4]
			           , local_temp_space[i+params.local_ncols].speeds[8]);
			    }
			}*/

// MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

			// Check if should alter for col stipes
			collision_local(params,local_work_space,local_temp_space,local_obstacles);


			if (params.my_rank==0)
			{
				for (int i = 4; i < params.local_ncols-96; ++i)
				{
					printf("Collision_After:Cell:%d\n N_W=%f N=%f N_E=%f\nC_W=%f C_C=%f C_E=%f\ns_W=%f s=%f s_E=%f\n\n\n"
					, i, local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[6]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[2]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[5]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[3]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[0]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[1]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[7]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[4]
					   , local_temp_space[i+(params.local_ncols*params.local_nrows)].speeds[8]);
				}
				printf("%d down.\n\n", ii+1);
			}
// if (ii>0)
// MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

			// Could skip on last it? Do I need grid halos at all?
			// Exchange grid halo back
			// MPI_Sendrecv(&local_work_space[params.local_ncols], 1, mpi_row, prev, HALO_CELLS,
			//  &local_work_space[params.local_ncols*(params.local_nrows+1)], 1, mpi_row, next, HALO_CELLS,
			//   MPI_COMM_WORLD, &status);

			// // Exchange grid halo forward
			// MPI_Sendrecv(&local_work_space[params.local_ncols*params.local_nrows], 1, mpi_row, next, HALO_CELLS,
			//  local_work_space, 1, mpi_row, prev, HALO_CELLS,
			//    MPI_COMM_WORLD, &status);

			// grab vels and div

			//av_vels[ii] = av_velocity(params, cells, obstacles);

			#ifdef DEBUG
			printf("==timestep: %d==\n", ii);
			printf("av velocity: %.12E\n", av_vels[ii]);
			printf("tot density: %.12E\n", total_density(params, cells));
			#endif
		}

		// MPI_Barrier(MPI_COMM_WORLD);

		MPI_Gather(&local_temp_space[params.local_ncols], (params.local_nrows), mpi_row
			, cells, (params.local_nrows), mpi_row,
			 MASTER, MPI_COMM_WORLD);
	}
	else
	{       // colwise stripes
		
	}

	if (params.my_rank==MASTER)
	{

		gettimeofday(&timstr,NULL);
		toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
		getrusage(RUSAGE_SELF, &ru);
		timstr=ru.ru_utime;
		usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
		timstr=ru.ru_stime;
		systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

		printf("==done==\n");
		printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,cells,obstacles));
		printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
		printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
		printf("Elapsed system CPU time:\t%.6f (s)\n", systim);

		write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels);
	}
	free(send_buff);
	free(read_buff);
	free(local_work_space);
	free(local_temp_space);
	free(local_obstacles);
	MPI_Type_free(&mpi_param);
	MPI_Type_free(&mpi_accel_area);
	MPI_Type_free(&mpi_speed_t);
	MPI_Type_free(&mpi_row);
	MPI_Type_free(&mpi_vels_North);
	MPI_Type_free(&mpi_vels_South);

	// MPI_Barrier(MPI_COMM_WORLD);
	
	// MPI_Type_free(&mpi_last_row);
	if (params.my_rank==MASTER)
	{
		finalise(&cells, &tmp_cells, &obstacles, &av_vels);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

void write_values(const char * final_state_file, const char * av_vels_file,
	const param_t params, speed_t* cells, int* obstacles, double* av_vels)
{
	FILE* fp;                     /* file pointer */
	int ii,jj,kk;                 /* generic counters */
	const double c_sq = 1.0/3.0;  /* sq. of speed of sound */
	double local_density;         /* per grid cell sum of densities */
	double pressure;              /* fluid pressure in grid cell */
	double u_x;                   /* x-component of velocity in grid cell */
	double u_y;                   /* y-component of velocity in grid cell */
	double u;                     /* norm--root of summed squares--of u_x and u_y */

	fp = fopen(final_state_file, "w");

	if (fp == NULL)
	{
		DIE("could not open file output file");
	}

	for (ii = 0; ii < params.ny; ii++)
	{
		for (jj = 0; jj < params.nx; jj++)
		{
			/* an occupied cell */
			if (obstacles[ii*params.nx + jj])
			{
				u_x = u_y = u = 0.0;
				pressure = params.density * c_sq;
			}
			/* no obstacle */
			else
			{
				local_density = 0.0;

				for (kk = 0; kk < NSPEEDS; kk++)
				{
					local_density += cells[ii*params.nx + jj].speeds[kk];
				}

				/* compute x velocity component */
				u_x = (cells[ii*params.nx + jj].speeds[1] +
						cells[ii*params.nx + jj].speeds[5] +
						cells[ii*params.nx + jj].speeds[8]
					- (cells[ii*params.nx + jj].speeds[3] +
						cells[ii*params.nx + jj].speeds[6] +
						cells[ii*params.nx + jj].speeds[7]))
					/ local_density;

				/* compute y velocity component */
				u_y = (cells[ii*params.nx + jj].speeds[2] +
						cells[ii*params.nx + jj].speeds[5] +
						cells[ii*params.nx + jj].speeds[6]
					- (cells[ii*params.nx + jj].speeds[4] +
						cells[ii*params.nx + jj].speeds[7] +
						cells[ii*params.nx + jj].speeds[8]))
					/ local_density;

				/* compute norm of velocity */
				u = sqrt((u_x * u_x) + (u_y * u_y));

				/* compute pressure */
				pressure = local_density * c_sq;
			}

			/* write to file */
			fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %d\n",
				jj,ii,u_x,u_y,u,pressure,obstacles[ii*params.nx + jj]);
		}
	}

	fclose(fp);

	fp = fopen(av_vels_file, "w");
	if (fp == NULL)
	{
		DIE("could not open file output file");
	}

	for (ii = 0; ii < params.max_iters; ii++)
	{
		fprintf(fp,"%d:\t%.12E\n", ii, av_vels[ii]);
	}

	fclose(fp);
}

double calc_reynolds(const param_t params, speed_t* cells, int* obstacles)
{
	const double viscosity = 1.0 / 6.0 * (2.0 / params.omega - 1.0);

	return av_velocity(params,cells,obstacles) * params.reynolds_dim / viscosity;
}

double total_density(const param_t params, speed_t* cells)
{
	int ii,jj,kk;        /* generic counters */
	double total = 0.0;  /* accumulator */

	for (ii = 0; ii < params.ny; ii++)
	{
		for (jj = 0; jj < params.ny; jj++)
		{
			for (kk = 0; kk < NSPEEDS; kk++)
			{
				total += cells[ii*params.nx + jj].speeds[kk];
			}
		}
	}

	return total;
}

void calculate_local_stripes(param_t* params, int com_size, speed_t** send_buff
	, speed_t** read_buff, speed_t** local_work_space, speed_t** local_temp_space, int** local_obstacles)
{
	if (params->grid_fat!=0)   // Row-wise
	{
		params->local_ncols = params->nx;     // Each row has all cols
		params->local_nrows = params->ny / com_size;

		if (params->my_rank == com_size-1)
		{
			params->local_nrows += (params->ny % com_size);
		}

		if (params->local_nrows <1)
		{
			fprintf(stderr, "Error: Too many processes, local_nrows<1\n");
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		
		// Allocate message bufferes and the local workspace with a halo.
		*send_buff = (speed_t*) malloc(sizeof(speed_t)*params->nx);
		if (*send_buff == NULL) DIE("Cannot allocate memory for the send buffer");
		*read_buff = (speed_t*) malloc(sizeof(speed_t)*params->nx);
		if (*read_buff == NULL) DIE("Cannot allocate memory for the read buffer");
		*local_work_space = (speed_t*) malloc(params->nx*(params->local_nrows+2)*sizeof(speed_t));
		if (*local_work_space == NULL) DIE("Cannot allocate memory for the local work space");
		*local_temp_space = (speed_t*) malloc(sizeof(speed_t)*params->nx*(params->local_nrows+2));
		if (*local_temp_space == NULL) DIE("Cannot allocate memory for the local temp space");
		*local_obstacles = (int*) malloc(sizeof(int)*params->nx*((params->local_nrows)*2));
		if (*local_obstacles == NULL) DIE("Cannot allocate memory for the local obstacles");
		
	}
	else            // Col-wise
	{
		params->local_nrows = params->ny;     // Each col has all rows
		params->local_ncols = params->nx / com_size;
		if (params->my_rank == com_size-1)
		{
			params->local_ncols += (params->nx % com_size);
		}

		if (params->local_ncols <1)
		{
			fprintf(stderr, "Error: Too many processes, local_ncols<1\n");
			MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		}
		
		// Allocate message bufferes and the local workspaces with a halo.
		*send_buff = (speed_t*) malloc(sizeof(speed_t)*params->ny);
		if (*send_buff == NULL) DIE("Cannot allocate memory for the send buffer");
		*read_buff = (speed_t*) malloc(sizeof(speed_t)*params->ny);
		if (*read_buff == NULL) DIE("Cannot allocate memory for the read buffer");
		*local_work_space = (speed_t*) malloc(sizeof(speed_t)*params->ny*(params->local_ncols+2));
		if (*local_work_space == NULL) DIE("Cannot allocate memory for the local work space");
		*local_temp_space = (speed_t*) malloc(sizeof(speed_t)*params->ny*(params->local_ncols+2));
		if (*local_temp_space == NULL) DIE("Cannot allocate memory for the local temp space");
		*local_obstacles = (int*) malloc(sizeof(int)*params->ny*((params->local_ncols)*2));
		if (*local_obstacles == NULL) DIE("Cannot allocate memory for the local obstacles");
	}
}
