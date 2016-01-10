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

#include "lbm.h"
#include "mpi.h"

#define MASTER 0

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int com_size;               // The number of ranks in the comunicator

    char * final_state_file = NULL;
    char * av_vels_file = NULL;
    char * param_file = NULL;

    accel_area_t accel_area;

    param_t  params;              /* struct to hold parameter values */
    speed_t* cells     = NULL;    /* grid containing fluid densities */
    speed_t* tmp_cells = NULL;    /* scratch space */
    int*     obstacles = NULL;    /* grid indicating which cells are blocked */
    double*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */

    int    ii;                    /*  generic counter */
    int    tcount =0;             // A count to calculate what to divide av_vels by
                                  // nx*ny - (tcount-9).
    struct timeval timstr;        /* structure to hold elapsed time */
    struct rusage ru;             /* structure to hold CPU time--system and user */
    double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
    double usrtim;                /* floating point number to record elapsed user CPU time */
    double systim;                /* floating point number to record elapsed system CPU time */

    // Check world
    MPI_Comm_size(MPI_COMM_WORLD, &com_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &params.my_rank);

    // Initialise
        parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);
        initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles/*, &av_vels*/, &tcount); 

    // Declare types
    MPI_Datatype mpi_param;
    MPI_Datatype mpi_accel_area;
    MPI_Datatype mpi_speed_t;
    // MPI_Datatype temp_t;
    // Declare buffers
    setup(&params, &accel_area, com_size, &mpi_param, &mpi_accel_area, &mpi_speed_t/*, &temp_t*/);

    speed_t* local_work_space = NULL;
    speed_t* local_temp_space = NULL;
    double*  send_buff = NULL;
    // double*  read_buff = NULL;
    int*     local_obstacles = NULL;
    int grid_shape = calculate_local_stripes(&params, com_size, &send_buff, /*&read_buff,*/
     &local_work_space, &local_temp_space, &local_obstacles);
    // printf("Rank:%d params->local_nrows=%d params->local_ncols=%d\n", params->my_rank, params->local_nrows, params->local_ncols);


// =======================================================================================
    
    int* bl_counts = NULL;
    int* disps = NULL;
    MPI_Datatype mpi_row;
    int acel_row_rank=-1;
    // scatter, shift and exchange based on striping
    if (grid_shape!=0)
    {               // Row stripes are contigous in C memory
        // create data type for
        MPI_Type_contiguous(params.local_ncols, mpi_speed_t, &mpi_row);
        MPI_Type_commit(&mpi_row);


        // MPI_Type_create_resized(speed_t, 0, 1*sizeof(speed_t), &cmpi_row_ype);

        // int MPI_Type_size(MPI_Datatype datatype, int *size);
        // int MPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size);


        if (grid_shape==-1)
        {

        MPI_Scatter(cells, params.local_nrows,  mpi_row,
         local_work_space, params.local_nrows, mpi_row,
          MASTER, MPI_COMM_WORLD);

        MPI_Scatter(obstacles, params.local_ncols*(params.local_nrows), MPI_INT,
         local_obstacles, params.local_ncols*(params.local_nrows), MPI_INT,
          MASTER, MPI_COMM_WORLD);
        }
        else
        {
            bl_counts = (int*) malloc(com_size*sizeof(int));
            disps = (int*) malloc(com_size*sizeof(int));

            disps[0] = 0;
            bl_counts[0] = (params.ny / com_size);
            for (int i = 1; i < com_size-1; ++i)
            {
                bl_counts[i] = (params.ny / com_size);
                disps[i] = disps[i-1] + (params.ny / com_size);
            }
            bl_counts[com_size-1] = (params.ny / com_size) +(params.ny % com_size);
            disps[com_size-1] = disps[com_size-2] + (params.ny / com_size);

            MPI_Scatterv(cells, bl_counts, disps, mpi_row,
             local_work_space, params.local_nrows, mpi_row,
              MASTER, MPI_COMM_WORLD);

            bl_counts[0] *= params.local_ncols;
            for (int i = 1; i < com_size; ++i)
            {
                bl_counts[i] *= params.local_ncols;
                disps[i] *= params.local_ncols;
            }

            MPI_Scatterv(obstacles, bl_counts,  disps, MPI_INT,
             local_obstacles, params.local_ncols*(params.local_nrows), MPI_INT,
              MASTER, MPI_COMM_WORLD);

            bl_counts[0] /= params.local_ncols;
            for (int i = 1; i < com_size; ++i)
            {
                bl_counts[i] /= params.local_ncols;
                disps[i] /= params.local_ncols;
            }
        }

        // Calculate the rank and nex index that the acel row exists in
        // the new row is the crrect row index, including the halo
        if (accel_area.col_or_row==ACCEL_ROW){
            int lrows = params.ny/com_size;
            // printf("idx=%d \n", accel_area.idx);
            acel_row_rank = (accel_area.idx) / lrows;
            if ((accel_area.idx %( lrows)) !=0 || params.ny%com_size !=0)  // If there's a remainder, check if the rank should be imcremented.
            { 
                acel_row_rank = ((accel_area.idx /(lrows))== com_size)? com_size-1 : acel_row_rank; 
                // Adjust the index for the local rank
                if (acel_row_rank>0)
                {
                    // printf("Cacl _acel row %d in rank %d\n", accel_area.idx, acel_row_rank);
                    accel_area.idx -= (acel_row_rank)*lrows;
                    // if (accel_area.idx<1) accel_area.idx = 1;
                }
            }
            else
            {
                if (acel_row_rank!=0)
                {
                    acel_row_rank--;
                    accel_area.idx -= (acel_row_rank)*params.local_nrows;
                }
            }

        }
        // accel_area.idx = 124;
        // accel_area.idx *= (params.local_nrows/BOX_Y_SIZE);
        // accel_area.idx = 79;
        // acel_row_rank = 3;
        // printf("acel row %d in rank %d\n", accel_area.idx, acel_row_rank);

    }
    else {printf("Not implemented\n");}

// ==========================================================

        /* iterate for max_iters timesteps */
        gettimeofday(&timstr,NULL);
        tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

        av_vels = (double*) malloc(sizeof(double)*(params.max_iters));
        if (av_vels == NULL) DIE("Cannot allocate memory for av_vels");

    int temp=-1;
    if (grid_shape!=0)
    { 
        // This only deals with the rowise stripes
        for (ii = 0; ii < params.max_iters; ii++)
        {
            //timestep(params, accel_area, cells, tmp_cells, obstacles);
            //av_vels[ii] = av_velocity(params, cells, obstacles);

            if(accel_area.col_or_row==ACCEL_COLUMN)    // Can send through all ranks
            { 
                // don't need to acel or exchane halo here, as prop only pushes non halo to temp
                accelerate_flow_Colum_RW(params,accel_area,local_work_space,local_obstacles);
            }
            else                                // The row is in one rank only, uses the calculated rank and index
            {
                // printf("Rank %d here2\n", params.my_rank);
                if (params.my_rank==acel_row_rank)
                 { 
                    // printf("Rank %d here4\n=================\n", params.my_rank);
                    accelerate_flow_Row_RW(params,accel_area,local_work_space,local_obstacles); 
                }
            }

            // Includes the two buffered halo steps
            propagate_row_wise2(params, local_work_space, local_temp_space, local_obstacles, send_buff/*, read_buff*/);
            // Versions that skip inner obstacles
            // propagate_row_wise3(params, local_work_space, local_temp_space, local_obstacles, send_buff/*, read_buff*/);

            collision_local(params, local_work_space, local_temp_space, local_obstacles);
            // Versions that skip inner obstacles
            // collision_local2(params, local_work_space, local_temp_space, local_obstacles);
            double hold;
            av_velocity_local(params, local_work_space, local_obstacles, &hold, &temp/*, double* av_buff*/);
                // if (params.my_rank==0 /*&& ii ==params.max_iters-1*/)
                // printf("Rank:%d temp=%d tcount=%d hold=%f\n", params.my_rank, temp, tcount, hold);
                av_vels[ii] = hold / (double) temp;

    // if (ii>=0)
    // MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // printf("\n");
            #ifdef DEBUG
            printf("==timestep: %d==\n", ii);
            printf("av velocity: %.12E\n", av_vels[ii]);
            printf("tot density: %.12E\n", total_density(params, cells));
            #endif
        }
            // if (params.my_rank==0)
            // printf("Rank:%d tcount=%d\n", params.my_rank, tcount);
    // ================================================================

        if (grid_shape==-1)
        {
            MPI_Gather(local_temp_space, params.local_nrows, mpi_row
                    , cells, params.local_nrows, mpi_row,
                     MASTER, MPI_COMM_WORLD);
        }
        else
        {   
            MPI_Gatherv(local_temp_space, (params.local_nrows), mpi_row
                    , cells, bl_counts, disps, mpi_row,
                     MASTER, MPI_COMM_WORLD);
        }

    }
    else {printf("Not implemented\n");}

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
        finalise(&cells, &tmp_cells, &obstacles/*, &av_vels*/);
        free(av_vels);

    free(bl_counts);
    free(disps);
    MPI_Type_free(&mpi_param);
    MPI_Type_free(&mpi_accel_area);
    MPI_Type_free(&mpi_speed_t);
    // MPI_Type_free(&temp_t);

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

void setup(param_t* params, accel_area_t* accel_area, int com_size,
 MPI_Datatype* mpi_param, MPI_Datatype* mpi_accel_area, MPI_Datatype* mpi_speed_t/*, MPI_Datatype* temp_t*/)
{
    // Create mpi type for param
    int blocklengths[2] = {9,3};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Aint displacements[2];
    MPI_Get_address(&(params->nx), &displacements[0]);
    MPI_Get_address(&(params->density), &displacements[1]);
    displacements[1] -= displacements[0];
    displacements[0] = 0;

    MPI_Type_create_struct(2, blocklengths, displacements, types, mpi_param);
    MPI_Type_commit(mpi_param);

    // Broadcast params
    MPI_Bcast(params, 1, *mpi_param, MASTER, MPI_COMM_WORLD);

    // Restore rank numbers
    MPI_Comm_rank(MPI_COMM_WORLD, &(params->my_rank));

    MPI_Type_create_struct(1, (int[1]){2}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_INT}, mpi_accel_area);
    MPI_Type_commit(mpi_accel_area);

    // Broadcast accel_area
    MPI_Bcast(accel_area, 1, *mpi_accel_area, MASTER, MPI_COMM_WORLD);

    //int params->local_nrows;    // Number of rows in the current rank         Done in params
    //int params->local_ncols;    // Number of cols in the current rank
    // int last_nrows = -1;     // Number of rows in the last rank
    // int last_ncols = -1;     // Number of cols in the last rank
    params->prev = (params->my_rank == MASTER)? (com_size-1) : params->my_rank-1;
    params->next = (params->my_rank+1) % com_size;
    // printf("Ramk=%d, prev=%d, next=%d\n", params->my_rank, params->prev, params->next);

    // printf("Rank:%d\nnx=%d ny=%d maxIt=%d reyDim=%d lrows=%d lcols=%d density=%f accel=%f omega=%f\n\n"
    //  , params->my_rank, params->nx, params->ny, params->max_iters, params->reynolds_dim, params->local_nrows, params->local_ncols
    //  , params->density, params->accel, params->omega);
    // MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
// ======================================================================================
    // Create mpi type for speed_t
    // speed_t sp = (speed_t){0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    // MPI_Type_create_struct(1, (int[1]){9}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_DOUBLE}, temp_t);
    MPI_Type_create_struct(1, (int[1]){9}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_DOUBLE}, mpi_speed_t);
    // MPI_Type_commit(temp_t);
    // MPI_Aint base, start;
    // MPI_Get_address(&sp, &base);
    // MPI_Get_address(&sp.speeds[0], &start);
    // start -= base;
    // MPI_Type_create_resized(*temp_t, start, 1*sizeof(speed_t), mpi_speed_t);

    // int s;
    // MPI_Type_size( *temp_t, &s );
    // MPI_Type_create_resized(*temp_t, 0, s, mpi_speed_t);

    // MPI_Type_create_resized(*temp_t, 0, 1*sizeof(speed_t), mpi_speed_t);
    MPI_Type_commit(mpi_speed_t);
}

int calculate_local_stripes(param_t* params, int com_size, double** send_buff
    , /*double** read_buff,*/ speed_t** local_work_space, speed_t** local_temp_space, int** local_obstacles)
{
    // Set stripe and buffer direction and size
    int grid_shape = ((params->nx - params->ny) >= -200)? 1:0; /* 1 if the grid is square or fat, 0 if it's tall; -1 perfectly divisible.*/
    if (grid_shape!=0)   // Row-wise
    {
        params->local_ncols = params->nx;     // Each row has all cols
        params->local_nrows = params->ny / com_size;

        if (params->my_rank == com_size-1)
        {
            params->local_nrows += (params->ny % com_size);
        }

        if ((params->ny % com_size)==0)
        {
            grid_shape = -1;
        }

        if (params->local_nrows <1)
        {
            fprintf(stderr, "Error: Too many processes, local_nrows<1\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Allocate message bufferes and the local workspace with a halo.
        *send_buff = (double*) malloc(3*sizeof(double)*params->nx);
        if (*send_buff == NULL) DIE("Cannot allocate memory for the send buffer");
        // *read_buff = (double*) malloc(3*sizeof(double)*params->nx);
        // if (*read_buff == NULL) DIE("Cannot allocate memory for the read buffer");
        *local_work_space = (speed_t*) malloc(params->nx*(params->local_nrows)*sizeof(speed_t));
        if (*local_work_space == NULL) DIE("Cannot allocate memory for the local work space");
        *local_temp_space = (speed_t*) malloc(sizeof(speed_t)*params->nx*(params->local_nrows));
        if (*local_temp_space == NULL) DIE("Cannot allocate memory for the local temp space");
        *local_obstacles = (int*) malloc(sizeof(int)*params->nx*((params->local_nrows)));
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

        // Shouldn't be needed.
        if ((params->nx % com_size)==0)
        {
            grid_shape = -1;
        }

        if (params->local_ncols <1)
        {
            fprintf(stderr, "Error: Too many processes, local_ncols<1\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Allocate message bufferes and the local workspaces with a halo.
        *send_buff = (double*) malloc(3*sizeof(double)*params->ny);
        if (*send_buff == NULL) DIE("Cannot allocate memory for the send buffer");
        // *read_buff = (double*) malloc(3*sizeof(double)*params->ny);
        // if (*read_buff == NULL) DIE("Cannot allocate memory for the read buffer");
        *local_work_space = (speed_t*) malloc(sizeof(speed_t)*params->ny*(params->local_ncols+2));
        if (*local_work_space == NULL) DIE("Cannot allocate memory for the local work space");
        *local_temp_space = (speed_t*) malloc(sizeof(speed_t)*params->ny*(params->local_ncols));
        if (*local_temp_space == NULL) DIE("Cannot allocate memory for the local temp space");
        *local_obstacles = (int*) malloc(sizeof(int)*params->ny*((params->local_ncols)));
        if (*local_obstacles == NULL) DIE("Cannot allocate memory for the local obstacles");
    }

    return grid_shape;
}