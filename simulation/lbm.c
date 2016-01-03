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

    int my_rank;                // Current rank
    int com_size;               // The number of ranks in the comunicator
    MPI_Status status;          // Struct for recieve

    param_t params /*= {0,0,0,0,0.0,0.0,0.0}*/;            /* struct to hold parameter values */
    int    ii;                  /*  generic counter */
    struct timeval timstr;      /* structure to hold elapsed time */
    struct rusage ru;           /* structure to hold CPU time--system and user */
    double tic,toc;             /* floating point numbers to calculate elapsed wallclock time */
    double usrtim;              /* floating point number to record elapsed user CPU time */
    double systim;              /* floating point number to record elapsed system CPU time */
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
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == MASTER)
    {
        parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);
        initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels);         
    }

    // Create mpi type for param
    int blocklengths[2] = {4,3};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    MPI_Aint displacements[2];
    MPI_Get_address(&params.nx, &displacements[0]);
    MPI_Get_address(&params.density, &displacements[1]);
    displacements[1] -= displacements[0];
    displacements[0] = 0;

    MPI_Datatype mpi_param;
    MPI_Type_create_struct(2, blocklengths, displacements, types, &mpi_param);
    MPI_Type_commit(&mpi_param);

    // Broadcast parms
    MPI_Bcast(&params, 1, mpi_param, MASTER, MPI_COMM_WORLD);

    MPI_Datatype mpi_accel_area;
    MPI_Type_create_struct(1, (int[1]){2}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_INT}, &mpi_accel_area);
    MPI_Type_commit(&mpi_accel_area);

    // Broadcast accel_area
    MPI_Bcast(&accel_area, 1, mpi_accel_area, MASTER, MPI_COMM_WORLD);


    // Set stripe and buffer direction and size
    int grid_fat = ((params.nx - params.ny) >= -200)? 0:1; /* 0 if the grid is square or fat, 1 if it's tall */
    int local_nrows;    // Number of rows in the current rank
    int local_ncols;    // Number of cols in the current rank
    // int last_nrows = -1;     // Number of rows in the last rank
    // int last_ncols = -1;     // Number of cols in the last rank
    int prev = (my_rank == 0)? (com_size-1) : my_rank-1;
    int next = (my_rank+1) % com_size;
    speed_t* send_buff = NULL;
    speed_t* read_buff = NULL;
    speed_t* local_work_space = NULL;


    calculate_local_stripes(params, &local_nrows, &local_ncols, com_size, my_rank, grid_fat, &send_buff, &read_buff, &local_work_space);
    // printf("Rank:%d local_nrows=%d local_ncols=%d\n", my_rank, local_nrows, local_ncols);

    // Create mpi type for spee_t
    MPI_Datatype mpi_speed_t;
    MPI_Type_create_struct(1, (int[1]){9}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_DOUBLE}, &mpi_speed_t);
    MPI_Type_commit(&mpi_speed_t);

    // scatter based on striping
    if (grid_fat==0)
    {               // Row stripes are contigous in C memory
        MPI_Scatter(cells, local_ncols*local_nrows, mpi_speed_t, local_work_space, local_ncols*local_nrows, mpi_speed_t, MASTER, MPI_COMM_WORLD);
        
        // Shift out of the halo space
        for (int k = local_ncols*local_nrows; k >= 0; k--)
        {
            local_work_space[k+local_ncols]= local_work_space[k];
            // local_work_space[k] = (speed_t){0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        }
        // for (int i = 0; i < local_ncols; ++i)
        // {
        //     printf("Rank%d_%d\n", my_rank, local_work_space[i]);
        // }
        // printf("\n");

    }
    else
    {               // Columb stripes need a strided data type

    }

    // Initial halo

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    /* iterate for max_iters timesteps */
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    for (ii = 0; ii < params.max_iters; ii++)
    {
        // timestep(params, accel_area, cells, tmp_cells, obstacles);
        accelerate_flow(params,accel_area,cells,obstacles);
        propagate(params,cells,tmp_cells);
        collision(params,cells,tmp_cells,obstacles);

        // halo1
        // Halo2
        // grab vels and div


        av_vels[ii] = av_velocity(params, cells, obstacles);

        #ifdef DEBUG
        printf("==timestep: %d==\n", ii);
        printf("av velocity: %.12E\n", av_vels[ii]);
        printf("tot density: %.12E\n", total_density(params, cells));
        #endif
    }

    gettimeofday(&timstr,NULL);
    toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    getrusage(RUSAGE_SELF, &ru);
    timstr=ru.ru_utime;
    usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    timstr=ru.ru_stime;
    systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    MPI_Type_free(&mpi_param);
    MPI_Finalize();
    // Remeber to free buffers

    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,cells,obstacles));
    printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
    printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
    printf("Elapsed system CPU time:\t%.6f (s)\n", systim);

    write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels);
    free(send_buff);
    free(read_buff);
    free(local_work_space);
    finalise(&cells, &tmp_cells, &obstacles, &av_vels);

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

void calculate_local_stripes(param_t params, int* local_nrows, int* local_ncols, int com_size, int my_rank, int grid_fat, speed_t** send_buff, speed_t** read_buff, speed_t** local_work_space)
{
    if (grid_fat==0)   // Row-wise
    {
        *local_ncols = params.nx;     // Each row has all cols
        *local_nrows = params.ny / com_size;

        if (my_rank == com_size-1)
        {
            *local_nrows += (params.ny % com_size);
        }

        if (*local_nrows <1)
        {
            fprintf(stderr, "Error: Too many processes, local_row<1\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Allocate message bufferes and the local workspace with a halo.
        *send_buff = (speed_t*) malloc(sizeof(speed_t)*params.nx);
        if (*send_buff == NULL) DIE("Cannot allocate memory for the send buffer");
        *read_buff = (speed_t*) malloc(sizeof(speed_t)*params.nx);
        if (*read_buff == NULL) DIE("Cannot allocate memory for the read buffer");
        *local_work_space = (speed_t*) malloc(sizeof(speed_t)*params.nx*(*local_nrows+2));
        if (*local_work_space == NULL) DIE("Cannot allocate memory for the local work space");
        
    }
    else            // Col-wise
    {
        *local_nrows = params.ny;     // Each col has all rows
        *local_ncols = params.nx / com_size;
        if (my_rank == com_size-1)
        {
            *local_ncols += (params.nx % com_size);
        }

        if (*local_ncols <1)
        {
            fprintf(stderr, "Error: Too many processes, local_col<1\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Allocate message bufferes and the local workspace with a halo.
        *send_buff = (speed_t*) malloc(sizeof(speed_t)*params.ny);
        if (*send_buff == NULL) DIE("Cannot allocate memory for the send buffer");
        *read_buff = (speed_t*) malloc(sizeof(speed_t)*params.ny);
        if (*read_buff == NULL) DIE("Cannot allocate memory for the read buffer");
        *local_work_space = (speed_t*) malloc(sizeof(speed_t)*params.ny*(*local_ncols+2));
        if (*local_work_space == NULL) DIE("Cannot allocate memory for the local work space");
    }
}
