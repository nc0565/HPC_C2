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
    MPI_Comm_rank(MPI_COMM_WORLD, &params.my_rank);

    // Initialise
    if (params.my_rank == MASTER)
    {
        parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);
        initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels);         
    }

    // Create mpi type for param
    int blocklengths[2] = {7,3};
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
    int grid_fat = ((params.nx - params.ny) >= -200)? 1:0; /* 1 if the grid is square or fat, 0 if it's tall */
    //int params.local_nrows;    // Number of rows in the current rank         Done in params
    //int params.local_ncols;    // Number of cols in the current rank
    // int last_nrows = -1;     // Number of rows in the last rank
    // int last_ncols = -1;     // Number of cols in the last rank
    int prev = (params.my_rank == 0)? (com_size-1) : params.my_rank-1;
    int next = (params.my_rank+1) % com_size;
    // printf("Ramk=%d, prev=%d, next=%d\n", params.my_rank, prev, next);


    speed_t* send_buff = NULL;
    speed_t* read_buff = NULL;
    speed_t* local_work_space = NULL;
    speed_t* local_temp_space = NULL;
    int*     local_obstacles = NULL;


    calculate_local_stripes(&params, com_size, grid_fat,
     &send_buff, &read_buff, &local_work_space, &local_temp_space, &local_obstacles);
    // printf("Rank:%d params.local_nrows=%d params.local_ncols=%d\n", params.my_rank, params.local_nrows, params.local_ncols);

    // Create mpi type for spee_t
    MPI_Datatype mpi_speed_t;
    MPI_Type_create_struct(1, (int[1]){9}, (MPI_Aint[1]){0}, (MPI_Datatype[1]){MPI_DOUBLE}, &mpi_speed_t);
    MPI_Type_commit(&mpi_speed_t);

    // Refactor with if(fat_grid)
    // Create mpi type for rows
    // if (params.my_rank!=com_size-1)
    // {
        MPI_Datatype mpi_row;
        MPI_Type_contiguous(params.local_ncols, mpi_speed_t, &mpi_row);
        MPI_Type_commit(&mpi_row);
    // }
    // else
    // {
    //     MPI_Datatype mpi_last_row;
    //     MPI_Type_contiguous(params.local_ncols, mpi_speed_t, &mpi_last_row);
    //     MPI_Type_commit(&mpi_last_row);
    // }

    // scatter, shift and exchange based on striping
    if (grid_fat!=0)
    {               // Row stripes are contigous in C memory
        MPI_Scatter(cells, params.local_nrows*params.local_ncols, mpi_speed_t,
         local_work_space, params.local_nrows*params.local_ncols, mpi_speed_t,
          MASTER, MPI_COMM_WORLD);

        // Shift out of the halo space
        for (int k = (params.local_ncols*params.local_nrows)-1; k >= params.local_ncols; k--)
        {
            local_work_space[params.local_ncols+k]= local_work_space[k];
            // local_work_space[k] = (speed_t){0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
        }    

        MPI_Scatter(obstacles, params.local_ncols*params.local_nrows, MPI_INT,
         local_obstacles, params.local_ncols*params.local_nrows, MPI_INT,
          MASTER, MPI_COMM_WORLD);


        // if(params.my_rank==0)
        // {
        //     printf("Rank%d nx=%d ny= %d\n\n\n\n", params.my_rank, params.local_ncols, params.local_nrows);
        //     for (int i = 0; i < params.local_ncols*(params.local_nrows+2); ++i)
        //     {
        //         for (int k = 0; k < 8; ++k)
        //         printf("%f", params.my_rank, local_work_space[i].speeds[k]);

        //         printf("%f_ _", local_work_space[i].speeds[8]);
        //     }
        //     printf("working here\n\n\n\n\n\n\n\n\n\n");
        // }
        
        // if (params.my_rank!=com_size-1)
        // {


            MPI_Sendrecv(&local_work_space[params.local_ncols], 1, mpi_row, prev, HALO,
             &local_work_space[params.local_ncols*(params.local_nrows+1)], 1, mpi_row, next, HALO,
              MPI_COMM_WORLD, &status);


for (int i = 0; i < params.local_ncols; ++i)
    {
        if (params.my_rank==0)
        {
            for (int k = 0; k < 9; ++k)
                printf("%f", local_work_space[(params.local_ncols*params.local_nrows)+i].speeds[k]);

                printf("_ _");
        }
    }
    if (params.my_rank==0)
        {printf("working jgjgjgv  here\n\n\n\n\n\n\n\n\n\n");}

MPI_Barrier(MPI_COMM_WORLD);

            MPI_Sendrecv(&local_work_space[params.local_ncols*params.local_nrows], 1, mpi_row, next, HALO,
             local_work_space, 1, mpi_row, prev, HALO,
              MPI_COMM_WORLD, &status);


for (int i = 0; i < params.local_ncols; ++i)
    {
        if (params.my_rank==1)
        {
            for (int k = 0; k < 9; ++k)
                printf("%f", local_work_space[i].speeds[k]);

                printf("_ _");
        }
    }
    if (params.my_rank==1)
        {printf("working dgrhthtjh  here\n\n\n\n\n\n\n\n\n\n");}

MPI_Barrier(MPI_COMM_WORLD);

    /*for (int i = 0; i < params.local_ncols; ++i)
    {
        send_buff[i] = local_work_space[(params.local_ncols*params.local_nrows)+i];
        if (params.my_rank==0)
        {
            for (int k = 0; k < 9; ++k)
                printf("%f", send_buff[i].speeds[k]);

                printf("_ _");
        }
    }
    if (params.my_rank==0)
        {printf("working jgjgjgv  here\n\n\n\n\n\n\n\n\n\n");}

    MPI_Sendrecv(send_buff, 128, mpi_speed_t, next, HALO,
             read_buff, 128, mpi_speed_t, prev, HALO,
              MPI_COMM_WORLD, &status);

    MPI_Barrier(MPI_COMM_WORLD);
            printf("Rank:%d params.local_nrows=%d params.local_ncols=%d\n", params.my_rank, params.local_nrows, params.local_ncols);


    for (int i = 0; i < params.local_ncols; ++i)
    {
        local_work_space[i] = read_buff[i];
        if (params.my_rank==1)
        {
            for (int k = 0; k < 9; ++k)
                printf("%f", read_buff[i].speeds[k]);

                printf("_ _");
        }
    }
    if (params.my_rank==1)
        {printf("working dgrhthtjh  here\n\n\n\n\n\n\n\n\n\n");}

MPI_Barrier(MPI_COMM_WORLD);



MPI_Barrier(MPI_COMM_WORLD);

printf("lll\n");
MPI_Barrier(MPI_COMM_WORLD);*/
        // }
        // else
        // {
            
        // }

        // if(params.my_rank==0)
        // {
        //     printf("Rank%d nx=%d ny= %d\n\n\n\n", params.my_rank, params.local_ncols, params.local_nrows);
        //     for (int i = 0; i < params.local_ncols*(params.local_nrows+2); ++i)
        //     {
        //         for (int k = 0; k < 9; ++k)
        //         printf("%f", local_work_space[i].speeds[k]);

        //         //printf("_ _");
        //     }
        //     printf("working here\n\n\n\n\n\n\n\n\n\n");
        // }

        /*if(params.my_rank==0)
        {
            printf("Rank%d nx=%d ny= %d\n\n\n\n", params.my_rank, params.local_ncols, params.local_nrows);
            for (int i = 0; i < 128; ++i)
            {
                for (int k = 0; k < 9; ++k)
                printf("%f", local_work_space[i].speeds[k]);

                printf("_ _");
            }
            printf("working jgjgjgv  here\n\n\n\n\n\n\n\n\n\n");
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if(params.my_rank==4)
        {
            printf("Rank%d nx=%d ny= %d\n\n\n\n", params.my_rank, params.local_ncols, params.local_nrows);
            for (int i = 26*128; i < (27*128)-1; ++i)
            {
                for (int k = 0; k < 9; ++k)
                printf("%f", local_work_space[i].speeds[k]);

                printf("_ _");
            }
            printf("working jgjgjgv  here\n\n\n\n\n\n\n\n\n\n");
        }*/
            
    }
    else
    {               // Columb stripes need a strided data type

    }




    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    /* iterate for max_iters timesteps */
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    for (ii = 0; ii < params.max_iters; ii++)
    {
        // timestep(params, accel_area, cells, tmp_cells, obstacles);
        accelerate_flow(params,accel_area,local_work_space,local_obstacles);
        propagate(params,local_work_space,local_temp_space);
        collision(params,local_work_space,local_temp_space,local_obstacles);

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

    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,cells,obstacles));
    printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
    printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
    printf("Elapsed system CPU time:\t%.6f (s)\n", systim);

    write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels);
    free(send_buff);
    free(read_buff);
    free(local_work_space);
    free(local_temp_space);
    free(local_obstacles);
    finalise(&cells, &tmp_cells, &obstacles, &av_vels);

    MPI_Type_free(&mpi_param);
    MPI_Type_free(&mpi_accel_area);
    MPI_Type_free(&mpi_speed_t);
    MPI_Type_free(&mpi_row);
    // MPI_Type_free(&mpi_last_row);
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

void calculate_local_stripes(param_t* params, int com_size, int grid_fat, speed_t** send_buff
    , speed_t** read_buff, speed_t** local_work_space, speed_t** local_temp_space, int** local_obstacles)
{
    if (grid_fat!=0)   // Row-wise
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
        *local_obstacles = (int*) malloc(sizeof(int)*params->nx*(params->local_nrows));
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
        *local_obstacles = (int*) malloc(sizeof(int)*params->ny*(params->local_ncols));
        if (*local_obstacles == NULL) DIE("Cannot allocate memory for the local obstacles");
    }
}
