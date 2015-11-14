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
**   ./lbm -a av_vels.dat -f final_state.dat -p ../inputs/box.params
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "lbm.h"

/*
** main program:
** initialise, timestep loop, finalise
*/
int main(int argc, char* argv[])
{
    char * final_state_file = NULL;
    char * av_vels_file = NULL;
    char * param_file = NULL;

    accel_area_t accel_area;

    param_t  params;              /* struct to hold parameter values */
    speed_t* cells     = NULL;    /* grid containing fluid densities */
    speed_t* tmp_cells = NULL;    /* scratch space */
    int*     obstacles = NULL;    /* grid indicating which cells are blocked */
    float*  av_vels   = NULL, *av_out = NULL;    /* a record of the av. velocity computed for each timestep */

    int    ii;                    /*  generic counter */
    struct timeval timstr;        /* structure to hold elapsed time */
    struct rusage ru;             /* structure to hold CPU time--system and user */
    float tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
    float usrtim;                /* floating point number to record elapsed user CPU time */
    float systim;                /* floating point number to record elapsed system CPU time */

    int device_id;
    lbm_context_t lbm_context;

    parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file, &device_id);
    int obstacles_count=0;
    initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels, &av_out,&obstacles_count);
    opencl_initialise(device_id, params, accel_area, &lbm_context, cells, obstacles, av_out);

    //==============================================
    cl_int err = clEnqueueWriteBuffer(lbm_context.queue, lbm_context.h_cells_buff, CL_TRUE, 0, (sizeof(speed_t)*params.nx*params.ny), cells, 0, NULL, NULL);
    if (err != CL_SUCCESS) DIE("OpenCL error %d writing to h_cells_buff", err);
    err = clEnqueueWriteBuffer(lbm_context.queue, lbm_context.h_obstacles_buff, CL_TRUE, 0, (sizeof(int)*params.ny*params.nx), obstacles, 0, NULL, NULL);
    if (err != CL_SUCCESS) DIE("OpenCL error %d writing to h_cells_buff", err); 
    // const size_t global[2] = {params.ny, params.nx};
    const int num_groups = 4;
    const int group_size = (params.nx *params.ny)/num_groups;
    float *h_p_sum = malloc(sizeof(float)*group_size);

    /* iterate for max_iters timesteps */
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    for (ii = 0; ii < params.max_iters; ii++)
    {
        timestep(params, accel_area, lbm_context, cells, tmp_cells, obstacles);

        err = clEnqueueReadBuffer(lbm_context.queue, lbm_context.h_cells_buff, CL_TRUE, 0, (sizeof(speed_t)*params.nx*params.ny), cells, 0, NULL, NULL);
        // if (err != CL_SUCCESS) DIE("OpenCL error %d Reading back h_cells_buff", err); 
        av_vels[ii] = av_velocity(params, cells, obstacles);
        
        /*err = clEnqueueNDRangeKernel(lbm_context.queue, lbm_context.k_av_vel, 1, NULL, params.ny*params.nx, \
            num_groups, 0, NULL, NULL); // hardcode for testing fast.params
        // err = clEnqueueTask(lbm_context.queue, lbm_context.k_av_vel, 0, NULL, NULL)
        if (err != CL_SUCCESS) DIE("OpenCL error %d: failed to execute av_vel kernel!", err); 
        err = clEnqueueReadBuffer(lbm_context.queue, lbm_context.h_p_sum_buff, CL_TRUE, 0, \
            (sizeof(float)*group_size), h_p_sum, 0, NULL, NULL);
        if (err != CL_SUCCESS) DIE("OpenCL error %d Reading back h_p_sum_buff", err); 

        for (int jj = 0; jj < group_size; jj++) // hardcode for testing fast.params
        {
            av_vels[ii] += h_p_sum[jj];
        }
        av_vels[ii] /= (float)(obstacles_count);*/

        #ifdef DEBUG
        printf("==timestep: %d==\n", ii);
        printf("av velocity: %.12E\n", av_vels[ii]);
        printf("tot density: %.12E\n", total_density(params, cells));
        #endif
    }


    /*err = */clEnqueueReadBuffer(lbm_context.queue, lbm_context.h_cells_buff, CL_TRUE, 0, (sizeof(speed_t)*params.nx*params.ny), cells, 0, NULL, NULL);
    // if (err != CL_SUCCESS) DIE("OpenCL error %d Reading back h_cells_buff", err); 
    //=============================================

    // Do not remove this, or the timing will be incorrect!
    clFinish(lbm_context.queue);

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

    free(h_p_sum);
    write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels);
    finalise(&cells, &tmp_cells, &obstacles, &av_vels);
    opencl_finalise(lbm_context);

    return EXIT_SUCCESS;
}

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, speed_t* cells, int* obstacles, float* av_vels)
{
    FILE* fp;                     /* file pointer */
    int ii,jj,kk;                 /* generic counters */
    const float c_sq = 1.0/3.0;  /* sq. of speed of sound */
    float local_density;         /* per grid cell sum of densities */
    float pressure;              /* fluid pressure in grid cell */
    float u_x;                   /* x-component of velocity in grid cell */
    float u_y;                   /* y-component of velocity in grid cell */
    float u;                     /* norm--root of summed squares--of u_x and u_y */

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

float calc_reynolds(const param_t params, speed_t* cells, int* obstacles)
{
    const float viscosity = 1.0 / 6.0 * (2.0 / params.omega - 1.0);

    return av_velocity(params,cells,obstacles) * params.reynolds_dim / viscosity;
}

float total_density(const param_t params, speed_t* cells)
{
    int ii,jj,kk;        /* generic counters */
    float total = 0.0;  /* accumulator */

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

