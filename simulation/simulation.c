/* Functions pertinent to the outer simulation steps */

#include <math.h>

#include "lbm.h"
#include "mpi.h"

void timestep(const param_t params, const accel_area_t accel_area,
    speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    accelerate_flow(params,accel_area,cells,obstacles);
    propagate(params,cells,tmp_cells);
    rebound(params,cells,tmp_cells,obstacles);
    collision(params,cells,tmp_cells,obstacles);
}

void accelerate_flow(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles)
{
    int ii,jj;     /* generic counters */
    double w1,w2;  /* weighting factors */

    /* compute weighting factors */
    w1 = params.density * params.accel / 9.0;
    w2 = params.density * params.accel / 36.0;

    if (accel_area.col_or_row == ACCEL_COLUMN)
    {
        jj = accel_area.idx;

        for (ii = 0; ii < params.ny; ii++)
        {
            /* if the cell is not occupied and
            ** we don't send a density negative */
            if (!obstacles[ii*params.nx + jj] &&
            (cells[ii*params.nx + jj].speeds[4] - w1) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[7] - w2) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[8] - w2) > 0.0 )
            {
                /* increase 'north-side' densities */
                cells[ii*params.nx + jj].speeds[2] += w1;
                cells[ii*params.nx + jj].speeds[5] += w2;
                cells[ii*params.nx + jj].speeds[6] += w2;
                /* decrease 'south-side' densities */
                cells[ii*params.nx + jj].speeds[4] -= w1;
                cells[ii*params.nx + jj].speeds[7] -= w2;
                cells[ii*params.nx + jj].speeds[8] -= w2;
            }
        }
    }
    else
    {
        ii = accel_area.idx;

        for (jj = 0; jj < params.nx; jj++)
        {
            /* if the cell is not occupied and
            ** we don't send a density negative */
            if (!obstacles[ii*params.nx + jj] &&
            (cells[ii*params.nx + jj].speeds[3] - w1) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[6] - w2) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[7] - w2) > 0.0 )
            {
                /* increase 'east-side' densities */
                cells[ii*params.nx + jj].speeds[1] += w1;
                cells[ii*params.nx + jj].speeds[5] += w2;
                cells[ii*params.nx + jj].speeds[8] += w2;
                /* decrease 'west-side' densities */
                cells[ii*params.nx + jj].speeds[3] -= w1;
                cells[ii*params.nx + jj].speeds[6] -= w2;
                cells[ii*params.nx + jj].speeds[7] -= w2;
            }
        }
    }
}

void propagate(const param_t params, speed_t* cells, speed_t* tmp_cells)
{
    int ii,jj;            /* generic counters */

    /* loop over _all_ cells */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            y_n = (ii + 1) % params.ny;
            x_e = (jj + 1) % params.nx;
            y_s = (ii == 0) ? (ii + params.ny - 1) : (ii - 1);
            x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);
            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii *params.nx + jj].speeds[0]  = cells[ii*params.nx + jj].speeds[0]; /* central cell, */
                                                     /* no movement   */
            tmp_cells[ii *params.nx + x_e].speeds[1] = cells[ii*params.nx + jj].speeds[1]; /* east */
            tmp_cells[y_n*params.nx + jj].speeds[2]  = cells[ii*params.nx + jj].speeds[2]; /* north */
            tmp_cells[ii *params.nx + x_w].speeds[3] = cells[ii*params.nx + jj].speeds[3]; /* west */
            tmp_cells[y_s*params.nx + jj].speeds[4]  = cells[ii*params.nx + jj].speeds[4]; /* south */
            tmp_cells[y_n*params.nx + x_e].speeds[5] = cells[ii*params.nx + jj].speeds[5]; /* north-east */
            tmp_cells[y_n*params.nx + x_w].speeds[6] = cells[ii*params.nx + jj].speeds[6]; /* north-west */
            tmp_cells[y_s*params.nx + x_w].speeds[7] = cells[ii*params.nx + jj].speeds[7]; /* south-west */
            tmp_cells[y_s*params.nx + x_e].speeds[8] = cells[ii*params.nx + jj].speeds[8]; /* south-east */
        }
    }
}

void rebound(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj;  /* generic counters */

    /* loop over the cells in the grid */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* if the cell contains an obstacle */
            if (obstacles[ii*params.nx + jj])
            {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[ii*params.nx + jj].speeds[1] = tmp_cells[ii*params.nx + jj].speeds[3];
                cells[ii*params.nx + jj].speeds[2] = tmp_cells[ii*params.nx + jj].speeds[4];
                cells[ii*params.nx + jj].speeds[3] = tmp_cells[ii*params.nx + jj].speeds[1];
                cells[ii*params.nx + jj].speeds[4] = tmp_cells[ii*params.nx + jj].speeds[2];
                cells[ii*params.nx + jj].speeds[5] = tmp_cells[ii*params.nx + jj].speeds[7];
                cells[ii*params.nx + jj].speeds[6] = tmp_cells[ii*params.nx + jj].speeds[8];
                cells[ii*params.nx + jj].speeds[7] = tmp_cells[ii*params.nx + jj].speeds[5];
                cells[ii*params.nx + jj].speeds[8] = tmp_cells[ii*params.nx + jj].speeds[6];
            }
        }
    }
}

void collision(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj,kk;                 /* generic counters */
    const double c_sq = 1.0/3.0;  /* square of speed of sound */
    const double w0 = 4.0/9.0;    /* weighting factor */
    const double w1 = 1.0/9.0;    /* weighting factor */
    const double w2 = 1.0/36.0;   /* weighting factor */

    double u_x,u_y;               /* av. velocities in x and y directions */
    double u_sq;                  /* squared velocity */
    double local_density;         /* sum of densities in a particular cell */
    double u[NSPEEDS];            /* directional velocities */
    double d_equ[NSPEEDS];        /* equilibrium densities */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* don't consider occupied cells */
            if (!obstacles[ii*params.nx + jj])
            {
                /* compute local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += tmp_cells[ii*params.nx + jj].speeds[kk];
                }

                /* compute x velocity component */
                u_x = (tmp_cells[ii*params.nx + jj].speeds[1] +
                        tmp_cells[ii*params.nx + jj].speeds[5] +
                        tmp_cells[ii*params.nx + jj].speeds[8]
                    - (tmp_cells[ii*params.nx + jj].speeds[3] +
                        tmp_cells[ii*params.nx + jj].speeds[6] +
                        tmp_cells[ii*params.nx + jj].speeds[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (tmp_cells[ii*params.nx + jj].speeds[2] +
                        tmp_cells[ii*params.nx + jj].speeds[5] +
                        tmp_cells[ii*params.nx + jj].speeds[6]
                    - (tmp_cells[ii*params.nx + jj].speeds[4] +
                        tmp_cells[ii*params.nx + jj].speeds[7] +
                        tmp_cells[ii*params.nx + jj].speeds[8]))
                    / local_density;

                /* velocity squared */
                u_sq = u_x * u_x + u_y * u_y;

                /* directional velocity components */
                u[1] =   u_x;        /* east */
                u[2] =         u_y;  /* north */
                u[3] = - u_x;        /* west */
                u[4] =       - u_y;  /* south */
                u[5] =   u_x + u_y;  /* north-east */
                u[6] = - u_x + u_y;  /* north-west */
                u[7] = - u_x - u_y;  /* south-west */
                u[8] =   u_x - u_y;  /* south-east */

                /* equilibrium densities */
                /* zero velocity density: weight w0 */
                d_equ[0] = w0 * local_density * (1.0 - u_sq / (2.0 * c_sq));
                /* axis speeds: weight w1 */
                d_equ[1] = w1 * local_density * (1.0 + u[1] / c_sq
                    + (u[1] * u[1]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[2] = w1 * local_density * (1.0 + u[2] / c_sq
                    + (u[2] * u[2]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[3] = w1 * local_density * (1.0 + u[3] / c_sq
                    + (u[3] * u[3]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[4] = w1 * local_density * (1.0 + u[4] / c_sq
                    + (u[4] * u[4]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                /* diagonal speeds: weight w2 */
                d_equ[5] = w2 * local_density * (1.0 + u[5] / c_sq
                    + (u[5] * u[5]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[6] = w2 * local_density * (1.0 + u[6] / c_sq
                    + (u[6] * u[6]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[7] = w2 * local_density * (1.0 + u[7] / c_sq
                    + (u[7] * u[7]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));
                d_equ[8] = w2 * local_density * (1.0 + u[8] / c_sq
                    + (u[8] * u[8]) / (2.0 * c_sq * c_sq)
                    - u_sq / (2.0 * c_sq));

                /* relaxation step */
                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    cells[ii*params.nx + jj].speeds[kk] = 
                        (tmp_cells[ii*params.nx + jj].speeds[kk] + params.omega * 
                        (d_equ[kk] - tmp_cells[ii*params.nx + jj].speeds[kk]));
                }
            }
        }
    }
}

double av_velocity(const param_t params, speed_t* cells, int* obstacles)
{
    int    ii,jj,kk;       /* generic counters */
    int    tot_cells = 0;  /* no. of cells used in calculation */
    double tot_u;          /* accumulated magnitudes of velocity for each cell */

    double local_density;  /* total density in cell */
    double u_x;            /* x-component of velocity for current cell */
    double u_y;            /* y-component of velocity for current cell */

    /* initialise */
    tot_u = 0.0;

    /* loop over all non-blocked cells */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* ignore occupied cells */
            if (!obstacles[ii*params.nx + jj])
            {
                /* local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += cells[ii*params.nx + jj].speeds[kk];
                }

                /* x-component of velocity */
                u_x = (cells[ii*params.nx + jj].speeds[1] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[8]
                    - (cells[ii*params.nx + jj].speeds[3] +
                        cells[ii*params.nx + jj].speeds[6] +
                        cells[ii*params.nx + jj].speeds[7])) /
                    local_density;

                /* compute y velocity component */
                u_y = (cells[ii*params.nx + jj].speeds[2] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[6]
                    - (cells[ii*params.nx + jj].speeds[4] +
                        cells[ii*params.nx + jj].speeds[7] +
                        cells[ii*params.nx + jj].speeds[8])) /
                    local_density;

                /* accumulate the norm of x- and y- velocity components */
                tot_u += sqrt(u_x*u_x + u_y*u_y);
                /* increase counter of inspected cells */
                ++tot_cells;
            }
        }
    }

    return tot_u / (double)tot_cells;
}

void propagate_row_wise2(const param_t params, speed_t* cells, speed_t* tmp_cells, double* send_buff, double* recv_buff)
{
    int ii=1,jj,addr;            /* generic counters */
    int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
    MPI_Status status;
    
        addr = params.local_ncols;
    for (jj = 0; jj < params.local_ncols; jj++, addr++)
    { 
        y_n = (ii + 1) % (params.local_nrows+2);
        x_e = (jj + 1) % (params.nx);
        x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);
        /* propagate densities to neighbouring cells, following
        ** appropriate directions of travel and writing into
        ** scratch space grid */
        tmp_cells[ii*params.local_ncols + jj].speeds[0]  = cells[addr].speeds[0]; /* central cell, */
                                                 /* no movement   */
        tmp_cells[ii *params.local_ncols + x_e].speeds[1] = cells[addr].speeds[1]; /* east */
        tmp_cells[y_n*params.local_ncols + jj].speeds[2]  = cells[addr].speeds[2]; /* north */
        tmp_cells[ii *params.local_ncols + x_w].speeds[3] = cells[addr].speeds[3]; /* west */
        send_buff[jj]  = cells[addr].speeds[4]; /* south */
        tmp_cells[y_n*params.local_ncols + x_e].speeds[5] = cells[addr].speeds[5]; /* north-east */
        tmp_cells[y_n*params.local_ncols + x_w].speeds[6] = cells[addr].speeds[6]; /* north-west */
        send_buff[jj+1] = cells[addr].speeds[7]; /* south-west */
        send_buff[jj+2] = cells[addr].speeds[8]; /* south-east */
    }

    // Exchange temp halo back
    MPI_Sendrecv(send_buff, params.local_ncols*3, MPI_DOUBLE, params.prev, HALO_VELS,
     recv_buff, params.local_ncols*3, MPI_DOUBLE, params.next, HALO_VELS,
      MPI_COMM_WORLD, &status);

    addr = params.local_nrows*params.local_ncols;
    for (jj = 0; jj < params.local_ncols; jj++, addr++)
    {
        tmp_cells[addr].speeds[4] = recv_buff[jj];
        tmp_cells[addr].speeds[7] = recv_buff[jj+1];
        tmp_cells[addr].speeds[8] = recv_buff[jj+2];
    }

     if (params.my_rank==1)
    {
        for (jj = 0; jj < params.local_ncols/*-96*/; ++jj)
        {
            printf("South_Source:Cell:%d\ns_W=%f s=%f s_E=%f\n\n\n\n\n"
               , jj ,tmp_cells[jj+(params.local_ncols*params.local_nrows)].speeds[7]
               , tmp_cells[jj+(params.local_ncols*params.local_nrows)].speeds[4]
               , tmp_cells[jj+(params.local_ncols*params.local_nrows)].speeds[8]);
        }
    }

    /* loop over _middle_ cells */
    for (ii = 2; ii < params.local_nrows; ii++)
    {
        addr = ii*params.local_ncols;
        for (jj = 0; jj < params.local_ncols; jj++, addr++)
        { 
            /*  Add onstacle consider*/
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            // jj=127; ii=25;
            // jj=0; ii=20;
            // jj=0; ii=1;
            // jj=0; ii=25;
            // jj=0; ii=4;
            y_n = (ii + 1) % (params.local_nrows+2);
            x_e = (jj + 1) % (params.nx);
            y_s = /*(ii == 0) ? (ii + params.local_nrows - 1) :*/ (ii - 1);
            x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);

            // if (params.my_rank==1)
            // {
            //   printf("Rank: %d\nCell:%d ii=%d, jj=%d,\n\ty_n=%d\nx_w=%d\t\tx_e=%d\n\ty_s=%d\n\n",params.my_rank,jj+(ii)*params.local_ncols,ii,jj
            //     ,y_n*params.local_ncols + jj,ii *params.local_ncols + x_w,ii *params.local_ncols + x_e,y_s*params.local_ncols + jj);
            // }
            // MPI_Barrier(MPI_COMM_WORLD);
            // MPI_Abort(MPI_COMM_WORLD, 0);

            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii*params.local_ncols + jj].speeds[0]  = cells[addr].speeds[0]; /* central cell, */
                                                     /* no movement   */
            tmp_cells[ii *params.local_ncols + x_e].speeds[1] = cells[addr].speeds[1]; /* east */
            tmp_cells[y_n*params.local_ncols + jj].speeds[2]  = cells[addr].speeds[2]; /* north */
            tmp_cells[ii *params.local_ncols + x_w].speeds[3] = cells[addr].speeds[3]; /* west */
            tmp_cells[y_s*params.local_ncols + jj].speeds[4]  = cells[addr].speeds[4]; /* south */
            tmp_cells[y_n*params.local_ncols + x_e].speeds[5] = cells[addr].speeds[5]; /* north-east */
            tmp_cells[y_n*params.local_ncols + x_w].speeds[6] = cells[addr].speeds[6]; /* north-west */
            tmp_cells[y_s*params.local_ncols + x_w].speeds[7] = cells[addr].speeds[7]; /* south-west */
            tmp_cells[y_s*params.local_ncols + x_e].speeds[8] = cells[addr].speeds[8]; /* south-east */
        }
    }

    addr = (params.local_nrows)*params.local_ncols;
    ii = params.local_nrows;
    for (jj = 0; jj < params.local_ncols; jj++, addr++)
    { 
        x_e = (jj + 1) % (params.nx);
        y_s = /*(ii == 0) ? (ii + params.local_nrows - 1) :*/ (ii - 1);
        x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);
        /* propagate densities to neighbouring cells, following
        ** appropriate directions of travel and writing into
        ** scratch space grid */
        tmp_cells[ii*params.local_ncols + jj].speeds[0]  = cells[addr].speeds[0]; /* central cell, */
                                                 /* no movement   */
        tmp_cells[ii *params.local_ncols + x_e].speeds[1] = cells[addr].speeds[1]; /* east */
        send_buff[jj]  = cells[addr].speeds[2]; /* north */
        tmp_cells[ii *params.local_ncols + x_w].speeds[3] = cells[addr].speeds[3]; /* west */
        tmp_cells[y_s*params.local_ncols + jj].speeds[4]  = cells[addr].speeds[4]; /* south */
        send_buff[jj+1] = cells[addr].speeds[5]; /* north-east */
        send_buff[jj+2] = cells[addr].speeds[6]; /* north-west */
        tmp_cells[y_s*params.local_ncols + x_w].speeds[7] = cells[addr].speeds[7]; /* south-west */
        tmp_cells[y_s*params.local_ncols + x_e].speeds[8] = cells[addr].speeds[8]; /* south-east */
    }

    MPI_Sendrecv(send_buff, params.local_ncols*3, MPI_DOUBLE, params.next, HALO_VELS,
                 recv_buff, params.local_ncols*3, MPI_DOUBLE, params.prev, HALO_VELS,
                   MPI_COMM_WORLD, &status);

    addr = params.local_ncols;
    for (jj = 0; jj < params.local_ncols; jj++, addr++)
    {
        tmp_cells[addr].speeds[2] = recv_buff[jj];
        tmp_cells[addr].speeds[5] = recv_buff[jj+1];
        tmp_cells[addr].speeds[6] = recv_buff[jj+2];
    }
}

void accelerate_flow_Row_RW(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles)
{
    int jj;     /* generic counters */
    double w1,w2;  /* weighting factors */

    /* compute weighting factors */
    w1 = params.density * params.accel / 9.0;
    w2 = params.density * params.accel / 36.0;

    int addr = (accel_area.idx)*params.local_ncols;
    // printf("row=%d, col=%d,  start at index=%d\n", ii, params.local_ncols, addr);

    // The new index accounts for the halo row
    for (jj = 0; jj < params.local_ncols; jj++, addr++)
    {
        /* if the cell is not occupied and
        ** we don't send a density negative */
        if (!obstacles[addr] &&
        (cells[addr].speeds[3] - w1) > 0.0 &&
        (cells[addr].speeds[6] - w2) > 0.0 &&
        (cells[addr].speeds[7] - w2) > 0.0 )
        {
            /* increase 'east-side' densities */
            cells[addr].speeds[1] += w1;
            cells[addr].speeds[5] += w2;
            cells[addr].speeds[8] += w2;
            /* decrease 'west-side' densities */
            cells[addr].speeds[3] -= w1;
            cells[addr].speeds[6] -= w2;
            cells[addr].speeds[7] -= w2;
        }
    }
}

void collision_local(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj,kk, addr;                 /* generic counters */
    const double c_sq = 1.0/3.0;  /* square of speed of sound */
    const double w0 = 4.0/9.0;    /* weighting factor */
    const double w1 = 1.0/9.0;    /* weighting factor */
    const double w2 = 1.0/36.0;   /* weighting factor */
    const double tmp1 = 2.0 * c_sq;
    const double tmp2 = tmp1 * c_sq;

    double u_x,u_y;               /* av. velocities in x and y directions */
    double u_sq;                  /* squared velocity */
    double local_density;         /* sum of densities in a particular cell */
    double u[NSPEEDS];            /* directional velocities */
    double d_equ[NSPEEDS];        /* equilibrium densities */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (ii = 1; ii <=  params.local_nrows; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            addr = ii*params.nx + jj;
            
            /* don't consider occupied cells */
            if (obstacles[ii*params.nx + jj])
            {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[addr].speeds[1] = tmp_cells[addr].speeds[3];
                cells[addr].speeds[2] = tmp_cells[addr].speeds[4];
                cells[addr].speeds[3] = tmp_cells[addr].speeds[1];
                cells[addr].speeds[4] = tmp_cells[addr].speeds[2];
                cells[addr].speeds[5] = tmp_cells[addr].speeds[7];
                cells[addr].speeds[6] = tmp_cells[addr].speeds[8];
                cells[addr].speeds[7] = tmp_cells[addr].speeds[5];
                cells[addr].speeds[8] = tmp_cells[addr].speeds[6];
            }
            else
            {
                /* compute local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += tmp_cells[addr].speeds[kk];
                }
// 
                // if (local_density == 0.0)
                // {
                   // u[1] =   0.0;        /* east */
                   // u[2] =         0.0;  /* north */
                   // u[3] =  0.0;        /* west */
                   // u[4] =       0.0;  /* south */
                   // u[5] =   0.0 ;  /* north-east */
                   // u[6] =  0.0 ;  /* north-west */
                   // u[7] =  0.0 ;  /* south-west */
                   // u[8] =   0.0 ;  /* south-east */
                   // printf("This has Run in collision\n");
                // }     
                // else
                // {
                /*compute x velocity component*/ 
                u_x = (tmp_cells[addr].speeds[1] +
                        tmp_cells[addr].speeds[5] +
                        tmp_cells[addr].speeds[8]
                    - (tmp_cells[addr].speeds[3] +
                        tmp_cells[addr].speeds[6] +
                        tmp_cells[addr].speeds[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (tmp_cells[addr].speeds[2] +
                        tmp_cells[addr].speeds[5] +
                        tmp_cells[addr].speeds[6]
                    - (tmp_cells[addr].speeds[4] +
                        tmp_cells[addr].speeds[7] +
                        tmp_cells[addr].speeds[8]))
                    / local_density;

                /* velocity squared */
                u_sq = u_x * u_x + u_y * u_y;

                /* directional velocity components */
                u[1] =   u_x;        /* east */
                u[2] =         u_y;  /* north */
                u[3] = - u_x;        /* west */
                u[4] =       - u_y;  /* south */
                u[5] =   u_x + u_y;  /* north-east */
                u[6] = - u_x + u_y;  /* north-west */
                u[7] = - u_x - u_y;  /* south-west */
                u[8] =   u_x - u_y;  /* south-east */
                // }

                /* equilibrium densities */
                /* zero velocity density: weight w0 */
                d_equ[0] = w0 * local_density * (1.0 - u_sq / (tmp1));
                /* axis speeds: weight w1 */
                d_equ[1] = w1 * local_density * (1.0 + u[1] / c_sq
                    + (u[1] * u[1]) / (tmp2)
                    - u_sq / (tmp1));
                d_equ[2] = w1 * local_density * (1.0 + u[2] / c_sq
                    + (u[2] * u[2]) / (tmp2)
                    - u_sq / (tmp1));
                d_equ[3] = w1 * local_density * (1.0 + u[3] / c_sq
                    + (u[3] * u[3]) / (tmp2)
                    - u_sq / (tmp1));
                d_equ[4] = w1 * local_density * (1.0 + u[4] / c_sq
                    + (u[4] * u[4]) / (tmp2)
                    - u_sq / (tmp1));
                /* diagonal speeds: weight w2 */
                d_equ[5] = w2 * local_density * (1.0 + u[5] / c_sq
                    + (u[5] * u[5]) / (tmp2)
                    - u_sq / (tmp1));
                d_equ[6] = w2 * local_density * (1.0 + u[6] / c_sq
                    + (u[6] * u[6]) / (tmp2)
                    - u_sq / (tmp1));
                d_equ[7] = w2 * local_density * (1.0 + u[7] / c_sq
                    + (u[7] * u[7]) / (tmp2)
                    - u_sq / (tmp1));
                d_equ[8] = w2 * local_density * (1.0 + u[8] / c_sq
                    + (u[8] * u[8]) / (tmp2)
                    - u_sq / (tmp1));

                /* relaxation step */
                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    cells[addr].speeds[kk] = 
                        (tmp_cells[addr].speeds[kk] + params.omega * 
                        (d_equ[kk] - tmp_cells[addr].speeds[kk]));
                }
            }
        }
    }
}

void accelerate_flow_Colum_RW(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles)
{
    int ii;     /* generic counters */
    double w1,w2;  /* weighting factors */

    /* compute weighting factors */
    w1 = params.density * params.accel / 9.0;
    w2 = params.density * params.accel / 36.0;

    // for (ii = 0; ii < params.local_nrows+1; ii++)
    for (ii = 1; ii <= params.local_nrows; ii++)
    {
        int addr = ii*params.local_ncols + accel_area.idx;
        /* if the cell is not occupied and
        ** we don't send a density negative */
        if (!obstacles[addr] &&
        (cells[addr].speeds[4] - w1) > 0.0 &&
        (cells[addr].speeds[7] - w2) > 0.0 &&
        (cells[addr].speeds[8] - w2) > 0.0 )
        {
            /* increase 'north-side' densities */
            cells[addr].speeds[2] += w1;
            cells[addr].speeds[5] += w2;
            cells[addr].speeds[6] += w2;
            /* decrease 'south-side' densities */
            cells[addr].speeds[4] -= w1;
            cells[addr].speeds[7] -= w2;
            cells[addr].speeds[8] -= w2;
        }
    }

}

// for buffer-less, exchange on the grid
void propagate_row_wise(const param_t params, speed_t* cells, speed_t* tmp_cells)
{
    int ii,jj,addr;            /* generic counters */

    /* loop over _all_ cells */
    for (ii = 1; ii <= params.local_nrows; ii++)
    {
        addr = ii*params.local_ncols;
        for (jj = 0; jj < params.local_ncols; jj++, addr++)
        { 
            /*  Add onstacle consider*/
            int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            // jj=127; ii=25;
            // jj=0; ii=20;
            // jj=0; ii=1;
            // jj=0; ii=25;
            // jj=0; ii=4;
            y_n = (ii + 1) % (params.local_nrows+2);
            x_e = (jj + 1) % (params.nx);
            y_s = /*(ii == 0) ? (ii + params.local_nrows - 1) :*/ (ii - 1);
            x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);

            // if (params.my_rank==1)
            // {
            //   printf("Rank: %d\nCell:%d ii=%d, jj=%d,\n\ty_n=%d\nx_w=%d\t\tx_e=%d\n\ty_s=%d\n\n",params.my_rank,jj+(ii)*params.local_ncols,ii,jj
            //     ,y_n*params.local_ncols + jj,ii *params.local_ncols + x_w,ii *params.local_ncols + x_e,y_s*params.local_ncols + jj);
            // }
            // MPI_Barrier(MPI_COMM_WORLD);
            // MPI_Abort(MPI_COMM_WORLD, 0);

            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii*params.local_ncols + jj].speeds[0]  = cells[addr].speeds[0]; /* central cell, */
                                                     /* no movement   */
            tmp_cells[ii *params.local_ncols + x_e].speeds[1] = cells[addr].speeds[1]; /* east */
            tmp_cells[y_n*params.local_ncols + jj].speeds[2]  = cells[addr].speeds[2]; /* north */
            tmp_cells[ii *params.local_ncols + x_w].speeds[3] = cells[addr].speeds[3]; /* west */
            tmp_cells[y_s*params.local_ncols + jj].speeds[4]  = cells[addr].speeds[4]; /* south */
            tmp_cells[y_n*params.local_ncols + x_e].speeds[5] = cells[addr].speeds[5]; /* north-east */
            tmp_cells[y_n*params.local_ncols + x_w].speeds[6] = cells[addr].speeds[6]; /* north-west */
            tmp_cells[y_s*params.local_ncols + x_w].speeds[7] = cells[addr].speeds[7]; /* south-west */
            tmp_cells[y_s*params.local_ncols + x_e].speeds[8] = cells[addr].speeds[8]; /* south-east */
        }
    }
}