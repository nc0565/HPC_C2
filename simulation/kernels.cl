#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define NSPEEDS 9

/* struct to hold the 'speed' values */
typedef struct {
    double speeds[NSPEEDS];
} speed_t;

/* struct to hold the parameter values */
typedef struct {
    int nx;              /* no. of cells in x-direction */
    int ny;              /* no. of cells in y-direction */
    int max_iters;       /* no. of iterations */
    int reynolds_dim;    /* dimension for Reynolds number */
    double density;      /* density per link */
    double accel;        /* density redistribution */
    double omega;        /* relaxation parameter */
} param_t;

typedef enum { ACCEL_ROW=0, ACCEL_COLUMN=1 } accel_e;
typedef struct {
    int col_or_row;
    int idx;
} accel_area_t;

/*
*   TODO
*   Write OpenCL kernels
*/

__kernel
void collision(param_t params, __global speed_t* cells, __global speed_t* tmp_cells, __global int* obstacles)
{
     double c_sq = (1.0/3.0);  /* square of speed of sound */
     double tmp1 = (2.0 * c_sq);
     double tmp2 = (tmp1 * c_sq);
     double w0 = (4.0/9.0);    /* weighting factor */
     double w1 = (1.0/9.0);    /* weighting factor */
     double w2 = (1.0/36.0);   /* weighting factor */

    double u_x,u_y;               /* av. velocities in x and y directions */
    double u_sq;                  /* squared velocity */
    double local_density;         /* sum of densities in a particular cell */
    double u[NSPEEDS];            /* directional velocities */
    double d_equ[NSPEEDS];        /* equilibrium densities */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    int addr = get_global_id(0)*params.nx + get_global_id(1);

    /*
    // Select to remove obsticle condition.
    long test = (obstacles[addr]);
    cells[addr].speeds[1] = select(cells[addr].speeds[1], tmp_cells[addr].speeds[3], test);
    cells[addr].speeds[2] = select(cells[addr].speeds[2], tmp_cells[addr].speeds[4], test);
    cells[addr].speeds[3] = select(cells[addr].speeds[3], tmp_cells[addr].speeds[1], test);
    cells[addr].speeds[4] = select(cells[addr].speeds[4], tmp_cells[addr].speeds[2], test);
    cells[addr].speeds[5] = select(cells[addr].speeds[5], tmp_cells[addr].speeds[7], test);
    cells[addr].speeds[6] = select(cells[addr].speeds[6], tmp_cells[addr].speeds[8], test);
    cells[addr].speeds[7] = select(cells[addr].speeds[7], tmp_cells[addr].speeds[5], test);
    cells[addr].speeds[8] = select(cells[addr].speeds[8], tmp_cells[addr].speeds[6], test);
    */
    /* called after propagate, so taking values from scratch space
    ** mirroring, and writing into main grid */

    if (obstacles[addr])
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
    /* don't consider occupied cells */
    else
    {


    /* compute local density total */
    local_density = 0.0;

    local_density += tmp_cells[addr].speeds[0];
    local_density += tmp_cells[addr].speeds[1];
    local_density += tmp_cells[addr].speeds[2];
    local_density += tmp_cells[addr].speeds[3];
    local_density += tmp_cells[addr].speeds[4];
    local_density += tmp_cells[addr].speeds[5];
    local_density += tmp_cells[addr].speeds[6];
    local_density += tmp_cells[addr].speeds[7];
    local_density += tmp_cells[addr].speeds[8];


    /* compute x velocity component */
    u_x = (tmp_cells[addr].speeds[1] + tmp_cells[addr].speeds[5] + tmp_cells[addr].speeds[8]   \
        - (tmp_cells[addr].speeds[3] + tmp_cells[addr].speeds[6] + tmp_cells[addr].speeds[7])) \
        / local_density;

    /* compute y velocity component */
    u_y = (tmp_cells[addr].speeds[2] + tmp_cells[addr].speeds[5] + tmp_cells[addr].speeds[6]   \
        - (tmp_cells[addr].speeds[4] + tmp_cells[addr].speeds[7] + tmp_cells[addr].speeds[8])) \
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
    d_equ[0] = w0 * local_density * (1.0 - u_sq / tmp1);
    /* axis speeds: weight w1 */
    d_equ[1] = w1 * local_density * (1.0 + u[1] / c_sq + (u[1] * u[1]) / tmp2 - u_sq / (tmp1));
    d_equ[2] = w1 * local_density * (1.0 + u[2] / c_sq + (u[2] * u[2]) / tmp2 - u_sq / (tmp1));
    d_equ[3] = w1 * local_density * (1.0 + u[3] / c_sq + (u[3] * u[3]) / tmp2 - u_sq / (tmp1));
    d_equ[4] = w1 * local_density * (1.0 + u[4] / c_sq + (u[4] * u[4]) / tmp2 - u_sq / (tmp1));
    /* diagonal speeds: weight w2 */
    d_equ[5] = w2 * local_density * (1.0 + u[5] / c_sq + (u[5] * u[5]) / tmp2 - u_sq / (tmp1));
    d_equ[6] = w2 * local_density * (1.0 + u[6] / c_sq + (u[6] * u[6]) / tmp2 - u_sq / (tmp1));
    d_equ[7] = w2 * local_density * (1.0 + u[7] / c_sq + (u[7] * u[7]) / tmp2 - u_sq / (tmp1));
    d_equ[8] = w2 * local_density * (1.0 + u[8] / c_sq + (u[8] * u[8]) / tmp2 - u_sq / tmp1);

    /* relaxation step */
    
    cells[addr].speeds[0] = (tmp_cells[addr].speeds[0] + params.omega * (d_equ[0] - tmp_cells[addr].speeds[0]));
    cells[addr].speeds[1] = (tmp_cells[addr].speeds[1] + params.omega * (d_equ[1] - tmp_cells[addr].speeds[1]));
    cells[addr].speeds[2] = (tmp_cells[addr].speeds[2] + params.omega * (d_equ[2] - tmp_cells[addr].speeds[2]));
    cells[addr].speeds[3] = (tmp_cells[addr].speeds[3] + params.omega * (d_equ[3] - tmp_cells[addr].speeds[3]));
    cells[addr].speeds[4] = (tmp_cells[addr].speeds[4] + params.omega * (d_equ[4] - tmp_cells[addr].speeds[4]));
    cells[addr].speeds[5] = (tmp_cells[addr].speeds[5] + params.omega * (d_equ[5] - tmp_cells[addr].speeds[5]));
    cells[addr].speeds[6] = (tmp_cells[addr].speeds[6] + params.omega * (d_equ[6] - tmp_cells[addr].speeds[6]));
    cells[addr].speeds[7] = (tmp_cells[addr].speeds[7] + params.omega * (d_equ[7] - tmp_cells[addr].speeds[7]));
    cells[addr].speeds[8] = (tmp_cells[addr].speeds[8] + params.omega * (d_equ[8] - tmp_cells[addr].speeds[8]));

    }
    /*
    cells[addr].speeds[0] = select((tmp_cells[addr].speeds[0] + params.omega * (d_equ[0] - tmp_cells[addr].speeds[0])), cells[addr].speeds[0], test);
    cells[addr].speeds[1] = select((tmp_cells[addr].speeds[1] + params.omega * (d_equ[1] - tmp_cells[addr].speeds[1])), cells[addr].speeds[1], test);
    cells[addr].speeds[2] = select((tmp_cells[addr].speeds[2] + params.omega * (d_equ[2] - tmp_cells[addr].speeds[2])), cells[addr].speeds[2], test);
    cells[addr].speeds[3] = select((tmp_cells[addr].speeds[3] + params.omega * (d_equ[3] - tmp_cells[addr].speeds[3])), cells[addr].speeds[3], test);
    cells[addr].speeds[4] = select((tmp_cells[addr].speeds[4] + params.omega * (d_equ[4] - tmp_cells[addr].speeds[4])), cells[addr].speeds[4], test);
    cells[addr].speeds[5] = select((tmp_cells[addr].speeds[5] + params.omega * (d_equ[5] - tmp_cells[addr].speeds[5])), cells[addr].speeds[5], test);
    cells[addr].speeds[6] = select((tmp_cells[addr].speeds[6] + params.omega * (d_equ[6] - tmp_cells[addr].speeds[6])), cells[addr].speeds[6], test);
    cells[addr].speeds[7] = select((tmp_cells[addr].speeds[7] + params.omega * (d_equ[7] - tmp_cells[addr].speeds[7])), cells[addr].speeds[7], test);
    cells[addr].speeds[8] = select((tmp_cells[addr].speeds[8] + params.omega * (d_equ[8] - tmp_cells[addr].speeds[8])), cells[addr].speeds[8], test);
    */
}


__kernel
void prop(param_t params, __global speed_t* cells, __global speed_t* tmp_cells)
{
    int ii=get_global_id(0),jj=get_global_id(1),nx=params.nx, ny=params.ny;            /* generic counters */

    int addr = ii*nx + jj;

    int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
    /* determine indices of axis-direction neighbours
    ** respecting periodic boundary conditions (wrap around) */
    y_n = (ii + 1) % ny;
    x_e = (jj + 1) % nx;
    y_s = (ii != 0)? (ii - 1): (ii + ny - 1);
    x_w = (jj != 0)? (jj - 1): (jj + nx - 1);
    /* propagate densities to neighbouring cells, following
    ** appropriate directions of travel and writing into
    ** scratch space grid */
    tmp_cells[addr].speeds[0]  = cells[addr].speeds[0]; /* central cell, */
                                             /* no movement   */
    tmp_cells[ii*nx + x_e].speeds[1] = cells[addr].speeds[1]; /* east */
    tmp_cells[y_n*nx + jj].speeds[2]  = cells[addr].speeds[2]; /* north */
    tmp_cells[ii*nx + x_w].speeds[3] = cells[addr].speeds[3]; /* west */
    tmp_cells[y_s*nx + jj].speeds[4]  = cells[addr].speeds[4]; /* south */
    tmp_cells[y_n*nx + x_e].speeds[5] = cells[addr].speeds[5]; /* north-east */
    tmp_cells[y_n*nx + x_w].speeds[6] = cells[addr].speeds[6]; /* north-west */
    tmp_cells[y_s*nx + x_w].speeds[7] = cells[addr].speeds[7]; /* south-west */
    tmp_cells[y_s*nx + x_e].speeds[8] = cells[addr].speeds[8]; /* south-east */
}

__kernel
void accel(param_t params, accel_area_t accel_area, __global speed_t* cells, __global int* obstacles)
{
    /* compute weighting factors */
    double w1 = params.density * params.accel / 9.0, w2 = params.density * params.accel / 36.0;

    if (accel_area.col_or_row == ACCEL_COLUMN)
    {
        //int addr = get_global_id(0)*params.nx + accel_area.idx;

        if (get_global_id(1) != accel_area.idx) return;
        int addr = get_global_id(0)*params.nx + accel_area.idx;

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
    else
    {
        //int addr = accel_area.idx*params.nx + get_global_id(1);

        if (get_global_id(0) != accel_area.idx) return;
        int addr = accel_area.idx*params.nx + get_global_id(1);

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

__kernel
void av_vel(param_t params, __global speed_t* cells, __global int* obstacles, __global double *av_output)
{
    //int    tot_cells = 0;  /* no. of cells used in calculation */
    /* initialise */
    //double tot_u = 0.0;    /* accumulated magnitudes of velocity for each cell */


    /* loop over all non-blocked cells */
    int addr = get_global_id(0)*params.nx+get_global_id(1);
    /* ignore occupied cells */
    if (!obstacles[addr])
    {
        /* local density total */
        double local_density = 0.0;

        local_density += cells[addr].speeds[1];
        local_density += cells[addr].speeds[2];
        local_density += cells[addr].speeds[3];
        local_density += cells[addr].speeds[4];
        local_density += cells[addr].speeds[5];
        local_density += cells[addr].speeds[6];
        local_density += cells[addr].speeds[7];
        local_density += cells[addr].speeds[8];

        /* x-component of velocity */
        double u_x = (cells[addr].speeds[1]+cells[addr].speeds[5]+cells[addr].speeds[8]
            - (cells[addr].speeds[3]+cells[addr].speeds[6]+cells[addr].speeds[7])) /
            local_density;

        /* compute y velocity component */
        double u_y = (cells[addr].speeds[2]+cells[addr].speeds[5]+cells[addr].speeds[6]
            - (cells[addr].speeds[4]+cells[addr].speeds[7]+cells[addr].speeds[8])) /
            local_density;

        /* accumulate the norm of x- and y- velocity components */
        //tot_u = sqrt(u_x*u_x + u_y*u_y);
        /* increase counter of inspected cells */
        //++tot_cells;
        av_output[get_group_id(0)] += sqrt(u_x*u_x + u_y*u_y);
    }

}

/*  10 tips
*    Unroll loops: Comparison operations are expensive – if you know how many iterations you need simply perform them individually.
*    Disable processing of denormalised numbers: The processing of denormalised numbers takes time. Disable them using the -cl-denorms-are-zero option. You can also disable processing of infinite values and NaNs using -cl-finite-math-only.
*    Transfer constant primitive values to the kernel with compiler defines instead of private memory parameters. Example: clBuildProgram(program, 0, NULL, "-DSIZE=128", NULL, NULL);.
*    Store small variable values in private memory instead of local memory: Private memory is faster and should be used unless you need to share data between work items in which case use local memory.
*    Avoid local memory bank conflicts by accessing local memory sequentially: Successive 32-bit elements are stored in successive local memory banks. Sequential bank access is parallel whereas contending on same bank access is serial.
*    Avoid using modulo operator: This operator is slow and alternatives should be used instead.
*    Reuse private variables throughout the kernel but use macros to distinguish each separate use.
*    For multiply-and-add operations use the fma function if available.
*    Inline non-kernel functions in a program file: This uses more memory but removes the context switches and stack operations associated with regular function calls.
*    Avoid branch miss penalties by coding conditional statements to be true more often than false: Processors often optimise for the true case and this can result in a penalty for the false case known as the branch miss penalty. Code your conditionals to evaluate to true as often as possible.
*/

