#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

//Custom
#include "define_env.h"
#include "orbits.h"
#include "richardson.h"
#include "custom.h"
#include "differential_correction.h"
#include "halo_orbit_computation.h"
#include "energy.h"
#include "manifold.h"


//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>
#include <gsl_blas.h>

//Gnuplot
#include "gnuplot_i.h"



void init_manifold_branch(Manifold_branch *branch, Orbit orbit, int stability, int way, struct value_params val_par)
{
    int i;
    branch->orbit = orbit;
    branch->stability = stability;
    branch->way = way;
    branch->val_par = val_par;
    branch->final_time = 0.0;
    branch->direction = gsl_vector_calloc(6);  //all the elements to zero
    branch->events_mat = (double**) malloc(6 * sizeof(double));
    for (i = 0; i < 6; i++)
        branch->events_mat[i]= (double*) malloc((val_par.max_events+1) * sizeof(double));  //+1 to allowed the last state to be stored

    branch->fvalue.val_par = &val_par;
    branch->fvalue.value = &linear_intersection;
}

void free_manifold_branche(Manifold_branch branch)
{
    free(branch.direction);
    free(branch.events_mat);
}


int manifold_branch_computation(Manifold_branch *branch,
                                double theta,
                                double tf,
                                custom_ode_structure *ode_s,
                                custom_ode_structure *ode_s_6)

{
    //------------------------------------------------
    //Init
    //------------------------------------------------
    int i;
    int dim = ode_s->d->sys->dimension;
    int dim_6 = ode_s_6->d->sys->dimension;
    double yv[dim];
    double t;

    //------------------------------------------------
    //Integration until theta is reached on the orbit
    //------------------------------------------------
    for(i=0;i<dim;i++) yv[i] = branch->orbit.y0[i];


    t = 0.0;
    if(theta >0) gsl_odeiv2_driver_apply(ode_s->d, &t, theta*branch->orbit.T, yv);


    //------------------------------------------------
    //Recovering the STM & vectors
    //------------------------------------------------
    gsl_matrix *STM = gsl_matrix_alloc(6,6);
    custom_vectorToMatrix(STM,yv,6, 6, 6);


    //Product with stable or unstable gsl vector
    switch(branch->stability)
    {
        case 1: //stable manifold
        gsl_blas_dgemv (CblasNoTrans , 1.0 , STM , branch->orbit.stable_direction , 0.0 , branch->direction);
        break;

        case -1: //unstable manifold
        gsl_blas_dgemv (CblasNoTrans , 1.0 , STM , branch->orbit.unstable_direction , 0.0 , branch->direction);
        break;
    }
    //Normalization
    gsl_vector_scale (branch->direction , 1/gsl_blas_dnrm2(branch->direction));

    //------------------------------------------------
    //Integration
    //------------------------------------------------
    //Scaling of direction vector
    gsl_vector_scale (branch->direction , branch->way*gsl_vector_get(branch->direction,0)*branch->orbit.cr3bp.d_man);

    //Initial state
    double y0b[dim_6];
    t = 0.0;
    for(i=0;i<dim_6;i++) y0b[i] = yv[i] + gsl_vector_get(branch->direction,i);

    //Initial position storage
    for(i=0;i<dim_6;i++) branch->initial_position[i] = y0b[i];

    //Backwards or forward integration?
    if(branch->stability==1) tf*=-1;

    //Termination scheme & integration
    int last_indix;
    //Starting in the right direction
    ode_s_6->d->h = (tf>t) ? fabs(ode_s_6->d->h) : -fabs(ode_s_6->d->h);
    ode_s_6->h = (tf>t) ? fabs(ode_s_6->h) : -fabs(ode_s_6->h);
    switch(branch->val_par.dim)
    {
        case -1: // "free" integration until the time tf is reached
        //Integration
        gsl_odeiv2_driver_apply(ode_s_6->d, &t, tf, y0b);
        //Output storage
        for(i=0;i<dim_6;i++) branch->final_position[0] = y0b[i];
        break;

        default: //termination by event triggering is desired
        //Integration
        last_indix = custom_odezero_2(y0b, branch->events_mat, &t, tf, ode_s_6, branch->fvalue);
        //Output storage
        for(i=0;i<dim_6;i++) branch->final_position[i] = branch->events_mat[i][last_indix];
        break;
    }

    //Final time
    branch->final_time = t;

    //Reset
    reset_ode_structure(ode_s);
    reset_ode_structure(ode_s_6);

    return GSL_SUCCESS;
}



void manifold_branch_plot(Manifold_branch *branch, gnuplot_ctrl *h1, int type, int points, int scale, custom_ode_structure *ode_s_6)
{
    //Checking initialization of h1
    if(h1 == NULL) h1 = gnuplot_init();

    //Init
    double t, t1, ti;
    double x1[points+1];
    double x2[points+1];
    double y[6];
    int i;
    for(i=0; i<6; i++) y[i]  = branch->initial_position[i];
    double scale_factor;

    //Scale
    switch(scale)
    {
    case ADIM:
        scale_factor = 1;
        break;
    case DIM:
        scale_factor = branch->orbit.cr3bp.L;
        break;
    default: //default: ADIM
        scale_factor = 1;
        break;
    }

    //First value
    switch(type)
    {
    case XY:
        x1[0] = y[0]*scale_factor;
        x2[0] = y[1]*scale_factor;
        break;
    case XZ:
        x1[0] = y[0]*scale_factor;
        x2[0] = y[2]*scale_factor;
        break;
    case YZ:
        x1[0] = y[1]*scale_factor;
        x2[0] = y[2]*scale_factor;
        break;
    }

    //Integration over a period
    t = 0;
    t1 = branch->final_time;
    //Starting in the right direction
    ode_s_6->d->h = (t1>t) ? fabs(ode_s_6->d->h) : -fabs(ode_s_6->d->h);
    for(i =1; i<=points; i++)
    {
        ti = i * t1 / points;
        gsl_odeiv2_driver_apply (ode_s_6->d, &t, ti, y);

        switch(type)
        {
        case XY:
            x1[i] = y[0]*scale_factor;
            x2[i] = y[1]*scale_factor;
            break;
        case XZ:
            x1[i] = y[0]*scale_factor;
            x2[i] = y[2]*scale_factor;
            break;
        case YZ:
            x1[i] = y[1]*scale_factor;
            x2[i] = y[2]*scale_factor;
            break;
        }
    }

    //Plotting
    switch(type)
    {
    case XY:
        gnuplot_set_xlabel (h1,"X");
        gnuplot_set_ylabel (h1,"Y");
        break;
    case XZ:
        gnuplot_set_xlabel (h1,"X");
        gnuplot_set_ylabel (h1,"Z");
        break;
    case YZ:
        gnuplot_set_xlabel (h1,"Y");
        gnuplot_set_ylabel (h1,"Z");
        break;
    }
    gnuplot_setstyle(h1, "lines");
    gnuplot_setcolor(h1, 2);
    gnuplot_plot_xy(h1, x1, x2, points+1, "manifold");

    //Reset
    reset_ode_structure(ode_s_6);


}
