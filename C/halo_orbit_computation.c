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

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>

//Gnuplot
#include "gnuplot_i.h"



int halo_orbit_computation(Orbit *halo,
                           custom_ode_structure *ode_s,
                           diff_corr_parameters diff_corr)
{
    int status;
    //-----------------------------------------------------------------------------------------------------
    // Determination of a good starting point for the given parameters Az, m and Li
    //-----------------------------------------------------------------------------------------------------
    //Initialization of the initial guess
    double ystart[42];

    //Third order approximation of the state
    third_order_halo_orbit(*halo, 0.0, ystart);

    //STM(0) is concatenated to the state
    gsl_matrix *STM0 = gsl_matrix_alloc(6,6);
    gsl_matrix_set_identity (STM0);
    custom_matrixToVector(ystart, STM0, 6, 6, 6);
    gsl_matrix_free(STM0);

    //-----------------------------------------------------------------------------------------------------
    // Refinement
    //-----------------------------------------------------------------------------------------------------
    status = halo_orbit_refinement(halo, ystart, ode_s, diff_corr);


    if(status != GSL_SUCCESS)
    {
        printf("WARNING: halo_orbit_refinement failed to converge in halo_orbit_computation. Premature ending.\n");
        return GSL_FAILURE;
    }

    //Updating halo status
    halo->status = HALO_REAL;

    return GSL_SUCCESS;

}


int halo_orbit_refinement(Orbit *halo,
                          double ystart[],
                          custom_ode_structure *ode_s,
                          diff_corr_parameters diff_corr)
{


    //-----------------------------------------------------------------------------------------------------
    // Differential correction
    //-----------------------------------------------------------------------------------------------------
    int i, status;
    //STM init
    gsl_matrix *STM0 = gsl_matrix_alloc(6,6);
    gsl_matrix_set_identity (STM0);
    //Copy of the initial state for next step
    double y[42], yhalf[42];
    for(i=0; i<42; i++) y[i] = ystart[i];

    //Maximum integration time (never reached)
    double t1 = 10;
    //Initialization of crossing time;
    double tcross = 0;

    //Differential correction procedure
    status = differential_correction(halo->cr3bp, y, yhalf, &tcross, t1, ode_s, diff_corr);

    if(status != GSL_SUCCESS)
    {
        printf("WARNING: differential_correction failed to converge in halo_orbit_refinement. Premature ending.\n");
        return GSL_FAILURE;
    }

    //Storage of good initial point
    for (i=0; i<6; i++)halo->y0[i] =y[i];
    custom_matrixToVector(halo->y0, STM0, 6, 6, 6);
    //Storage of half-period point
    for (i=0; i<6; i++)halo->yhalf[i]=yhalf[i];
    custom_matrixToVector(halo->yhalf, STM0, 6, 6, 6);


    //Storage of orbital period
    halo->T = 2*tcross;

    //Az
    halo->Az_T = (fabs(y[2]) > fabs(yhalf[2])) ? y[2] : yhalf[2];
    halo->Azdim_T = halo->cr3bp.L * (halo->Az_T);

    //-----------------------------------------------------------------------------------------------------
    // Computation of the monodromy matrix and corresponding eigenvectors
    //-----------------------------------------------------------------------------------------------------
    //STM concatenation
    custom_matrixToVector(y, STM0, 6, 6, 6);

    //Integration over a period
    double t = 0;
    gsl_odeiv2_driver_reset(ode_s->d);
    gsl_odeiv2_driver_apply (ode_s->d, &t, 2*tcross, y);

    //Monodromy matrix
    gsl_matrix *monodromy = gsl_matrix_alloc(6,6);
    custom_vectorToMatrix(monodromy,y,6, 6, 6);

    //Output
    halo->monodromy = monodromy;

    //Eigenvectors
    gsl_vector_complex *eval = gsl_vector_complex_alloc (6);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (6, 6);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (6);
    gsl_eigen_nonsymmv (monodromy, eval, evec, w);


    //Sort in descending order
    gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    //printf("Eigenvectors:\n");
    //for(i = 0;i<6;i++) printf("%i   %f + i%f\n", i, GSL_REAL(gsl_vector_complex_get(eval, i)), GSL_IMAG(gsl_vector_complex_get(eval, i)));

    //Output (stable)
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, 5); //last one
    gsl_vector_view revec_i = gsl_vector_complex_real(&evec_i.vector);
    halo->stable_direction = gsl_vector_alloc(6);
    gsl_vector_memcpy (halo->stable_direction , &revec_i.vector);

    //Output (unstable)
    evec_i = gsl_matrix_complex_column (evec, 0); //first one
    revec_i = gsl_vector_complex_real(&evec_i.vector);
    halo->unstable_direction = gsl_vector_alloc(6);
    gsl_vector_memcpy (halo->unstable_direction , &revec_i.vector);

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_eigen_nonsymmv_free (w);
    gsl_matrix_free(STM0);

    return GSL_SUCCESS;

}


void halo_orbit_plot(Orbit *halo, gnuplot_ctrl *h1, int type, int points, int scale, double duration,custom_ode_structure *ode_s_6)
{
    //Checking initialization of h1
    if(h1 == NULL) h1 = gnuplot_init();

    //Init
    double t, t1, ti;
    double x1[points+1];
    double x2[points+1];
    double y[6];
    int i;
    for(i=0; i<6; i++) y[i]  = halo->y0[i];
    double scale_factor;

    //Scale
    switch(scale)
    {
    case ADIM:
        scale_factor = 1;
        break;
    case DIM:
        scale_factor = halo->cr3bp.L;
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
    t1 = duration;//halo.T;
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
    gnuplot_plot_xy(h1, x1, x2, points+1, "");

    //Reset
    reset_ode_structure(ode_s_6);

}

void jacobi_halo_orbit_plot(Orbit *halo, gnuplot_ctrl *h1, int points, double duration, custom_ode_structure *ode_s_6)
{
    //Checking initialization of h1
    if(h1 == NULL) h1 = gnuplot_init();

    //Init
    double t, t1, ti;
    double x1[points+1];
    double x2[points+1];
    double y[6];
    int i;
    for(i=0; i<6; i++) y[i]  = halo->y0[i];


    //First value
    x1[0] = 0;
    x2[0] = jacobi(y, halo->cr3bp.mu);


    //Integration over a period
    t = 0;
    t1 = duration;//halo.T;
    for(i =1; i<=points; i++)
    {
        ti = i * t1 / points;
        gsl_odeiv2_driver_apply (ode_s_6->d, &t, ti, y);

        //First value
        x1[i] = t;
        x2[i] = (jacobi(y, halo->cr3bp.mu)- x2[0])/x2[0]*100;
    }
    x2[0] = 0;

    //Plotting
    gnuplot_set_xlabel (h1,"time [adim]");
    gnuplot_set_ylabel (h1,"jacobi constant error [%]");


    gnuplot_setstyle(h1, "lines");
    gnuplot_plot_xy(h1, x1, x2, points+1, "");

    //Reset
    reset_ode_structure(ode_s_6);

}



void third_order_halo_orbit_plot(CR3BP cr3bp, Orbit halo, gnuplot_ctrl *h1, int type, int points, int scale)
{
    //Init
    double t1, ti;
    int i;
    double x1[points+1];
    double x2[points+1];
    double y[6];
    double scale_factor;

    //Scale
    switch(scale)
    {
    case ADIM:
        scale_factor = 1;
        break;
    case DIM:
        scale_factor = cr3bp.L;
        break;
    default: //default: ADIM
        scale_factor = 1;
        break;
    }

    t1 = halo.T;
    for(i =0; i<=points; i++)
    {
        ti = i * t1 / points;
        //Third order approximation of the state
        third_order_halo_orbit(halo, ti, y);

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
    gnuplot_plot_xy(h1, x1, x2, points+1, "THIRD ORDER");

}

