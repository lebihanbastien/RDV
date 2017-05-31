/*=========================================================================
 * ode78.c 
 * 
 * C source file that contains the basic routines that appear in the mex files
 *  - ode78_cr3bp.c
 *  - ode78_cr3bp_event.c
 *
 * author:  BLB
 * year:    2015
 * version: 1.0
 *=======================================================================*/

//-------------------------------------------------------------------------
// Headers
//-------------------------------------------------------------------------
//Custom
#include "ode78.h"

//Precisions
double eps_AbsTol=1e-13;
double eps_RelTol=1e-13;
double eps_Diff=1e-13;
double eps_Root=1e-13;
/*-------------------------------------------------------------------------
 * Integration of the CRTBP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. 
 * 
 * The inputs are:
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *
 *  Remarks:
 *   - The number of variables is either 6 (full state) 
 *          or 42 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      cr3bp_derivatives_6 or cr3bp_derivatives_42
 * ----------------------------------------------------------------------*/
void ode78_cr3bp(double *t,                  //current time
                 double *y,                  //current state
                 double const *y0,           //initial condition
                 double const t0,            //initial time
                 double tf,                  //final time
                 int nvar,                   //number of state variables
                 double *mu)                 //cr3bp mass ratio
{
    //---------------------------------------------------------------------
    // Initialize the integration structures
    //---------------------------------------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method

    //General ode structure
    custom_ode_structure ode_s;
    switch(nvar)
    {
        case 6:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 6, 1e-6, cr3bp_derivatives_6, NULL, mu);
            break;
        }
        case 42:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 42, 1e-6, cr3bp_derivatives_42, NULL, mu);
            break;
        }
        default:   //if nvar !=6 && nvar != 42, return without integration
        {
            mexPrintf("ode78_cr3bp. Error: wrong number nvar of state variables. nvar must 6 or 42. return.");
            return;
        }
    
    }
    
    
    //---------------------------------------------------------------------
    // Initalization of the state
    //---------------------------------------------------------------------
    *t = t0;
    for(int i = 0; i < nvar; i++) y[i] = y0[i];    
    
    //---------------------------------------------------------------------
    //Starting in the right direction
    //---------------------------------------------------------------------
    ode_s.d->h = (tf>t0) ? fabs(ode_s.d->h) : -fabs(ode_s.d->h);
    ode_s.h = (tf>t0) ? fabs(ode_s.h) : -fabs(ode_s.h);
    
    //---------------------------------------------------------------------
    // Integration until t = tf
    //---------------------------------------------------------------------
    //Apply ode
    gsl_odeiv2_driver_apply (ode_s.d, t, tf, y); 
}


/*-------------------------------------------------------------------------
 * Integration of the CRTBP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. The state is stored along the
 * trajectory on a given grid.
 * 
 * The inputs are:
 *  - int nGrid,         nGrid+1 is the number of point on the grid
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *  - double *tv, the time on the grid [0, ..., nGrid]
 *  - double *yv, the state on the grid [0, ..., nGrid]
 *
 *  Remarks:
 *   - The number of variables is either 6 (full state) 
 *          or 42 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      cr3bp_derivatives_6 or cr3bp_derivatives_42
 * ----------------------------------------------------------------------*/
void ode78_cr3bp_vec(double *t,                  //current time
                     double *y,                  //current state
                     double *tv,                 //time on a given grid of size [nGrid+1]
                     double **yv,                //state on a given grid of size [nGrid+1]
                     int    nGrid,               //size of the grid
                     double const *y0,           //initial condition
                     double const t0,            //initial time
                     double tf,                  //final time
                     int nvar,                   //number of state variables
                     double *mu)                 //cr3bp mass ratio
{
    //---------------------------------------------------------------------
    // Initialize the integration structures
    //---------------------------------------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method

    //General ode structure
    custom_ode_structure ode_s;
    switch(nvar)
    {
        case 6:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 6, 1e-6, cr3bp_derivatives_6, NULL, mu);
            break;
        }
        case 42:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 42, 1e-6, cr3bp_derivatives_42, NULL, mu);
            break;
        }
        default:   //if nvar !=6 && nvar != 42, return without integration
        {
            mexPrintf("ode78_cr3bp. Error: wrong number nvar of state variables. nvar must 6 or 42. return.");
            return;
        }
    
    }
    
    
    //---------------------------------------------------------------------
    // Initialization of the state
    //---------------------------------------------------------------------
    *t = t0;
    for(int i = 0; i < nvar; i++) y[i] = y0[i];   
    
    //---------------------------------------------------------------------
    //Starting in the right direction
    //---------------------------------------------------------------------
    ode_s.d->h = (tf>t0) ? fabs(ode_s.d->h) : -fabs(ode_s.d->h);
    ode_s.h = (tf>t0) ? fabs(ode_s.h) : -fabs(ode_s.h);
    
    //---------------------------------------------------------------------
    // Integration until t = t1
    // TODO: case when t0 > t1 !!!!
    //---------------------------------------------------------------------
    //Loop on the time grid
    double ti;
    for(int ki = 0; ki <= nGrid; ki++)
    {
         //Current time on the grid
         ti = t0 + (double) ki*(tf-t0)/nGrid;
         //Apply ode
         gsl_odeiv2_driver_apply (ode_s.d, t, ti, y); 
         //Store the current state
         tv[ki] = *t;
         for(int i = 0; i < nvar; i++) yv[i][ki] = y[i];
    }
}




/*-------------------------------------------------------------------------
 * Integration of the BCP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. 
 * 
 * The inputs are:
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *
 *  Remarks:
 *   - The number of variables is either 6 (full state) 
 *          or 42 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      bcp_derivatives_6 or bcp_derivatives_42
 * ----------------------------------------------------------------------*/
void ode78_bcp(double *t,                    //current time
                 double *y,                  //current state
                 double const *y0,           //initial condition
                 double const t0,            //initial time
                 double tf,                  //final time
                 int nvar,                   //number of state variables
                 double *param)              //bcp parameters
{
    //---------------------------------------------------------------------
    // Initialize the integration structures
    //---------------------------------------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method

    //General ode structure
    custom_ode_structure ode_s;
    switch(nvar)
    {
        case 6:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 6, 1e-6, bcp_derivatives_6, NULL, param);
            break;
        }
        case 42:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 42, 1e-6, bcp_derivatives_42, NULL, param);
            break;
        }
        default:   //if nvar !=6 && nvar != 42, return without integration
        {
            mexPrintf("ode78_bcp. Error: wrong number nvar of state variables. nvar must 6 or 42. return.");
            return;
        }
    
    }
    
    
    //---------------------------------------------------------------------
    // Initalization of the state
    //---------------------------------------------------------------------
    *t = t0;
    for(int i = 0; i < nvar; i++) y[i] = y0[i];    
    
    //---------------------------------------------------------------------
    //Starting in the right direction
    //---------------------------------------------------------------------
    ode_s.d->h = (tf>t0) ? fabs(ode_s.d->h) : -fabs(ode_s.d->h);
    ode_s.h = (tf>t0) ? fabs(ode_s.h) : -fabs(ode_s.h);
    
    //---------------------------------------------------------------------
    // Integration until t = tf
    //---------------------------------------------------------------------
    //Apply ode
    gsl_odeiv2_driver_apply (ode_s.d, t, tf, y); 
}


/*-------------------------------------------------------------------------
 * Integration of the BCP vector field of mass ratio mu, from the initial
 * state y0(t0) to the final time tf. The state is stored along the
 * trajectory on a given grid.
 * 
 * The inputs are:
 *  - int nGrid,         nGrid+1 is the number of point on the grid
 *  - double const *y0,  the initial conditions
 *  - double const t0,   initial time
 *  - double tf,         the final time
 *  - int nvar,          the number of state variables
 *  - double *mu,        the cr3bp mass ratio
 *
 * The outputs are:
 *  - double *t,  the current time (should be equal to tf at the end).
 *  - double *y,  the current state, equal to y(tf) at the end of the process
 *  - double *tv, the time on the grid [0, ..., nGrid]
 *  - double *yv, the state on the grid [0, ..., nGrid]
 *
 *  Remarks:
 *   - The number of variables is either 6 (full state) 
 *          or 42 (full state + State Transition Matrix)
 *
 *   - The vector field is computed through either the routine 
 *      bcp_derivatives_6 or bcp_derivatives_42
 * ----------------------------------------------------------------------*/
void ode78_bcp_vec(double *t,                  //current time
                     double *y,                  //current state
                     double *tv,                 //time on a given grid of size [nGrid+1]
                     double **yv,                //state on a given grid of size [nGrid+1]
                     int    nGrid,               //size of the grid
                     double const *y0,           //initial condition
                     double const t0,            //initial time
                     double tf,                  //final time
                     int nvar,                   //number of state variables
                     double *param)              //bcp parameters
{
    //---------------------------------------------------------------------
    // Initialize the integration structures
    //---------------------------------------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method

    //General ode structure
    custom_ode_structure ode_s;
    switch(nvar)
    {
        case 6:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 6, 1e-6, bcp_derivatives_6, NULL, param);
            break;
        }
        case 42:
        {
            init_ode_structure(&ode_s,   T, T_root, eps_AbsTol, eps_RelTol, eps_Root, eps_Diff, 42, 1e-6, bcp_derivatives_42, NULL, param);
            break;
        }
        default:   //if nvar !=6 && nvar != 42, return without integration
        {
            mexPrintf("ode78_cr3bp. Error: wrong number nvar of state variables. nvar must 6 or 42. return.");
            return;
        }
    
    }
    
    
    //---------------------------------------------------------------------
    // Initialization of the state
    //---------------------------------------------------------------------
    *t = t0;
    for(int i = 0; i < nvar; i++) y[i] = y0[i];   
    
    //---------------------------------------------------------------------
    //Starting in the right direction
    //---------------------------------------------------------------------
    ode_s.d->h = (tf>t0) ? fabs(ode_s.d->h) : -fabs(ode_s.d->h);
    ode_s.h = (tf>t0) ? fabs(ode_s.h) : -fabs(ode_s.h);
    
    //---------------------------------------------------------------------
    // Integration until t = t1
    // TODO: case when t0 > t1 !!!!
    //---------------------------------------------------------------------
    //Loop on the time grid
    double ti;
    for(int ki = 0; ki <= nGrid; ki++)
    {
         //Current time on the grid
         ti = t0 + (double) ki*(tf-t0)/nGrid;
         //Apply ode
         gsl_odeiv2_driver_apply (ode_s.d, t, ti, y); 
         //Store the current state
         tv[ki] = *t;
         for(int i = 0; i < nvar; i++) yv[i][ki] = y[i];
    }
}