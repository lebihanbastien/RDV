#ifndef CUSTOM_ODEZERO_H_INCLUDED
#define CUSTOM_ODEZERO_H_INCLUDED

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl_interp.h>
#include <gsl/gsl_roots.h>

//Custom
#include "define_env.h"
#include "custom_ode.h"
/*---------------------------------
            Structures
------------------------------------*/
//Set the parameters to trigger an event
struct value_params
{
    int max_events;
    int direction;
    int dim;
    double value;
    double *center;
    int type;
};

//Structure of the output of a value function for event procedure
struct value_output
{
    double val;
    //int max_events;
    //int direction; 
};

//Structure of a value function for event procedure
struct value_function
{
    struct value_params *val_par;
    double (*value)(double, double [], void *);
};

//Structure dedicated to ode parameters
struct ode_params
{
    int nvar;
    double t0;
    double* y0;
    gsl_odeiv2_driver *d;
    struct value_function fvalue;
};


//Structure dedicated to ode parameters
struct cr3bp_params
{
    double t0;
    double* y0;
    gsl_odeiv2_driver *d;
};

/*---------------------------------
            Routines
------------------------------------*/
/*int custom_odezero(double y[],
                   double *tcross,
                   double t1,
                   const gsl_odeiv2_step_type *T,
                   double eps_int,
                   double eps_root,
                   gsl_root_fsolver *s_root,
                   gsl_odeiv2_step *s,
                   gsl_odeiv2_control *c,
                   gsl_odeiv2_evolve *e,
                   gsl_odeiv2_driver * d,
                   gsl_odeiv2_system sys,
                   CR3BP cr3bp);*/

int custom_odezero(double y[], double *tcross, double t1, custom_ode_structure *ode_s);

/*int custom_odezero_2(double y[],
                     double** ye,
                     double *tcross,
                     double t1,
                     const gsl_odeiv2_step_type *T,
                     double eps_int,
                     double eps_root,
                     gsl_root_fsolver *s_root,
                     gsl_odeiv2_step *s,
                     gsl_odeiv2_control *c,
                     gsl_odeiv2_evolve *e,
                     gsl_odeiv2_driver * d,
                     gsl_odeiv2_system sys,
                     CR3BP cr3bp,
                     struct value_function fvalue);*/

int custom_odezero_2(double y[],
                     double** ye,
                     double *tcross,
                     double t1,
                     custom_ode_structure *ode_s,
                     struct value_function fvalue);

double linear_intersection(double t, double yv[], void *params);
double angle_intersection(double t, double yv[], void *params);
double null_flight_path_angle(double t, double yv[], void *params);
double cr3bp_event (double t, void *params);
double cr3bp_y (double t, void *params);

#endif // CUSTOM_ODEZERO_H_INCLUDED
