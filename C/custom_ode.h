#ifndef CUSTOM_ODE_H_INCLUDED
#define CUSTOM_ODE_H_INCLUDED

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>

typedef struct custom_ode_structure custom_ode_structure;
struct custom_ode_structure
{
    //Stepper
    const gsl_odeiv2_step_type *T;
    gsl_odeiv2_step *s;
    //Control
    gsl_odeiv2_control *c;
    //Evolution
    gsl_odeiv2_evolve *e;
    //System
    gsl_odeiv2_system sys;
    //Driver
    gsl_odeiv2_driver * d;
    //Solver for root finding
    gsl_root_fsolver *s_root;

    //Precisions
    double eps_int_rel;  //for integration (relative precision)
    double eps_int_abs;  //for integration (absolute precision)
    double eps_root;     //for root finding
    double eps_diff;     //for differential correction

    //Initial step
    double h;

    //Dimension
    int dim;
};


void update_ode_structure(custom_ode_structure *ode_s,
                        const gsl_odeiv2_step_type *T,
                        gsl_odeiv2_step *s,
                        gsl_odeiv2_control *c,
                        gsl_odeiv2_evolve *e,
                        gsl_odeiv2_system sys,
                        gsl_odeiv2_driver * d,
                        gsl_root_fsolver *s_root,
                        double eps_root,
                        double eps_diff);

int init_ode_structure(custom_ode_structure *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        double eps_int_abs,
                        double eps_int_rel,
                        double eps_root,
                        double eps_diff,
                        size_t dim,
                        double h,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params),
                        void *params);

void reset_ode_structure(custom_ode_structure *ode_s);
void free_ode_structure(custom_ode_structure *ode_s);


#endif // CUSTOM_ODE_H_INCLUDED
