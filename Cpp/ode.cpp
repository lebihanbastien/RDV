//Custom
#include "ode.h"

#include <iostream>
using namespace std;

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

void update_ode_structure(OdeStruct *ode_s,
                         const gsl_odeiv2_step_type *T,
                         gsl_odeiv2_step *s,
                         gsl_odeiv2_control *c,
                         gsl_odeiv2_evolve *e,
                         gsl_odeiv2_system sys,
                         gsl_odeiv2_driver * d,
                         gsl_root_fsolver *s_root,
                         double eps_root,
                         double eps_diff)
{
    ode_s->dim = d->sys->dimension;
    ode_s->T = T;
    ode_s->s = s;
    ode_s->c = c;
    ode_s->e = e;
    ode_s->sys = sys;
    ode_s->d = d;
    ode_s->s_root = s_root;

    ode_s->eps_root = eps_root;
    ode_s->eps_diff = eps_diff;
    ode_s->h = d->h;

    //Integration tolerances are set to -1 because they are hidden in d and e (be careful during initialization!)
    ode_s->eps_int_abs = -1;
    ode_s->eps_int_rel = -1;
}

void init_ode_structure(OdeStruct *ode_s,
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
                        void *params)
{
    //Precisions and initial step
    ode_s->eps_root = eps_root;
    ode_s->eps_diff = eps_diff;
    ode_s->h = h;
    ode_s->eps_int_abs = eps_int_abs;
    ode_s->eps_int_rel = eps_int_rel;
    ode_s->dim = dim;

    //Stepper
    ode_s->T = T;
    ode_s->s = gsl_odeiv2_step_alloc(T, dim);
    if (ode_s->s == NULL)
    {
        cout << "failed to allocate stepper object in init_ode_structure" << endl;
    }

    //Control
    ode_s->c = gsl_odeiv2_control_y_new(eps_int_abs, eps_int_rel);
    if (ode_s->c == NULL)
    {
        cout << "failed to allocate control object in init_ode_structure"  << endl;
    }

    //Evolution
    ode_s->e = gsl_odeiv2_evolve_alloc(dim);
    if (ode_s->e == NULL)
    {
        cout << "failed to allocate evolution object in init_ode_structure"  << endl;
    }


    //System
    ode_s->sys.function = func;
    ode_s->sys.jacobian = jacobian;
    ode_s->sys.dimension = dim;
    ode_s->sys.params = params;

    ode_s->d = gsl_odeiv2_driver_alloc_y_new (&ode_s->sys, ode_s->T, ode_s->h, ode_s->eps_int_abs, ode_s->eps_int_rel);
    if (ode_s->d == NULL)
    {
        cout << "failed to allocate driver object in init_ode_structure"  << endl;
    }

    ode_s->s_root = gsl_root_fsolver_alloc (T_root);
    if (ode_s->s_root == NULL)
    {
        cout << "failed to allocate root_fsolver object in init_ode_structure"  << endl;
    }

}

void reset_ode_structure(OdeStruct *ode_s)
{
    gsl_odeiv2_step_reset(ode_s->s);
    gsl_odeiv2_evolve_reset(ode_s->e);
    gsl_odeiv2_driver_reset(ode_s->d);
}


void free_ode_structure(OdeStruct *ode_s)
{
    gsl_odeiv2_step_free(ode_s->s);
    gsl_odeiv2_evolve_free(ode_s->e);
    gsl_odeiv2_driver_free(ode_s->d);
}
