#ifndef DIFFERENTIAL_CORRECTION_H_INCLUDED
#define DIFFERENTIAL_CORRECTION_H_INCLUDED


#define X0_FIXED 1
#define Z0_FIXED 2

#include "define_env.h"
#include "cr3bp_derivatives.h"
#include "custom_odezero.h"
#include "user_parameters.h"

/*int differential_correction(double ystart[], double *tcross, double t1, const gsl_odeiv2_step_type *T,
                                             double eps_int, double eps_root, double eps_diff, CR3BP cr3bp, int type);*/

/*int differential_correction(CR3BP cr3bp,
                            double ystart[],
                            double yhalf[],
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
                            double eps_diff,
                            int diff_corr_type);*/

int differential_correction(CR3BP cr3bp,
                            double ystart[],
                            double yhalf[],
                            double *tcross,
                            double t1,
                            custom_ode_structure *ode_s,
                            diff_corr_parameters diff_corr);

#endif // DIFFERENTIAL_CORRECTION_H_INCLUDED
