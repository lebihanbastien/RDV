#ifndef HALO_ORBIT_COMPUTATION_H_INCLUDED
#define HALO_ORBIT_COMPUTATION_H_INCLUDED

//GSL
#include <gsl/gsl_odeiv2.h>
//Gnuplot
#include "gnuplot_i.h"

//Plot type
#define XY 1
#define XZ 2
#define YZ 3

//scale
#define ADIM 0
#define DIM  1

//void halo_orbit_computation(CR3BP cr3bp, Orbit *halo, const gsl_odeiv2_step_type *T, int diff_corr_type);
/*int halo_orbit_computation(CR3BP cr3bp,
                           Orbit *halo,
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
int halo_orbit_computation(Orbit *halo,
                           custom_ode_structure *ode_s,
                           diff_corr_parameters diff_corr);
//void halo_orbit_refinement(CR3BP cr3bp, Orbit *halo, const gsl_odeiv2_step_type *T, double ystart[], int diff_corr_type);
/*int halo_orbit_refinement(CR3BP cr3bp,
                           Orbit *halo,
                           const gsl_odeiv2_step_type *T,
                           double ystart[],
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
int halo_orbit_refinement(Orbit *halo,
                          double ystart[],
                          custom_ode_structure *ode_s,
                          diff_corr_parameters diff_corr);
void halo_orbit_plot(Orbit *halo, gnuplot_ctrl *h1, int type, int points, int scale, double duration,custom_ode_structure *ode_s_6);
void jacobi_halo_orbit_plot(Orbit *halo, gnuplot_ctrl *h1, int points, double duration, custom_ode_structure *ode_s_6);
void third_order_halo_orbit_plot(CR3BP cr3bp, Orbit halo, gnuplot_ctrl *h1, int type, int points, int scale);

#endif // HALO_ORBIT_COMPUTATION_H_INCLUDED
