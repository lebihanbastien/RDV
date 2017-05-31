#ifndef CR3BP_DERIVATIVES_H_INCLUDED
#define CR3BP_DERIVATIVES_H_INCLUDED

#include <gsl/gsl_matrix.h>


int cr3bp_derivatives_42 (double t, const double y[], double f[], void *params);
int cr3bp_derivatives_6 (double t, const double y[], double f[], void *params);

int bcp_derivatives_6 (double t, const double y[], double f[], void *params);
int bcp_derivatives_42 (double t, const double y[], double f[], void *params);

#endif // CR3BP_DERIVATIVES_H_INCLUDED
