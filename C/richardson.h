#ifndef RICHARDSON_COEFFICIENTS_H_INCLUDED
#define RICHARDSON_COEFFICIENTS_H_INCLUDED

#include "define_env.h"
#include "libration_points.h"
#include "orbits.h"

//Halo orbit family
#define NORTHERN 1
#define SOUTHERN 3

/* Structures */
typedef struct RichardsonCoeffs RichardsonCoeffs;

struct RichardsonCoeffs{

double a21, a22, a23, a24, b21, b22, d21, a31, a32, b31, b32, d31, d32, s1, s2, l1, l2;
double Delta;
double lambda, omega_p, omega_v, kappa;
double c2, c3, c4;
double d1, d2;

};

/* Routines */

//Basic
void richardson_coefficients(CR3BP cr3bp, int pointNumber, RichardsonCoeffs *RC);
double cn(double mu, double liM2dist, int ss, int n);

//Period
double third_order_period(Orbit orbit);

//Third order approximation (various methods)
void third_order_halo_orbit(Orbit halo, double t, double yv_li[]);  //results in yv_li
void apply_halo_orbit_third_order(Orbit *orbit);                //results in orbit.y0

void print_richardson_coefficients(CR3BP cr3bp, int pointNumber);

#endif // RICHARDSON_COEFFICIENTS_H_INCLUDED
