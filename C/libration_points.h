/* -------------------------------------------------------------------  *

*        Position of the libration points in the CR3BP                  *

* --------------------------------------------------------------------  */
#ifndef LIBRATION_POINTS_H_INCLUDED
#define LIBRATION_POINTS_H_INCLUDED

//#include "define_env.h"

/*---------------------------------
            Structures
------------------------------------*/
typedef struct LibrationPoint LibrationPoint;
struct LibrationPoint
{
    int number;           //number
    double position[3];   //position [adim] in the corresponding CR3BP-frame
    double gamma_i;       //distance to closest primary (only for l1,l2,l3);
    double Ei;            //energy
    double Ci;            //jacobi constant
};
/*---------------------------------
            Routines
------------------------------------*/
/**
* \brief Using the Newton-Raphson method, find the root of a function known to lie close to x1 The root
*  /rtnewt will be refined until its accuracy is known within ± xacc.
*  funcd is a user-supplied routine that returns both the function value and the first derivative of the
*  function at the point x.
**/
double rtnewt(void (*funcd)(double, int, double, double *, double *), double x1, double xacc, double mu, int number);

/**
* \brief Provides the function value and its first derivative for the newton-raphson method.
* f corresponds to the equation satisfied by the Li-m2 distance for the L1/L2 cases
* and by 1-(Li-m1 distance) for the L3 case
**/
void polynomialLi(double mu, int number, double y, double *f, double *df);

/**
* \brief compute the energy and jacobi constant of a given li point
**/

#endif // LIBRATION_POINTS_H_INCLUDED
