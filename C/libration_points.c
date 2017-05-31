#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Custom
#include "libration_points.h"
#include "define_env.h"

/* -------------------------------------------------------------------  *

*        Position of the libration points in the CR3BP                  *

* --------------------------------------------------------------------  */
/**
* \brief Using the Newton-Raphson method, find the root of a function known to lie close to x1 The root
*  /rtnewt will be refined until its accuracy is known within ± xacc.
*  funcd is a user-supplied routine that returns both the function value and the first derivative of the
*  function at the point x.
**/
double rtnewt(void (*funcd)(double, int, double, double *, double *), double x1, double xacc, double mu, int number)

{
    void nrerror(char error_text[]);
    int j;
    double df,dx,f,rtn;

    rtn=x1;   //initial guess

    for (j=1; j<=JMAX; j++)
    {
        (*funcd)(mu, number, rtn,&f,&df);
        dx=f/df;
        rtn -= dx;

        if (fabs(dx) < xacc) return rtn;  //Convergence
    }
    printf("WARNING: Maximum number of iterations exceeded in rtnewt");
    return 0.0;   //Never get here.
}


/**
* \brief Provides the function value and its first derivative for the newton-raphson method.
* f corresponds to the equation satisfied by the Li-m2 distance for the L1/L2 cases
* and by 1-(Li-m1 distance) for the L3 case
**/
void polynomialLi(double mu, int number, double y, double *f, double *df)
{
    switch(number)
    {

    case 1:
        *f =  pow(y,5.0)   - (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) +  2*mu*y - mu;
        *df = 5*pow(y,4.0) - 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    +  2*mu;
        break;

    case 2:
        *f =  pow(y,5.0)   + (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) -  2*mu*y - mu;
        *df = 5*pow(y,4.0) + 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    -  2*mu;
        break;

    case 3:
        //*f =  pow(y,5.0) + (2.0+mu)*pow(y,4.0) + (1+2*mu)*pow(y,3.0) + (1+mu)*pow(y,2.0) +  2*(1-mu)*y + 1-mu;
        *f= pow(y,5.0) + (7+mu)*pow(y,4.0) + (19+6*mu)*pow(y,3.0) -(24+13*mu)*pow(y,2.0) +  (12+14*mu)*y -7*mu;
        *df= 5*pow(y,4.0) + 4*(7+mu)*pow(y,3.0) + 3*(19+6*mu)*pow(y,2.0) -2*(24+13*mu)*pow(y,1.0) +  (12+14*mu);
        //*df = 5*pow(y,4.0) + 4*(2.0+mu)*pow(y,3.0) + 3*(1+2*mu)*pow(y,2.0) + 2*(1+mu)*y +  2*(1-mu);
        break;

    }

}
