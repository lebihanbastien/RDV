#include <math.h>


/**
* \brief Compute the potential energy for the given state and CR3BP (mu)
**/
double energy(double y[], double mu)
{
    double r1 = sqrt( pow(y[0]+mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );
    double r2 = sqrt( pow(y[0]- 1 + mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );

    return - ( 1/2.0*(pow(y[0],2) + pow(y[1],2)) + (1-mu)/r1 + mu/r2 + 1/2.0*mu*(1-mu) );
}


/**
* \brief Compute the Jacobi constant for the given state and CR3BP (mu)
**/
double jacobi(double y[], double mu)
{
    return -2*energy(y, mu) - (pow(y[3],2.0)+pow(y[4],2.0)+pow(y[5],2.0));
}


