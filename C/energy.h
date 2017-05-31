#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

/**
* \brief Compute the potential energy for the given state and CR3BP (mu)
**/
double energy(double y[], double mu);

/**
* \brief Compute the Jacobi constant for the given state and CR3BP (mu)
**/
double jacobi(double y[], double mu);

#endif // ENERGY_H_INCLUDED
