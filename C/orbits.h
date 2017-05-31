#ifndef ORBITS_H_INCLUDED
#define ORBITS_H_INCLUDED

#include <gsl/gsl_matrix.h>

//Custom
#include "define_env.h"

#define HALO_EMPTY                -1
#define HALO_APPROXIMATED          0
#define HALO_REAL                  1

/*---------------------------------
            Structures
------------------------------------*/
// Orbit structure
typedef struct Orbit Orbit;
struct Orbit
{

//----------------------------------------------------------------------------
// Basic parameters
//----------------------------------------------------------------------------
    int m;               //m parameter (1 for Class I, 3 for Class II)
    int dm;              //corresponding dm to m. BEWARE: Northern and Southern classes are opposite to Class I and Class II formulation (see Richardson 1980)
    CR3BP cr3bp;         //CR3BP system
    LibrationPoint Li;   //Libration point

//status = HALO_EMPTY if its starting point was not initialized or if there is a discrepancy between its parameters (e.g. Az_R and y0)
//status = HALO_APPROXIMATED if a third-order approximation was used to generate the orbit from its Az value
//status = HALO_REAL if a differential correction scheme was applied
    int status;

//Approximation for third order computation (R stands for RICHARDSON)
    double Az_R;           //Vertical extension (adim)
    double Azdim_R;        //Vertical extension (dim)

//----------------------------------------------------------------------------
// Require computation
//---------------------------------------------------------------------------
    double y0[42];       //starting point @ y = 0;
    double yhalf[42];     //half-period point
    double T;            //orbital period (adim)

//TRUE VALUES
    double Az_T;           //Vertical extension (adim)
    double Azdim_T;        //Vertical extension (dim)

//Dynamic parameters
    gsl_matrix* monodromy;           //monodromy matrix
    gsl_vector* stable_direction;    //stable eigenvector;
    gsl_vector* unstable_direction;  //unstable eigenvector;


};

/*---------------------------------
            Routines
------------------------------------*/
/**
* \brief Initialization of halo orbit
* WARNING: the fields;
                      y0,
                      yhalf,
                      T,
                      monodromy,
                      stable_direction,
                      unstable_direction,
                      Az_T,
                      Azdim_T
* are not initialized since they require additionnal computation
**/
void init_halo_orbit(Orbit *orbit, CR3BP cr3bp , LibrationPoint Li,  int m, double Azdim);

/**
* \brief Change the Azdim_R and Az_R value
**/
void set_Az_R(Orbit *orbit, double Azdim);


#endif // ORBITS_H_INCLUDED
