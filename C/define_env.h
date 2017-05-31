/* -------------------------------------------------------------------  *

*  Define the working environment (the Sun, planets and the moon)       *

*  All values taken from JPL and Goddard Space Flight Center websites   *

* --------------------------------------------------------------------  */
#ifndef DEFINE_ENV_H_INCLUDED
#define DEFINE_ENV_H_INCLUDED

#include "libration_points.h"

// Consistent with JPL's HORIZON numerotation
#define SUN 10
#define MERCURY 199
#define VENUS 299
#define EARTH 399
#define MARS 499
#define JUPITER 599
#define SATURN 699
#define URANUS 799
#define NEPTUNE 899
#define PLUTO 999
#define MOON 301
#define EARTH_AND_MOON 700 //Custom indix for the Sun+Earth system

#define JMAX 50                          //Set to maximum number of iterations.
#define LIBRATION_POINT_PRECISION 1e-12  //Precision on the position of the librations points L1/L2/L3

/*---------------------------------
            Structures
------------------------------------*/
// Body
typedef struct Body Body;
struct Body
{
    // Physical parameters
    double M;    //mass [kg]
    double GM;   //Gravitationnal parameter [km^3/s^2]
    double Req;  //Equatorial Radius [km]
    double Rm;   //Mean radius [km]

    // Orbital parameters
    double T;    //Sidereal orbit period [s]
    double a;    //Semi-major axis [km]

    //Name
    char name[50];
};


// Environment of the CR3BP
typedef struct CR3BP CR3BP;
struct CR3BP
{

    Body m1;   //First primary
    Body m2;   //Second primary

    double mu; // Gravitational constant [-]
    double L;  // Distance parameter [km]
    double T;  // Time parameter [s]
    double R1; // Radius of the first primary  [km]
    double R2; // Radius of the second primary [km]
    double rh; // Hill's radius [-]

    //Libration points (adim)
    LibrationPoint l1;
    LibrationPoint l2;
    LibrationPoint l3;
    LibrationPoint l4;
    LibrationPoint l5;

    //precision on L1/L2/L3 position
    double gprecision;

    //Name
    char name[50];

    //Distance to manifold approximation
    double d_man;
};


/*---------------------------------
            Routines
------------------------------------*/
/**
* \brief Initialize the Circular Restricted 3-Body Problem
* \param cr3bp pointer on the CR3BP
* \param n1 name of the first primary
* \param n2 name of the second primary
**/
void init_body(Body *body, int name);

/**
* \brief Initialize one celestial body
* \param body pointer on the current body
* \param name the name of the body in integer format (consistent with HORIZON numerotation)
**/
void init_CR3BP(CR3BP *cr3bp, int n1, int n2);

/**
* \brief Initializes a libration point
* **/
void init_libp(LibrationPoint *libp, CR3BP cr3bp, int number);


#endif // DEFINE_ENV_H_INCLUDED
