#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include <complex.h>

//------------------------------------------------------------------------------------
// Environment
//------------------------------------------------------------------------------------
/// Numerotation of the Solar System planets and objects, consistent with JPL's HORIZON numerotation
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
/// Custom indix for the Sun+Earth system
#define EARTH_AND_MOON 700
/// Precision on the position of the librations points L1/L2/L3 in define_env.h
#define LIBRATION_POINT_PRECISION 1e-16


//------------------------------------------------------------------------------------
//   Global constants
//------------------------------------------------------------------------------------
extern int OFTS_ORDER;
extern int OFS_ORDER;
extern int OTS_ORDER;
extern int MODEL_TYPE;
extern int REDUCED_NV;

//------------------------------------------------------------------------------------
//   typedef
//------------------------------------------------------------------------------------
typedef complex double cdouble;


//------------------------------------------------------------------------------------
//   ORDER AND NUMBER OF VARIABLES FOR OFS AND OFTS OBJECTS
//------------------------------------------------------------------------------------
//Model
#define M_RTBP  0 // RTBP model indix
#define M_QBCP  1 // QBCP model indix
#define M_BCP   2 // BCP model indix
#define M_ERTBP 3 // ERTBP model indix
//Frameworks
#define F_EM 0
#define F_SEM 1
#define F_SE 2

//Order of the potential of the primaries
#define POTENTIAL_ORDER 60

//Number of variables
#define NV 6
/// Number of variables in the OFS object (a priori always 1)
#define OFS_NV 1


#endif // PARAMETERS_H_INCLUDED
