#ifndef SINGLE_ORBIT_H_INCLUDED
#define SINGLE_ORBIT_H_INCLUDED

#include "env.h"
#include "ode.h"
#include "pmcoc.h"

extern "C"{
 #include "nrutil.h"
}

/**
 *  \struct SingleOrbit
 *  \brief  Defines a given orbit with proper arrays to store results.
 **/
typedef struct SingleOrbit SingleOrbit;
struct SingleOrbit
{
    //-----------
    //Parent
    //-----------
    //-----------
    QBCP_L   *qbcp_l;              //QBCP around a given Li point (parent)

    //-----------
    //Parameterization (common to all orbits)
    //-----------
    vector<Oftsc>*  W;             //z(t) = W(s(t), t)
    vector<Oftsc>*  Wh;            //zh(t) = Wh(s(t), t)
    Ofsc* ofs;                     //Auxiliary Ofs object
    double  n;                     //Pulsation of the QBCP
    int order;                     //order of the pm
    int ofs_order;                 //order of the Fourier coefficients
    bool isGS;                     //was the pm obtained through graph style?

    //-----------
    //COC (common to all orbits)
    //-----------
    matrix<Ofsc>* PC;               //COC matrix
    matrix<Ofsc>* CQ;               //COC matrix
    vector<Ofsc>* V;                //COC vector

    //Characteristics
    //-----------
    double   *z0;                    //Initial position in NC coordinates dim = 6
    double   *si;                    //Initial RCM configuration dim = REDUCED_NV
    double   *s0d;                   //Initial position in CCM8 coordinates (real+imag part) dim = 2*REDUCED_NV
    cdouble  *s0;                    //Initial position in CCM8 coordinates (real+imag part) dim = 4
    double   *xf;                    //Final position dim = 6
    double    tf;                    //final time after computation
    double    t0;                    //initial time
    double    tproj;                 //default time between each projection
    double    tprojmin;              //minimum time between each projection
    
    //-----------
    //ODE integration
    //-----------
    OdeStruct *driver;              //NC ode struct
};


/**
    \brief Initialize one SingleOrbit structure
 **/
void init_orbit(SingleOrbit &orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Ofsc>*  PC,
                matrix<Ofsc>*  CQ,
                vector<Ofsc>*  V,
                Ofsc* orbit_ofs,
                int order,
                int ofs_order,
                int isGS,
                double t0,
                double tf,
                double tproj,
                OdeStruct *driver,
                QBCP_L *qbcp_l);

/**
    \brief Free one orbit
 **/
void free_orbit(SingleOrbit *orbit);

/**
    \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], double t0);

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn);

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid of maximum size N
 **/
int trajectory_integration_variable_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn);;

/**
 * \brief Projection on the center manifold
 **/
void NCprojCCMtoCUS(double *yv, double tv, double sti[5], SingleOrbit &orbit, double epsilon, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc);
#endif // SINGLE_ORBIT_H_INCLUDED
