//Custom
#include "orbits.h"

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
void init_halo_orbit(Orbit *orbit, CR3BP cr3bp , LibrationPoint Li,  int m, double Azdim_R)
{
    //CR3BP
    orbit->cr3bp = cr3bp;

    //Class
    orbit->m = m;

    switch(Li.number)
    {
    case 1:
        orbit->dm = 2-m;
        break;
    case 2:
        orbit->dm = m-2;  //BEWARE: Northern and Southern classes are opposite to Class I and Class II formulation (see Richardson 1980)
        break;
    case 3:
        orbit->dm = 2-m;
        break;
    default:
        orbit->dm = 0;  //NO dm if Li.number!=1,2,3
        break;
    }

    //Az
    orbit->Azdim_R = Azdim_R;
    orbit->Az_R = Azdim_R/cr3bp.L;

    //Li
    orbit->Li = Li;

    //status
    orbit->status = HALO_EMPTY;
}

/**
* \brief Change the Azdim_R and Az_R value
**/
void set_Az_R(Orbit *orbit, double Azdim)
{
    //Azdim_R
    orbit->Azdim_R = Azdim;

    //Az_R
    orbit->Az_R = Azdim/orbit->cr3bp.L;

    //Updating status
    orbit->status = HALO_EMPTY;
}
