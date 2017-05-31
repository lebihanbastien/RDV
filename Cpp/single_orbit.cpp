#include "single_orbit.h"


/**
 *   \brief Initialize one SingleOrbit structure
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
                QBCP_L *qbcp_l)
{
    //-----------
    //Parameterization (common to all orbits)
    //-----------
    orbit.W          =  W;            //z(t) = W(s(t), t)
    orbit.Wh         =  Wh;           //zh(t) = Wh(s(t), t)

    //-----------
    //COC (common to all orbits)
    //-----------
    orbit.PC  = PC;            //COC matrix
    orbit.CQ  = CQ;            //inv COC matrix
    orbit.V   = V;             //COC vector

    //-----------
    //Orders
    //-----------
    orbit.ofs        =  orbit_ofs;    //Auxiliary Ofs object
    orbit.order      =  order;        //Order of the expansions
    orbit.ofs_order  =  ofs_order;    //Order of the Fourier coefficients
    orbit.isGS       =  isGS;         //Was the pm obtained through graph style?

    //-----------
    //Characteristics
    //-----------
    orbit.z0  = dvector(0, 5);                 //Initial position in NC coordinates dim = 6
    orbit.si  = dvector(0, REDUCED_NV-1);      //Initial RCM configuration dim = 4
    orbit.s0d = dvector(0, 2*REDUCED_NV-1);    //Initial position in CCM8 coordinates (real+imag part) dim = 8
    orbit.xf  = dvector(0, 5);                 //Final position NC dim = 6
    orbit.s0  = dcvector(0,REDUCED_NV-1);      //Initial position in CCM4 coordinates (real+imag part) dim = 4
    orbit.t0  = t0;                            //Initial time
    orbit.tf  = tf;                            //Final time after computation
    orbit.tproj  = tproj;                      //default time between each projection
    orbit.tprojmin  = 1e-2;                    //minimum time between each projection

    //-----------
    //ODE integration
    //-----------
    orbit.driver  = driver;              //NC ode struct

    //-----------
    //Parent
    //-----------
    orbit.qbcp_l = qbcp_l;  //QBCP around a given Li point (parent)

    //-----------
    //Pulsation
    //-----------
    orbit.n = qbcp_l->us.n;
}



/**
 *   \brief Free one orbit
 **/
void free_orbit(SingleOrbit *orbit)
{
    //-----------
    //Characteristics
    //-----------
    free_dvector(orbit->z0,  0, 5);
    free_dvector(orbit->si,  0, REDUCED_NV-1);
    free_dvector(orbit->s0d, 0, 2*REDUCED_NV-1);
    free_dvector(orbit->xf,  0, 5);
    free_dcvector(orbit->s0, 0, REDUCED_NV-1);
    //-----------
    //Ode
    //-----------
    free_ode_structure(orbit->driver);
}

/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], double t0)
{
    //------------------------------------------
    // 1. Update si
    //------------------------------------------
    for(int p = 0; p < REDUCED_NV; p++) orbit.si[p] = si[p];

    //------------------------------------------
    // 2. Update s0
    //------------------------------------------
    RCMtoCCM(si, orbit.s0, REDUCED_NV);

    //------------------------------------------
    // 2. Update s0d
    //------------------------------------------
    RCMtoCCM8(si, orbit.s0d);

    //------------------------------------------
    // 4. Update z0
    //------------------------------------------
    //z0 = W(si, 0.0)
    RCMtoNCbyTFC(si,
                 t0,
                 orbit.n,
                 orbit.order,
                 orbit.ofs_order,
                *orbit.Wh,
                *orbit.ofs,
                *orbit.PC,
                *orbit.V,
                 orbit.z0,
                 orbit.isGS);
}



/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(SingleOrbit &orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *ePm,
                   int *nreset,
                   int isResetOn)
{
    int status;
    double yvp[6], yvi[6];
    cdouble scp[REDUCED_NV];

    //----------------------
    //Projection tools
    //----------------------
    double omega1 = cimag(orbit.Wh->at(0).getCoef(1,0)->getCoef(0));
    double omega3 = cimag(orbit.Wh->at(2).getCoef(1,1)->getCoef(0));
       
    //----------------------
    //Evolve one step of z(t)
    //----------------------
    status = gsl_odeiv2_evolve_apply (orbit.driver->e, orbit.driver->c, orbit.driver->s, &orbit.driver->sys, t, t1, &orbit.driver->h, yv);
    if (status != 0)
    {
        cout << "error in gslc_dual_step: integration of z(t) has gone wrong. break." << endl;
        return -1;
    }

    //----------------------
    //Projection if necessary
    //----------------------
    if(isResetOn && fabs(*t-t0) > fabs(*nreset*orbit.tproj))
    {
        //----------------------
        // Projection on the center manifold
        //----------------------
        //Get the closest point on the center manifold
        NCprojCCM(yv, *t, SEML.us_em.n, OFS_ORDER, *orbit.CQ, *orbit.V, omega1, omega3, scp, REDUCED_NV);
        //Update the state
        CCMtoNCbyTFC(scp, *t, orbit.qbcp_l->us.n, orbit.order,  orbit.ofs_order,  *orbit.Wh,  *orbit.ofs, *orbit.PC, *orbit.V,  yvp,  orbit.isGS);
        //For comparison
        for(int i = 0; i <6; i++) yvi[i] = yv[i];
        // Copy of yvp in current state
        for(int i=0; i<6; i++) yv[i]  = yvp[i];

        //-----------------
        // Get the current projection error
        //-----------------
        //Get the current error
        *ePm = fabs(yvi[0] - yv[0]);
        for(int i = 1; i <6 ; i++)
        {
            if(fabs(yvi[i] - yv[i]) > *ePm) *ePm = fabs(yvi[i] - yv[i]);
        }

        if(*ePm > 1e-6)
        {
            cout << "Warning: Reset n° " << *nreset << ". ePm = " << *ePm << endl;
            return -2;
        }

        //-----------------
        //Reset ode structure for next step
        //-----------------
        reset_ode_structure(orbit.driver);

        //-----------------
        //One additional reset
        //-----------------
        *nreset = *nreset +1;
    }


    return 0;
}



/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(SingleOrbit &orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *ePm,
                     int *nreset,
                     int isResetOn)
{
    reset_ode_structure(orbit.driver);
    int status;
    do
    {
        status = gslc_proj_step(orbit, yv, t, t0, t1, ePm, nreset, isResetOn);

    }while(status == 0 && fabs(*t)<fabs(t1));

    return status;
}



/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------
    //Initialization
    //------------------------------------------
    int status;            //current status
    double yv[6], t;       //current state and time

    //Projection tools
    double ePm;
    int nreset, nt;

    //Plot
    double ti;

    //Change sign of step if necessary
    if((tf < t0 && orbit.driver->h>0) || (tf > t0 && orbit.driver->h<0)) orbit.driver->h *= -1;

    //------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------
    do
    {
        //Init the state & time
        for(int i = 0; i < 6; i++) yv[i] = orbit.z0[i];
        t = t0;
        nt = 0;
        nreset = 1;

        //Loop
        do
        {
            ti = t0 + (double) nt *(tf-t0)/N;
            status = gslc_proj_evolve(orbit, yv, &t, t0, ti, &ePm, &nreset, isResetOn);

            for(int k = 0; k < 6; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = ti;

            //Advance one step
            nt++;
        }
        while((nt<=N) && (status == 0) && (orbit.tproj > orbit.tprojmin));

        //If a new reset is necessary
        if (status == -2 && isResetOn)
        {
            cout << "Warning in trajectory_integration_grid: the interval of projection has to be reduced: ";
            cout << setprecision(3) << "orbit.tproj : " << orbit.tproj << " -> ";
            orbit.tproj *= 0.5;
            cout << orbit.tproj << setprecision(15) << endl;
        }

    }
    while(status!= 0 && orbit.tproj > orbit.tprojmin);

    if(orbit.tproj < orbit.tprojmin)
    {
        cout << "Error in trajectory_integration_grid: the interval of projection is too small." << endl;
        return -1;
    }

    return 0;
}


/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid of maximum size N. Return the last position that is filled on the grid.
 **/
int trajectory_integration_variable_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------
    //Initialization
    //------------------------------------------
    int status;            //current status
    double yv[6], t;       //current state and time

    //Projection tools
    double ePm;
    int nreset, nt;

    //Change sign of step if necessary
    if((tf < t0 && orbit.driver->h>0) || (tf > t0 && orbit.driver->h<0)) orbit.driver->h *= -1;

    //------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------
    do
    {
        //Reset ode structure.
        reset_ode_structure(orbit.driver);

        //Init the state & time
        for(int i = 0; i < 6; i++) yv[i] = orbit.z0[i];
        t = t0;
        nt = 1;
        nreset = 1;

        //First step
        for(int k = 0; k < 6; k++) yNCE[k][0] = orbit.z0[k];
        tNCE[0] = t0;

        //Loop
        do
        {
            status = gslc_proj_step(orbit, yv, &t, t0, tf, &ePm, &nreset, isResetOn);
            for(int k = 0; k < 6; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = t;
            //Advance one step
            nt++;
        }
        while((t < tf) && (nt <= N) && (status == 0) && (orbit.tproj > orbit.tprojmin));

        if (status == -2 && isResetOn)
        {
            cout << "Warning in trajectory_integration_variable_grid: the interval of projection has to be reduced: ";
            cout << setprecision(3) << "orbit.tproj : " << orbit.tproj << " -> ";
            orbit.tproj *= 0.5;
            cout << orbit.tproj << setprecision(15) << endl;
        }

        if(nt == N)
        {
            cout << "Warning in trajectory_integration_variable_grid: the final time was not reached because the maximum number of points is reached." << endl;

        }

    }
    while(status!= 0 && orbit.tproj > orbit.tprojmin);

    if(orbit.tproj < orbit.tprojmin)
    {
        cout << "Error in trajectory_integration_grid: the interval of projection is too small." << endl;
        return -1;
    }

    return nt-1;
}


void NCprojCCMtoCUS(double *yv, double tv, double sti[5], SingleOrbit &orbit, double epsilon, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc)
{
    //Projection tools
    double omega1 = cimag(orbit.Wh->at(0).getCoef(1,0)->getCoef(0));
    double omega3 = cimag(orbit.Wh->at(2).getCoef(1,1)->getCoef(0));
    cdouble scp[REDUCED_NV];

    //Get the closest point on the center manifold, in scp[4]
    NCprojCCM(yv, tv, SEML.us.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, REDUCED_NV);
    //Get the correspondance in RCM coordinates
    CCMtoRCM(scp, sti, REDUCED_NV);
    //Add a given quantity on the hyperbolic direction
    sti[4] = epsilon;
}