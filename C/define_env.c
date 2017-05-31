#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

//Custom
#include "define_env.h"
#include "libration_points.h"
#include "energy.h"

/* -------------------------------------------------------------------  *

*  Define the working environment (the Sun, planets and the moon)       *

*  All values taken from JPL and Goddard Space Flight Center websites   *

* --------------------------------------------------------------------  */

/**
* \brief Initialize the Circular Restricted 3-Body Problem
* \param cr3bp pointer on the CR3BP
* \param n1 name of the first primary
* \param n2 name of the second primary
**/
void init_CR3BP(CR3BP *cr3bp, int n1, int n2)
{
    //Body initialization
    init_body(&(*cr3bp).m1, n1);
    init_body(&(*cr3bp).m2, n2);

    cr3bp->mu =  (*cr3bp).m2.GM/( (*cr3bp).m1.GM + (*cr3bp).m2.GM );        // µ = m2/(m1 + m2)
    cr3bp->L = (*cr3bp).m2.a;                                               // Distance parameter = semi major axis of m2
    cr3bp->T = (*cr3bp).m2.T;                                               //Time parameter = sidereal period of m2
    cr3bp->R1 = (*cr3bp).m1.Req;
    cr3bp->R2 = (*cr3bp).m2.Req;
    cr3bp->rh = pow((*cr3bp).mu/3,1/3.0);                                   //Hill's radius adim formula
    cr3bp->gprecision = 1e-12;                                               //arbitrary

    //Li initialization
    init_libp(&cr3bp->l1, *cr3bp, 1);
    init_libp(&cr3bp->l2, *cr3bp, 2);
    init_libp(&cr3bp->l3, *cr3bp, 3);
    init_libp(&cr3bp->l4, *cr3bp, 4);
    init_libp(&cr3bp->l5, *cr3bp, 5);

    //Name
    strcpy(cr3bp->name, cr3bp->m1.name);
    strcat(cr3bp->name, "-");
    strcat(cr3bp->name, cr3bp->m2.name);

    //Distance to manifold approximation
    if(n1 == EARTH && n2 == MOON) cr3bp->d_man = 50/cr3bp->L; //50km for the Earth-Moon system
    else cr3bp->d_man = 100/cr3bp->L; //100km otherwise (to be changed for specific systems!)
}


/**
* \brief Initializes a libration point
* **/
void init_libp(LibrationPoint *libp, CR3BP cr3bp, int number)
{

    double gamma_i;

    //Number
    libp->number = number;

    switch(number)
    {
    case 1:
        //Gamma
        gamma_i = cr3bp.rh - 1/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);                      //initial guess
        gamma_i = rtnewt(polynomialLi, gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);   //newton-raphson method
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu - gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;
        break;

    case 2:
        //Gamma
        gamma_i = cr3bp.rh + 1/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);
        gamma_i = rtnewt(polynomialLi, gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu + gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;
        break;

    case 3:
        //Gamma
        gamma_i = 7/12.0*cr3bp.mu + pow(237,2.0)/pow(12,4.0)*pow(cr3bp.mu,3.0);
        gamma_i = rtnewt(polynomialLi, gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);
        libp->gamma_i = 1-gamma_i;  //BEWARE: for L3, gamma3 = L3-M1 distance != L3-M2


        //Position
        libp->position[0] = -1 - cr3bp.mu + gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;
        break;

    case 4:
        //Gamma
        libp->gamma_i = 1;

        //Position
        libp->position[0] = -cr3bp.mu + 0.5;
        libp->position[1] = sqrt(3)/2.0;
        libp->position[2] = 0;
        break;

    case 5:
        //Gamma
        libp->gamma_i = 1;

        //Position
        libp->position[0] = -cr3bp.mu + 0.5;
        libp->position[1] = -sqrt(3)/2.0;
        libp->position[2] = 0;
        break;
    }

    //Energy & Jacobi constant
    libp->Ei = energy(libp->position, cr3bp.mu);
    libp->Ci = -2*libp->Ei;
}


/**
* \brief Initialize one celestial body
* \param body pointer on the current body
* \param name the name of the body in integer format (consistent with HORIZON numerotation)
**/
void init_body(Body *body, int name)
{

    double days = 86400; //days to seconds

    switch(name)
    {

    case MERCURY:

        //Physical parameters
        body->Req = 2439.7;        //[km]
        body->Rm = 2439.7;         //[km]
        body->M = 0.330104e24;     //[kg]
        body->GM = 22032;          //[km^3.s^-2]

        //Orbital parameters
        body->a = 57.91e6;         //[kg]
        body->T = 87.9691*days;     //[s]

        strcpy(body->name, "Mercury");
        break;

    case VENUS:

        //Physical parameters
        body->Req = 6051.8;        //[km]
        body->Rm = 6501.8;         //[km]
        body->M = 4.86732e24;      //[kg]
        body->GM = 324858.63;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 108.21e6;        //[km]
        body->T = 224.701*days;     //[s]

        strcpy(body->name, "Venus");
        break;


    case EARTH:

        //Physical parameters
        body->Req = 6378.14;        //[km]
        body->Rm = 6371.00;         //[km]
        body->M = 5.97219e24;       //[kg]
        body->GM = 398600.440;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;     //[s]

        strcpy(body->name, "Earth");
        break;

    case MOON:

        //Physical parameters
        body->Req = 1737.5;       //[km]
        body->Rm = 1737.5;        //[km]
        body->M =  0.07342e24;    //[kg]
        body->GM = 4902.801;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 384400;           //[km]
        body->T = 27.321582*days;    //[s]

        strcpy(body->name, "Moon");
        break;

    case MARS:

        //Physical parameters
        body->Req = 3396.19;       //[km]
        body->Rm = 3389.50;        //[km]
        body->M = 0.641693e24;     //[kg]
        body->GM = 42828.3;        //[km^3.s^-2]

        //Orbital parameters
        body->a = 227.92e6;       //[kg]
        body->T = 686.98*days;     //[s]

        strcpy(body->name, "Mars");
        break;


    case JUPITER:

        //Physical parameters
        body->Req =  71492;      //[km]
        body->Rm = 69911;        //[km]
        body->M = 1898.13e24;     //[kg]
        body->GM = 126686511;       //[km^3.s^-2]

        //Orbital parameters
        body->a = 778.57e6;       //[kg]
        body->T = 4332.589*days;     //[s]

        strcpy(body->name, "Jupiter");
        break;

    case SATURN:

        //Physical parameters
        body->Req =  60268;    //[km]
        body->Rm = 58232;       //[km]
        body->M = 568.319e24;     //[kg]
        body->GM = 37931207.8;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 1433.53e6;       //[kg]
        body->T = 10759.22*days;     //[s]

        strcpy(body->name, "Saturn");
        break;

    case URANUS:

        //Physical parameters
        body->Req = 25559;      //[km]
        body->Rm = 25362;        //[km]
        body->M =  86.8103e24;    //[kg]
        body->GM =  5793966;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  2872.46e6;      //[kg]
        body->T =  30685.4*days;   //[s]

        strcpy(body->name, "Uranus");
        break;

    case NEPTUNE:

        //Physical parameters
        body->Req = 24764;      //[km]
        body->Rm = 24622;        //[km]
        body->M = 102.410e24;     //[kg]
        body->GM =  6835107;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  4495.06e6;      //[kg]
        body->T =  60189*days;    //[s]

        strcpy(body->name, "Neptune");
        break;

    case PLUTO:

        //Physical parameters
        body->Req =  1195;     //[km]
        body->Rm =  1195;       //[km]
        body->M = .01309e24;     //[kg]
        body->GM =  872.4;      //[km^3.s^-2]

        //Orbital parameters
        body->a =  5906.38e6;      //[kg]
        body->T =  90465*days;    //[s]

        strcpy(body->name, "Pluto");
        break;


    case SUN:

        //Physical parameters
        body->Req = 696342;                //[km]
        body->Rm =  696342;                //[km]
        body->M  = 1988500e24;             //[kg]
        body->GM = 1.3271244004193938e11;  //[km^3.s^-2]

        //Orbital parameters
        body->a = 0;    //[kg]
        body->T = 0;    //[s]

        strcpy(body->name, "Sun");
        break;

    case EARTH_AND_MOON:
        //Equivalent mass of the Earth+Moon system based at the center of mass
        //additionnal physical properties are those of the Earth for consistency)

        //Physical parameters
        body->Req = 6378.14;        //[km]
        body->Rm = 6371.00;         //[km]
        body->M = 5.97219e24+0.07342e24;       //[kg]
        body->GM = 398600.440+4902.801;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;     //[s]

        strcpy(body->name, "Earth+Moon");
        break;


    }

}
