#include "pmcoc.h"


/**
 * \file pmcoc.cpp
 * \brief Change of coordinates for the parameterization methods.
 * \author BLB.
 * \date 2016
 * \version 1.0
 *
 *  Two types of change of coordinates are implemented:
 *              1. Changes between different manifold coordinates (Real Center to Complex Center...).
 *              2. Evaluations of the parameterization (Real Center to Normalized-Centered and projection the other way around).
 *
 */


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Change of coordinates
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief from CCM to RCM coordinates
 **/
void CCMtoRCM(const cdouble s1[], double si[], int nv)
{
    si[0] = creal(1.0/sqrt(2)*(s1[0]   + s1[2]*I));
    si[2] = creal(1.0/sqrt(2)*(s1[0]*I + s1[2]));
    si[1] = creal(1.0/sqrt(2)*(s1[1]   + s1[3]*I));
    si[3] = creal(1.0/sqrt(2)*(s1[1]*I + s1[3]));
    if(nv == 5)
    {
        si[4] = creal(s1[4]);
    }
    else if(nv == 6)
    {
        si[4] = creal(s1[4]);
        si[5] = creal(s1[5]);
    }
}

/**
 *  \brief from RCM to CCM coordinates
 **/
void RCMtoCCM(const double si[], cdouble s1[], int nv)
{
    //From real to complex TFC
    s1[0] = 1.0/sqrt(2)*(si[0] - si[2]*I);
    s1[2] = 1.0/sqrt(2)*(si[2] - si[0]*I);
    s1[1] = 1.0/sqrt(2)*(si[1] - si[3]*I);
    s1[3] = 1.0/sqrt(2)*(si[3] - si[1]*I);
    if(nv == 5)
    {
        s1[4] = si[4]+I*0.0;
    }
    else if(nv == 6)
    {
        s1[4] = si[4]+I*0.0;
        s1[5] = si[5]+I*0.0;
    }
}

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void RCMtoCCM8(const double si[], double s0d[])
{
    //From real to complex TFC
    cdouble s1[REDUCED_NV];
    RCMtoCCM(si, s1, REDUCED_NV);

    //Store real and imag part separately
    CCMtoCCM8(s1, s0d);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void CCM8toRCM(const double s0d[], double si[])
{
    //CCM8 to CCM
    cdouble s1[REDUCED_NV];
    CCM8toCCM(s0d, s1);
    //CCM to RCM
    CCMtoRCM(s1, si, REDUCED_NV);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void CCM8toCCM(const double s0d[], cdouble s1[])
{
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        s1[p]  =   s0d[p2++];
        s1[p] += I*s0d[p2++];
    }
}

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void CCMtoCCM8(const cdouble s1[], double s0d[])
{
    //Store real and imag part separately
    int p2 = 0;
    for(int p = 0; p < REDUCED_NV; p++)
    {
        s0d[p2++] = creal(s1[p]);
        s0d[p2++] = cimag(s1[p]);
    }
}

/**
 *  \brief from TFC to TF coordinates
 **/
void TFCtoTF(const cdouble s1[6], double si[6])
{
    //First center
    si[0] = creal(1.0/sqrt(2)*(s1[0]   + s1[3]*I));
    si[3] = creal(1.0/sqrt(2)*(s1[0]*I + s1[3]));
    //Second center
    si[2] = creal(1.0/sqrt(2)*(s1[2]   + s1[5]*I));
    si[5] = creal(1.0/sqrt(2)*(s1[2]*I + s1[5]));
    //Hyperbolic dir
    si[1] = creal(s1[1]);
    si[4] = creal(s1[4]);
}


//---------------------------------------------------------------------------------------------------------------------------------------
// COC: NC <--> SYS
//---------------------------------------------------------------------------------------------------------------------------------------
//From Normalized-Centered coordinates to system coordinates
void NCtoSYS(double t, const double yNC[], double yEM[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}


//-----------------------------------------------------------------------------
// COC: NC <--> EM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Normalized-Centered coordinates to Earth-Moon coordinates
 **/
void NCtoEM(double t, const double yNC[], double yEM[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];

}

/**
 *  \brief COC: from Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void EMtoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Sun-Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void SEMtoNC(double t, const double ySEM[], double yNC[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -ySEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -ySEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +ySEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -ySEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -ySEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +ySEM[5]/gamma;
}

/**
 *  \brief COC: from Normalized-Centered coordinates to Sun-Earth-Moon coordinates
 **/
void NCtoSEM(double t, const double yNC[], double ySEM[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    ySEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    ySEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    ySEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    ySEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    ySEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    ySEM[5] = +gamma*yNC[5];

}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluation of the pm
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void RCMtoNC(const double st0[],
             const double t,
             const double n,
             const int order,
             const int ofs_order,
             vector<Oftsc> &W,
             Ofsc &ofs,
             double z1[],
             bool isGS)
{
    //------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------
    cdouble s0[REDUCED_NV];
    cdouble z0[6];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, REDUCED_NV);
    //------------------------------------------
    // 2. Update z0
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            if(p == 2 || p == 5)
            {
                //order 1 is sufficient for p = 2,5
                W[p].evaluate(s0, ofs, 1, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
                z1[p] = creal(z0[p]);
            }
            else
            {
                //For p = 0,1,3,4 normal computation
                W[p].evaluate(s0, ofs, order, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
                z1[p] = creal(z0[p]);
            }
        }
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            W[p].evaluate(s0, ofs, order, ofs_order);
            z0[p] = ofs.evaluate(n*t, ofs_order);
            z1[p] = creal(z0[p]);
        }
    }

}

/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector:  z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoNCbyTFC(const double st0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  vector<Oftsc> &Wh,
                  Ofsc &ofs,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // RCM to TFC
    //------------------------------------------
    RCMtoTFC(st0, order, ofs_order, Wh, zIn, isGS);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}

/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param s0 an array of 4 complex double which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoNCbyTFC(cdouble s0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  vector<Oftsc> &Wh,
                  Ofsc &ofs,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // CCM to TFC
    //------------------------------------------
    CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}


/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoTFC(const double st0[],
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              vector<Ofsc> &zIn,
              bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[REDUCED_NV];

    //------------------------------------------
    // 1. RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, REDUCED_NV);

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);
}


/**
 *  \brief Apply the change of variables in zIN/zOut. (OFS version).
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Ofsc> &zIn,
              vector<Ofsc> &zOut)
{
    //zeroing the target
    for(unsigned int i = 0; i < zOut.size(); i++) zOut[i].zero();
    //zOut = PC*zIn
    smvprod_ofs(PC, zIn, zOut);
    //zOut+=V(theta)
    for(int i = 0; i < (int) zOut.size(); i++) zOut[i] += V[i];
}


/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(s0), t)
 *   \param s0 an array of 4 complex which gives the configuration to input in complex CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoTFC(cdouble s0[],
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              vector<Ofsc> &zIn,
              bool isGS)
{
    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    if(isGS)
    {


        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        cdouble temp;
        // zIn[0]
        //---------------
        temp = Wh[0].getCoef(1,0)->getCoef(0);
        zIn[0].setCoef(temp*s0[0],0);

        // zIn[1]
        //---------------
        Wh[1].evaluate(s0, zIn[1], order, ofs_order);
        // zIn[2]
        //---------------
        temp = Wh[2].getCoef(1,1)->getCoef(0);
        zIn[2].setCoef(temp*s0[1],0);
        // zIn[3]
        //---------------
        temp = Wh[3].getCoef(1,2)->getCoef(0);
        zIn[3].setCoef(temp*s0[2],0);
        // zIn[4]
        //---------------
        Wh[4].evaluate(s0, zIn[4], order, ofs_order);
        // zIn[5]
        //---------------
        temp = Wh[5].getCoef(1,3)->getCoef(0);
        zIn[5].setCoef(temp*s0[3],0);
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            Wh[p].evaluate(s0, zIn[p], order, ofs_order);
        }
    }
}

/**
 *  \brief Projection of the current NC state on the central manifold, via CCM coordinates
 **/
void NCprojCCM(const double z[], const double t, const double n, const int ofs_order, matrix<Ofsc> &CQ, vector<Ofsc> &V, double omega1, double omega3, cdouble sc[], int nv)
{
    //-------------------
    //Wh: TFC coordinates
    //-------------------
    cdouble zh[6];

    //-------------------
    //z - V
    //-------------------
    cdouble zd[6];
    for(int p = 0; p < 6; p++) zd[p] = z[p] - V[p].evaluate(n*t, ofs_order);

    //-------------------
    //Update Wh = CQ*(z - V)
    //-------------------
    for(int k = 0; k <6; k++)
    {
        zh[k] = 0.0+0.0*I;
        for(int p = 0; p <6; p++)
        {
            zh[k] += zd[p]* CQ.getCoef(k, p).evaluate(n*t, ofs_order);
        }
    }

    //-------------------
    //Projection on the center manifold
    //-------------------
    TFCprojCCM(zh, omega1, omega3, sc, nv);
}

/**
 *  \brief Projection of the current TFC state on the central manifold, via CCM coordinates
 **/
void TFCprojCCM(const cdouble zh[], double omega1, double omega3, cdouble sc[], int nv)
 {
    //-------------------
    //Projection on the center manifold
    //-------------------
    //Wh1 = i*w1*s1 => s1 = -i/w1*Wh1
    sc[0] = -1.0*I/omega1*zh[0];

    //Wh3 = i*w3*s2 => s2 = -i/w3*Wh3
    sc[1] = -1.0*I/omega3*zh[2];

    //Wh4 = -i*w1*s3 => s3 = +i/w1*Wh4
    sc[2] = +1.0*I/omega1*zh[3];

    //Wh6 = -i*w3*s4 => s4 = +i/w3*Wh6
    sc[3] = +1.0*I/omega3*zh[5];

    if(nv > 4) sc[4] = 0.0;
    if(nv > 5) sc[5] = 0.0;
 }


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluation at time t
//
//---------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoef(double *alpha, double t, double omega, int order, double *params, int number)
{
    double *header;
    int l;
    double cR[order];
    double sR[order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    for(l = 0, header = params; l < number ; l++, header+=(order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13) alpha[l] = evaluateOdd(t, omega , order, header, sR);    //Odd funtions (alpha_2,5,8,10,12,14)
        else  alpha[l] = evaluateEven(t, omega, order, header, cR);                                                  //Even functions
    }
}

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoefDerivatives(double *alpha, double t, double omega, int order, double *params, int number)
{
    double *header;
    int l;
    double cR[order];
    double sR[order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    for(l = 0, header = params; l < number ; l++, header+=(order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13)
        {
            alpha[l] = evaluateOddDerivative(t, omega , order, header, cR); //Odd funtions (alpha_2,5,8,10,12,14)
        }
        else  alpha[l] = evaluateEvenDerivative(t, omega, order, header, sR);   //Even functions
    }
}


//-----------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double evaluateEven(double t, double omega, int order, double *coef, double *cR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += coef[i]*cR[i-1];//even type
    result += coef[0];
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double evaluateEvenDerivative(double t, double omega,  int order, double *coef, double *sR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += -omega*i*coef[i]*sR[i-1];//even type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double evaluateOdd(double t, double omega,  int order, double *coef, double *sR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += coef[i]*sR[i-1]; //odd type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double evaluateOddDerivative(double t, double omega,  int order, double *coef, double *cR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += omega*i*coef[i]*cR[i-1];//odd type
    return result;
}

//-----------------------------------------------------------------------------
// Evaluating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evz(Ofsc& zt, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*zt.evaluate(n*t);
}

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzdot(Ofsc& zt, Ofsc& ztdot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*(ztdot.evaluate(n*t) + I*ni*zt.evaluate(n*t));
}

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzddot(Ofsc& zt, Ofsc& ztdot, Ofsc& ztddot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*( 2*I*ni*ztdot.evaluate(n*t) - ni*ni*zt.evaluate(n*t) + ztddot.evaluate(n*t));
}

