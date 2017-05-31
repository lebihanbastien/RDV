#ifndef PMCOC_H_INCLUDED
#define PMCOC_H_INCLUDED

#include "Oftsc.h"
#include "env.h"
#include "matrix.h"

/**
 * \file pmcoc.h
 * \brief Change of coordinates for the parameterization methods.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  Two types of change of coordinates are implemented:
 *              1. Changes between different manifold coordinates (Real Center to Complex Center...).
 *              2. Evaluations of the parameterization (Real Center to Normalized-Centered and projection the other way around).
 *
 */

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Change of coordinates in the manifolds
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief from CCM to RCM coordinates
 **/
void CCMtoRCM(const cdouble s1[], double si[], int nv);

/**
 *  \brief from RCM to CCM coordinates
 **/
void RCMtoCCM(const double si[], cdouble s1[], int nv);

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void RCMtoCCM8(const double si[], double s0d[]);

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void CCM8toRCM(const double s0d[], double si[]);

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void CCM8toCCM(const double s0d[], cdouble s1[]);

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void CCMtoCCM8(const cdouble s1[], double s0d[]);

/**
 *  \brief from TFC to TF coordinates
 **/
void TFCtoTF(const cdouble s1[6], double si[6]);


//---------------------------------------------------------------------------------------------------------------------------------------
// Change of coordinates between the systems
//---------------------------------------------------------------------------------------------------------------------------------------
//From Normalized-Centered coordinates to system coordinates
void NCtoSYS(double t, const double yNC[], double yEM[], void *params_void);

//-----------------------------------------------------------------------------
// COC: NC <--> EM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Normalized-Centered coordinates to Earth-Moon coordinates
 **/
void NCtoEM(double t, const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void EMtoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Sun-Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void SEMtoNC(double t, const double ySEM[], double yNC[], QBCP_L *qbp);

/**
 *  \brief COC: from Normalized-Centered coordinates to Sun-Earth-Moon coordinates
 **/
void NCtoSEM(double t, const double yNC[], double ySEM[], QBCP_L *qbp);

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
             bool isGS);

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
                  bool isGS);

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
                  bool isGS);

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
              bool isGS);



/**
 *  \brief Apply the change of variables in zIN/zOut. (OFS version).
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC(matrix<Ofsc> &PC,
              vector<Ofsc> &V,
              vector<Ofsc> &zIn,
              vector<Ofsc> &zOut);


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
              bool isGS);

/**
 *  \brief Projection of the current NC state on the central manifold, via CCM coordinates
 **/
void NCprojCCM(const double z[], const double t, const double n, const int ofs_order, matrix<Ofsc> &CQ, vector<Ofsc> &V, double omega1, double omega3, cdouble sc[], int nv);

/**
 *  \brief Projection of the current TFC state on the central manifold, via CCM coordinates
 **/
void TFCprojCCM(const cdouble zh[], double omega1, double omega3, cdouble sc[], int nv);

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
void evaluateCoef(double *alpha, double t, double omega, int order, double *params, int number);

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoefDerivatives(double *alpha, double t, double omega, int order, double *params, int number);

//-----------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double evaluateEven(double t, double omega, int order, double *coef, double *cR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double evaluateEvenDerivative(double t, double omega,  int order, double *coef, double *sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double evaluateOdd(double t, double omega,  int order, double *coef, double *sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double evaluateOddDerivative(double t, double omega,  int order, double *coef, double *cR);

//-----------------------------------------------------------------------------
// Evaluate the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evz(Ofsc& zr, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzdot(Ofsc& zr, Ofsc& ztdot, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzddot(Ofsc& zr, Ofsc& ztdot, Ofsc& ztddot, double t, double n, double ni, double ai);
#endif // PMCOC_H_INCLUDED
