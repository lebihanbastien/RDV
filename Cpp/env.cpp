/**
 * \file env.cpp
 * \brief Define the working environment (the Sun, planets and the moon). All values taken from JPL and Goddard Space Flight Center websites.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 */

//std
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

//custom
#include "env.h"
#include "parameters.h"
#include "Ofsc.h"

int MODEL_TYPE;


//-----------------------------------------------------------------------------------------------------------------------------
//QBCP
//-----------------------------------------------------------------------------------------------------------------------------
QBCP SEM;              ///global variable that describes the Sun-Earth-Moon system
QBCP_L SEML;           ///global variable that describes the Sun-Earth-Moon system around Li (EM framework)

/**
 *   \brief Initialization of the environnement (Sun, Earth, Moon, Li...).
 *
 *    The global variables SEM and SEML are initialized in order to describe the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem around the Lagrangian point LI,
 *    LI being given as an integer in the file parameters.h. The default coordinates system is the normalized-centered (NC) coordinates of the Earth-Moon system.
 *    Note that, in order for the initialization to be complete - Sun-(Earth+Moon) motion given as Fourier series within SEML -
 *    the routine qbtbp(true) must have been run *once.
 **/
void init_env(int li_EM, int li_SEM, int isNormalized, int model, int fwrk, int pmStyle, int manType_EM, int manType_SEM)
{
    //Init the Sun-Earth-Moon problem
    init_QBCP(&SEM, SUN, EARTH, MOON, model);
    //Init the Sun-Earth-Moon problem focused on one libration point
    init_QBCP_L(&SEML, &SEM, isNormalized, li_EM, li_SEM, false, model, F_EM, pmStyle, manType_EM, manType_SEM);
}


//---------------------------------------
//Center manifold
//---------------------------------------
//vector<Oftsc>  CM;     ///center manifold in NC coordinates
//vector<Oftsc> CMh;     ///center manifold in TFC coordinates

/**
 *   \brief Initialization of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    Of course, these data files must have been previously computed through the use of the routine pm(int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *
 **/
void initCM(QBCP_L &qbcp)
{
    //Memory allocation
    //CM     = vector<Oftsc>(NV);
    //CMh    = vector<Oftsc>(NV);

    //Read from bin files
    updateCM(qbcp);
}

/**
 *   \brief Update of the parameterization center manifold around LI, LI being given as an integer in the file parameters.h.
 *
 *    The parameterization is retrieved from text file given in the folder F_PM, F_PM being a global string array defined in parameters.h.
 *    Of course, these data files must have been previously computed through the use of the routine pm(int, int).
 *    Moreover, as always, the Fourier-Taylor algebra must have been initialized by the routine FTDA::init in the main program.
 *    The global variables initialized by this routine are:
 *    - vector<Oftsc>  CM, center manifold in NC coordinates
 *    - vector<Oftsc> CMh, center manifold in TFC coordinates
 *    - matrix<Oftsc> JCM, jacobian of CM
 *    - vector<Oftsc>  Fh, reduced vector field
 *
 **/
void updateCM(QBCP_L &qbcp)
{
    //Read from bin files
    //readVOFTS_bin(CM,    qbcp.cs.F_PMS+"W/W",         OFS_ORDER);
    //readVOFTS_bin(CMh,   qbcp.cs.F_PMS+"W/Wh",        OFS_ORDER);
}


//---------------------------------------
//COC
//---------------------------------------
matrix<Ofsc>  Mcoc;    ///COC matrix
matrix<Ofsc>  Pcoc;    ///COC matrix (Mcoc = Pcoc*Complex matrix)
matrix<Ofsc>  MIcoc;    ///COC matrix = inv(Mcoc)
matrix<Ofsc>  PIcoc;    ///COC matrix = inv(Pcoc)
vector<Ofsc>  Vcoc;    ///COC vector

/**
 *  Main routine for the update of the COC
 **/
void updateCOC(QBCP_L &qbcp)
{
    //Read from files
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, qbcp);
}


/**
 *  Main routine for the initialization of the COC
 **/
void initCOC(QBCP_L &qbcp)
{
    //Memory allocation
    Pcoc  = matrix<Ofsc>(NV,NV);
    Mcoc  = matrix<Ofsc>(NV,NV);
    PIcoc = matrix<Ofsc>(NV,NV);
    MIcoc = matrix<Ofsc>(NV,NV);
    Vcoc  = vector<Ofsc>(NV);

    //Read from files
    updateCOC(qbcp);
}


//-----------------------------------------------------------------------------------------------------------------------------
//COC
//-----------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Update a complex Fourier series given as the component (k,p) of a matrix matrix, from a given txt file.
 *  \param xFFT: the \c ofs<cdouble> to update
 *  \param matrix: the beginning of the name of the source txt file. e.g. "alpha"
 *  \param k the lign indix of the desired component in the matrix \c matrix.
 *  \param p the column indix of the desired component in the matrix \c matrix.
 *
 *  As an example, the call readCOC(xFFT, "alpha", 2, 1), will update xFFT with the file "alpha21.txt".
 **/
void readCOC(Ofsc& xFFT, string matrix, int k, int p)
{
    //Init
    ifstream readStream;
    string ss1, ss2;
    ss1 = static_cast<ostringstream*>( &(ostringstream() << k) )->str();
    ss2 = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    //Reading an OFS from a text file
    readOFS_txt(xFFT, (matrix+ss1+ss2));
}

/**
 *  \brief Lighter version of the init routine, with only P, PC and V retrieved.
 **/
void initCOC(matrix<Ofsc> &P,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &V,
             QBCP_L& qbcp_l)
{
    //----------------------------
    //Retrieve folder
    //----------------------------
    string F_COC = qbcp_l.cs.F_COC;

    //----------------------------
    //Switch RTBP/QBCP/BCP
    //----------------------------
    if(qbcp_l.model == M_QBCP || qbcp_l.model == M_BCP)
    {
        //Recovering the data: matrix P
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*P.getCA(i,j), F_COC+"P",  i+1, j+1);
        //Recovering the data: matrix Q
        for(int i = 0; i < NV ; i++) for(int j = 0; j < NV ; j++) readCOC(*Q.getCA(i,j), F_COC+"Q",  i+1, j+1);

        //Recovering the data: vector V
        //V = [G1_11 G1_12 0 G1_21 G1_22 0]^T
        //Note:V[2] and V[5] are kept null
        readCOC(V[0], F_COC+"G1_",  1, 1);
        readCOC(V[1], F_COC+"G1_",  1, 2);
        readCOC(V[3], F_COC+"G1_",  2, 1);
        readCOC(V[4], F_COC+"G1_",  2, 2);
    }
    else //EM RTPB
    {
        //note that V is left untouched (set to zero by default)
        //--------------------

        //--------------------
        //Init double variables
        //--------------------
        double eta1, eta2, la1, om1, om2, dl1, do1, s1, s2, c2;
        c2 = qbcp_l.c2;
        eta1 = (c2 - 2.0 - sqrt(9*c2*c2 - 8*c2))/2.0;
        eta2 = (c2 - 2.0 + sqrt(9*c2*c2 - 8*c2))/2.0;
        om1 = sqrt(-eta1);
        la1 = sqrt(+eta2);
        om2 = sqrt(c2);
        dl1 = 2*la1*( (4.0+3*c2)*la1*la1 + 4 + 5*c2 - 6*c2*c2);
        do1 =   om1*( (4.0+3*c2)*om1*om1 - 4 - 5*c2 + 6*c2*c2);
        s1  = sqrt(dl1);
        s2  = sqrt(do1);

        //--------------------
        //Init P
        //--------------------
        P.setCoef(+2*la1/s1,                          0, 1);
        P.setCoef(+(la1*la1  - 2*c2 - 1)/s1,          1, 1);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 1);
        P.setCoef(+(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 1);

        P.setCoef(-(om1*om1 + 2*c2 + 1)/s2,           1, 0);
        P.setCoef(-(om1*om1 - 2*c2 - 1)/s2,           3, 0);

        P.setCoef(-2*la1/s1,                          0, 4);
        P.setCoef(+(la1*la1  - 2*c2 - 1)/s1,          1, 4);
        P.setCoef(+(la1*la1  + 2*c2 + 1)/s1,          3, 4);
        P.setCoef(-(la1*la1*la1 + (1 - 2*c2)*la1)/s1, 4, 4);

        P.setCoef(+2*om1/s2,                          0, 3);
        P.setCoef(-(om1*om1*om1 - (1 - 2*c2)*om1)/s2, 4, 3);

        P.setCoef(+1.0/sqrt(om2),                     2, 2);
        P.setCoef(+sqrt(om2),                         5, 5);

        //--------------------
        //Q = inv(P) (gsl lib)
        //--------------------
        int s;
        gsl_matrix* Pc   = gsl_matrix_alloc (NV, NV);
        gsl_matrix* Qc   = gsl_matrix_alloc (NV, NV);
        gsl_permutation * p6 = gsl_permutation_alloc (NV);

        //Init Pc
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) gsl_matrix_set(Pc, i, j, creal(P.getCA(i,j)->getCoef(0)));
        //Use of GSL library
        gsl_linalg_LU_decomp (Pc, p6, &s);
        gsl_linalg_LU_invert (Pc, p6, Qc);

        //--------------------
        // Init Q
        //--------------------
        for(int i =0; i < NV; i++) for(int j =0; j < NV; j++) Q.setCoef(gsl_matrix_get(Qc, i, j), i, j);
    }

    //----------------------------
    // Building PC
    //----------------------------
    Ofsc BUX(OFS_ORDER);
    //Keymap used to initialize PC and CQ
    vector<int> keyMap(4);
    keyMap[0] = 0;
    keyMap[1] = 1;
    keyMap[2] = 3;
    keyMap[3] = 4;
    keyMap[4] = 2;
    keyMap[5] = 5;


    int ii;
    //Init PC by rows
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,0),   1.0/sqrt(2)+0.0*I, P(ii,3), I*1.0/sqrt(2));
        PC.setCoef(BUX, ii, 0);
        BUX.ofs_fsum(P(ii,0), I*1.0/sqrt(2), P(ii,3),   1.0/sqrt(2)+0.0*I);
        PC.setCoef(BUX, ii, 3);
        PC.setCoef(P(ii,1), ii, 1);
        PC.setCoef(P(ii,4), ii, 4);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(P(ii,2),   1.0/sqrt(2)+0.0*I, P(ii,5), I*1.0/sqrt(2));
        PC.setCoef(BUX, ii, 2);
        BUX.ofs_fsum(P(ii,2), I*1.0/sqrt(2), P(ii,5),   1.0/sqrt(2)+0.0*I);
        PC.setCoef(BUX, ii, 5);
    }

    //Init CQ by columns
    for(int i = 0; i <= 3; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(0,ii),    1.0/sqrt(2)+0.0*I, Q(3,ii), -1.0/sqrt(2)*I);
        CQ.setCoef(BUX, 0, ii);
        BUX.ofs_fsum(Q(0,ii), -1.0/sqrt(2)*I, Q(3,ii),    1.0/sqrt(2)+0.0*I);
        CQ.setCoef(BUX, 3, ii);
        CQ.setCoef(Q(1,ii), 1, ii);
        CQ.setCoef(Q(4,ii), 4, ii);
    }

    for(int i = 4; i <= 5; i++)
    {
        ii = keyMap[i];
        BUX.ofs_fsum(Q(2,ii),    1.0/sqrt(2)+0.0*I, Q(5,ii), -1.0/sqrt(2)*I);
        CQ.setCoef(BUX, 2, ii);
        BUX.ofs_fsum(Q(2,ii), -1.0/sqrt(2)*I, Q(5,ii),    1.0/sqrt(2)+0.0*I);
        CQ.setCoef(BUX, 5, ii);
    }

}



//-----------------------------------------------------------------------------------------------------------------------------
//            Init routines
//-----------------------------------------------------------------------------------------------------------------------------
/**
* \brief Initialize a QBCP_L structure, i.e. a QBCP focused on two libration points: one the EM system and one of the SEM system.
*        The libration point must be L1 or L2 both for the EM and SEM systems.
* \param qbcp_l a pointer on the QBCP_L structure to init.
* \param qbcp a pointer on the QBCP parent structure.
* \param isNormalized: are the equations of motion normalized?
* \param li_EM number of the libration point considered for the EM system
* \param li_SEM number of the libration point considered for the SEM system
* \param isNew an integer: equal to 1 if no solution has been previously computed with the routine qbtbp(), 0 otherwise
* \param model: QBCP, BCP, CRTBP...
* \param fwrk: default coordinate system for this structure: for example:
*                      - if fwrk == F_EM,  the qbcp_l is focused on the li_EM  point of the EM  system.
*                      - if fwrk == F_SEM, the qbcp_l is focused on the li_SEM point of the SEM system.
*        The default focus can be change dynamically during computation.
*
*   Note that the QBCP structure is used only for the initialization of the coordinate systems. More precisely, it contains some parameters
*   specific to each libration point (gamma), via its CR3BP structures.
**/
void init_QBCP_L(QBCP_L *qbcp_l, QBCP *qbcp, int isNormalized, int li_EM, int li_SEM, int isNew, int model, int fwrk, int pmType, int manType_EM, int manType_SEM)
{
    //-------------------------------------------------------------------------------------------------
    //      Common to all models
    //      These settings may be needed to initialize the CSYS and USYS structures
    //      For this reason, they are initialized in priority
    //-------------------------------------------------------------------------------------------------
    //-----------------
    // - Hard-coded value of the order of the Fourier series. Equals zero for RTBP
    // - Note that for the bicircular models, the order of the Fourier series is still equals to OFS_ORDER
    // This may need to be changed when EFFECTIVE_OFS_ORDER will be used, to be able to manipulate Fourier series
    // with an order smaller than the value 30, hard coded in the data files.
    // - Used in the CSYS init functions.
    //-----------------
    switch(model)
    {
    case M_QBCP:
    case M_BCP:
    case M_ERTBP:
        qbcp_l->nf = OFS_ORDER;
        break;
    case M_RTBP:
        qbcp_l->nf = 0;
        break;
    default:
        cout << "init_QBCP_L. Warning: unknown model." << endl;
    }

    //-----------------
    //Normalization or not. True is always preferable - may be deprecated
    //Kept for compatibility with old code.
    //-----------------
    qbcp_l->isNormalized = isNormalized;

    //-----------------
    //Libration point
    //Note: passing the complete libration point structure as argument may be done in the near future.
    //not necessary for now.
    //-----------------
    qbcp_l->li_EM  = li_EM;
    qbcp_l->li_SEM = li_SEM;

    //-----------------
    // Model - used in the CSYS init functions.
    //-----------------
    qbcp_l->model   = model;

    //-----------------
    // Parameterization style - used in the CSYS init functions.
    //-----------------
    qbcp_l->pms   = pmType;

    //-----------------
    // Effective order of the Fourier series in RVF (reduced vector field)
    //-----------------
    switch(pmType)
    {
        case PMS_GRAPH:
        case PMS_MIXED:
            qbcp_l->eff_nf = qbcp_l->nf;    //for graph method, the reduced vector field is non autonomous
        break;
        case PMS_NORMFORM:
            qbcp_l->eff_nf = 0;             //for normal form, the reduced vector field is autonomous!
        break;
        default:
        cout << "init_QBCP_L. Warning: unknown pmType." << endl;
    }

    //-------------------------------------------------------------------------------------------------
    // Unit systems
    //-------------------------------------------------------------------------------------------------
    init_USYS(&qbcp_l->us_em,  USYS_EM,  model);
    init_USYS(&qbcp_l->us_sem, USYS_SEM, model);

    //-------------------------------------------------------------------------------------------------
    // Coordinate systems
    //-------------------------------------------------------------------------------------------------
    switch(model)
    {
    case M_QBCP:
    case M_BCP:
    case M_RTBP:
        qbcp_l->numberOfCoefs = 15;
        break;
    case M_ERTBP:
        qbcp_l->numberOfCoefs = 15;
        break;
    default:
        cout << "init_QBCP_L. Warning: unknown model." << endl;
    }

    //--------------------
    // Earth-Moon
    //--------------------
    init_CSYS(&qbcp_l->cs_em_l1, qbcp_l, qbcp,  F_EM, 1, qbcp_l->numberOfCoefs, isNew, pmType, manType_EM); //L1
    init_CSYS(&qbcp_l->cs_em_l2, qbcp_l, qbcp,  F_EM, 2, qbcp_l->numberOfCoefs, isNew, pmType, manType_EM); //L2
    init_CSYS(&qbcp_l->cs_em_l3, qbcp_l, qbcp,  F_EM, 3, qbcp_l->numberOfCoefs, isNew, pmType, manType_EM); //L3
    //--------------------
    // Sun-Earth
    //--------------------
    init_CSYS(&qbcp_l->cs_sem_l1, qbcp_l, qbcp,  F_SEM, 1, qbcp_l->numberOfCoefs, isNew, pmType, manType_SEM); //L1
    init_CSYS(&qbcp_l->cs_sem_l2, qbcp_l, qbcp,  F_SEM, 2, qbcp_l->numberOfCoefs, isNew, pmType, manType_SEM); //L2
    init_CSYS(&qbcp_l->cs_sem_l3, qbcp_l, qbcp,  F_SEM, 3, qbcp_l->numberOfCoefs, isNew, pmType, manType_SEM); //L3

    //-------------------------------------------------------------------------------------------------
    // DEFAULT SETTINGS
    // The flag fwrk determine the default focus; either on the Earth-Moon or Sun-Earth+Moon framework
    //-------------------------------------------------------------------------------------------------
    //Coord. syst. for the EM point
    //--------------------
    switch(li_EM)
    {
    case 1:
        qbcp_l->cs_em  = qbcp_l->cs_em_l1;
        break;
    case 2:
        qbcp_l->cs_em  = qbcp_l->cs_em_l2;
        break;
    case 3:
        qbcp_l->cs_em  = qbcp_l->cs_em_l3;
        break;
    }

    //--------------------
    //Coord. syst. for the SE point
    //--------------------
    switch(li_SEM)
    {
    case 1:
        qbcp_l->cs_sem  = qbcp_l->cs_sem_l1;
        break;
    case 2:
        qbcp_l->cs_sem  = qbcp_l->cs_sem_l2;
        break;
    case 3:
        qbcp_l->cs_sem  = qbcp_l->cs_sem_l3;
        break;
    }

    //--------------------
    //Default Li, Unit & Coord. system
    //--------------------
    qbcp_l->fwrk = fwrk;
    switch(fwrk)
    {
    case F_EM:
        qbcp_l->us = qbcp_l->us_em;
        qbcp_l->cs = qbcp_l->cs_em;
        qbcp_l->li = qbcp_l->li_EM;
        break;
    case F_SEM:
        qbcp_l->us = qbcp_l->us_sem;
        qbcp_l->cs = qbcp_l->cs_sem;
        qbcp_l->li = qbcp_l->li_SEM;
        break;
    case F_SE:
        qbcp_l->us = qbcp_l->us_sem;
        qbcp_l->cs = qbcp_l->cs_sem;
        qbcp_l->li = qbcp_l->li_SEM;
        break;
    default:
        cout << "init_QBCP_L. Warning: unknown fwrk." << endl;
    }

    //-------------------------------------------------------------------------------------------------
    // 6*6 matrix B used for Floquet transformation
    // Do we need it here?
    //-------------------------------------------------------------------------------------------------
    qbcp_l->B  = (double*) calloc(36, sizeof(double));

    //-------------------------------------------------------------------------------------------------
    // €€ TODO
    //-------------------------------------------------------------------------------------------------
//    //Moon "on" by default
    qbcp_l->epsilon = 1.0;
//    //c2 value --- A REVOIR!!!!!
    qbcp_l->c2 = cn(*qbcp_l, 2);
}

/**
 * \brief Initializes the Circular Restricted 3-Body Problem
 * \param cr3bp pointer on the CR3BP
 * \param n1 name of the first primary
 * \param n2 name of the second primary
 **/
void init_CR3BP(CR3BP *cr3bp, int n1, int n2)
{
    //Body initialization
    init_body(&(*cr3bp).m1, n1);
    init_body(&(*cr3bp).m2, n2);

    cr3bp->mu = (*cr3bp).m2.M/( (*cr3bp).m1.M + (*cr3bp).m2.M );  // µ = m2/(m1 + m2)
    cr3bp->L  = (*cr3bp).m2.a;                                      // Distance parameter = semi major axis of m2
    cr3bp->T  = (*cr3bp).m2.T;                                      //Time parameter = sidereal period of m2
    cr3bp->R1 = (*cr3bp).m1.Req;
    cr3bp->R2 = (*cr3bp).m2.Req;
    cr3bp->rh = pow((*cr3bp).mu/3,1/3.0);                          //Hill's radius adim formula
    cr3bp->gprecision = 1e-12;                                     //arbitrary

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
 * \brief Initializes a unit system in the form of a usys structure, such as Earth-Moon or Sun-(Earth+Moon) unit systems.
 * \param usys pointer on the usys structure to init.
 * \param label the type of unit system
 *
 * NEEDS TO BE MODIFIED TO GET RID OF THE HARD CODED VALUES?
 **/
void init_USYS(USYS *usys, int label, int model)
{
    //EM mass ratio
    usys->mu_EM  = +1.215058162343360e-02;
    //Lunar eccentricity
    usys->lecc   = +0.054900489;


    //Physical params in EM units
    //--------------------------
    double n_EM  = 0.925195985520347;
    double ns_EM = 1.0-n_EM;
    double as_EM = 388.81114;
    double ms_EM = ns_EM*ns_EM*pow(as_EM,3.0)-1;
    double me_EM = 1.0-usys->mu_EM;
    double mm_EM = usys->mu_EM;

    //Rest of the parameters
    switch(model)
    {
    case M_QBCP:
    case  M_BCP:
    {

        switch(label)
        {

        case USYS_EM :
        {
            //Physical params in EM units
            //--------------------------
            usys->n  = 0.925195985520347;
            usys->ns = 1.0-usys->n;
            usys->ni = 1.0;
            usys->as = 388.81114;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
            break;
        }


        case USYS_SEM :
        {
            //Physical params in SEM units
            //--------------------------
            usys->ns = 1.0;
            usys->ni = 1.0/ns_EM;
            usys->n  = usys->ni-usys->ns;
            usys->as = 1.0;
            usys->ai = 1.0/as_EM;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = ms_EM/(1.0+ms_EM);
            usys->me = me_EM/(1.0+ms_EM);
            usys->mm = mm_EM/(1.0+ms_EM);
            break;
        }

        default :  //Earth-Moon
        {
            cout << "init_USYS. Unrecognised unit system. Earth-Moon units are used by default. " << endl;
            //Physical params in EM units
            //--------------------------
            usys->n  = 0.925195985520347;
            usys->ns = 1.0-usys->n;
            usys->ni = 1.0;
            usys->as = 388.81114;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;
            usys->ms = usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;
        }
        }

        //SEM mass ratio
        usys->mu_SEM = (usys->me+usys->mm)/(usys->ms+usys->me+usys->mm);
        //SE mass ratio
        usys->mu_SE  = (usys->me)/(usys->ms+usys->me);
        break;

    }
    case M_RTBP:
    case M_ERTBP:
    {
        switch(label)
        {

        //SEM mass ratio
        usys->mu_SEM = +3.040277404912506e-06;
        //SE mass ratio
        usys->mu_SE  = usys->mu_EM*usys->mu_SEM/(1.0 + usys->mu_EM);

        case USYS_EM:
            usys->n  = 0.925195985520347;
            usys->ns = 1.0-usys->n;
            usys->ni = 1.0;
            usys->as = 388.81114;
            usys->ai = 1.0;
            usys->T  = 2*M_PI/usys->n;

            usys->ms = 0.0; //usys->ns*usys->ns*pow(usys->as,3.0)-1;
            usys->me = 1.0-usys->mu_EM;
            usys->mm = usys->mu_EM;

            //Physical params in EM units
            //--------------------------
            //            usys->ns  = 0.0;
            //            usys->ni  = 1.0;
            //            usys->n   = 1.0;
            //            usys->as  = 0.0;
            //            usys->ai  = 1.0;
            //            usys->T   = 2*M_PI;  //Note the choice is made to take T = 2pi here
            //            usys->ms  = 0.0;
            //            usys->me  = 1.0-usys->mu_EM;
            //            usys->mm  = usys->mu_EM;

            break;

        case USYS_SEM:

            //Physical params in SEM units
            //--------------------------
            usys->ns = 1.0;
            usys->ni = 1.0/ns_EM;
            usys->n  = usys->ni-usys->ns;
            usys->as = 1.0;
            usys->ai = 1.0/as_EM;
            usys->T  = 2*M_PI/usys->n;

            usys->ms = 1.0-usys->mu_SEM;
            usys->me = usys->mu_EM;
            usys->mm = 0.0;//mm_EM/(1.0+ms_EM);

            //Physical params in SEM units
            //--------------------------
            //usys->ns  = 1.0;
            //usys->ni  = 0.0;
            //usys->n   = usys->ni-usys->ns;
            //usys->as  = 1.0;
            //usys->ai  = 0.0;
            //usys->T   = 2*M_PI;  //Note the choice is made to take T = 2pi here
            //usys->ms  = 1.0-usys->mu_SEM;
            //usys->me  = usys->mu_EM;
            //usys->mm  = 0.0;

            break;

        default:  //Earth-Moon
            cout << "init_USYS. Unrecognised unit system. Earth-Moon units are used by default. " << endl;
            //Physical params in EM units
            //--------------------------
            usys->ns  = 0.0;
            usys->ni  = 1.0;
            usys->n   = 1.0;
            usys->as  = 0.0;
            usys->ai  = 1.0;
            usys->T   = 2*M_PI;  //Note the choice is made to take T = 2pi here
            usys->ms  = 0.0;
            usys->me  = 1.0-usys->mu_EM;
            usys->mm  = usys->mu_EM;
        }
    }
    }


}

/**
 *  \brief From EM to NC coordinates for the primaries. Used in qbtbp_ofs_fft_*
 */
void EMtoNC_prim(double Zc[3], double zc[3], double c1, double gamma)
{
    zc[0] = c1 - Zc[0]/gamma;
    zc[1] =    - Zc[1]/gamma;
    zc[2] =    + Zc[2]/gamma;
}

/**
 * \brief Initializes a coordinate systems (CSYS structure), with associated vector field coefficients, data folder names, and unit system.
 * \param csys pointer on the CSYS structure to initialize.
 * \param qbcp_l pointer on the QBCP_L structure that contains csys.
 * \param qbcp pointer on the QBCP structure that contains parameters specific to each libration points (namely, gamma)
 * \param fwrk indix of the coordinate system to use (F_EM, F_SEM).
 * \param li number of the libration point to focus on (L1, L2).
 * \param coefNumber the number of vector field coefficients to initialize. It has been set in the QBCP_init function.
 * \param isNew boolean. if true, the qbtbp has not been computed via the qbtbp() routine, so the vector field coefficients cannot be initialized.
 *
 *   Note that the QBCP structure is used only for the initialization of the coordinate systems. More precisely, it contains some parameters
 *   specific to each libration point (gamma), via its CR3BP structures.
 **/
void init_CSYS(CSYS *csys, QBCP_L *qbcp_l, QBCP *qbcp, int fwrk, int li, int coefNumber, int isNew, int pmType, int manType)
{
    int nf = qbcp_l->nf;
    int model = qbcp_l->model;

    //--------------------------------------
    //Complete folders
    //--------------------------------------
    csys->F_GS     = init_F_FOLDER("data/PMFFT",  model, fwrk, li);     //Graph style (PM)
    csys->F_NF     = init_F_FOLDER("data/NF",     model, fwrk, li);     //Normal form style(PM)
    csys->F_MS     = init_F_FOLDER("data/MS",     model, fwrk, li);     //Mixed style (PM)

    csys->F_CS     = init_F_FOLDER("data/CS",     model, fwrk, li);     //Center-stable (PM)
    csys->F_CU     = init_F_FOLDER("data/CU",     model, fwrk, li);     //Center-unstable (PM)
    csys->F_CUS    = init_F_FOLDER("data/CUS",    model, fwrk, li);     //Center-hyperbolic (PM)

    csys->F_COEF   = init_F_FOLDER("data/VF",     model, fwrk, li);     //For integration in a given coord. system
    csys->F_COC    = init_F_FOLDER("data/COC",    model, fwrk, li);     //For the change of coordinates of the PM
    csys->F_PLOT   = init_F_FOLDER("plot",        model, fwrk, li);     //For plot output (gnuplot)
    csys->F_PRINT  = init_F_FOLDER("fprint",      model, fwrk, li);     //For print output (data used postprocessed in R)

    //--------------------------------------
    //Parameterization folders
    //--------------------------------------
    csys->manType = manType;
    switch(manType)
    {
    case MAN_CENTER:
        //If the Center Manifold is selected,
        //we can choose between the styles
        switch(pmType)
        {
        case PMS_GRAPH:
            csys->F_PMS = csys->F_GS;
            break;
        case PMS_NORMFORM:
            csys->F_PMS = csys->F_NF;
            break;
        case PMS_MIXED:
            csys->F_PMS = csys->F_MS;
            break;

        }
        break;

    case MAN_CENTER_S:
        csys->F_PMS = csys->F_CS;
        break;

    case MAN_CENTER_U:
        csys->F_PMS = csys->F_CU;
        break;

    case MAN_CENTER_US:
        csys->F_PMS = csys->F_CUS;
        break;
    }

    //--------------------------------------
    // Unit system associated with csys
    //--------------------------------------
    switch(fwrk)
    {
    case F_EM:
        csys->us = qbcp_l->us_em;
        break;
    case F_SEM:
    case F_SE:
        csys->us = qbcp_l->us_sem;
        break;
    default:
        cout << "init_CSYS. Warning: unknown model." << endl;
    }

    //--------------------------------------
    // c1 and gamma, for normalized computation
    //--------------------------------------
    // First, select the right framework
    // Gives the assosicate CR3BP and mu
    //------------------------------
    CR3BP cr3bp_root;
    switch(fwrk)
    {
    case F_EM:
        cr3bp_root  = qbcp->cr3bp1;
        csys->cr3bp = qbcp->cr3bp1;
        csys->mu    = qbcp_l->us_em.mu_EM;
        break;
    case F_SEM:
    case F_SE:
        cr3bp_root  = qbcp->cr3bp2;
        csys->cr3bp = qbcp->cr3bp2;
        csys->mu          = qbcp_l->us_sem.mu_SEM;
        break;
    default:
        cout << "init_CSYS. Warning: unknown model." << endl;
        cr3bp_root  = qbcp->cr3bp1; //cr3bp_root is set here to avoid compilation warning, but the default case should NOT be used!
        csys->cr3bp = qbcp->cr3bp1;
        csys->mu    = qbcp_l->us_em.mu_EM;
    }

    //------------------------------
    //Then set the value of c1 and gamma
    //according to the selected libration point
    //------------------------------
    csys->li = li;
    switch(li)
    {
    case 1:
    {
        csys->gamma = cr3bp_root.l1.gamma_i;
        csys->c1    = (csys->mu-1+csys->gamma)/csys->gamma;
        break;
    }

    case 2:
    {
        csys->gamma =  cr3bp_root.l2.gamma_i;
        csys->c1    = (csys->mu-1-csys->gamma)/csys->gamma;
        break;
    }

    case 3: //same as L1/L2. Different from convention by Jorba & Masdemont 1999
    {
        csys->gamma =  cr3bp_root.l3.gamma_i;
        csys->c1    = (csys->mu + csys->gamma)/csys->gamma;
        break;
    }
    }


    //--------------------------------------
    // Creating the arrays of coefficients
    //--------------------------------------
    csys->coeffs = (double*) calloc(coefNumber*(qbcp_l->nf+1), sizeof(double)); //Default set of vector field coefficients
    csys->Ps  = (double*) calloc(3*(nf+1), sizeof(double));  //Sun   position in EM coordinates
    csys->Pe  = (double*) calloc(3*(nf+1), sizeof(double));  //Earth position in EM coordinates
    csys->Pm  = (double*) calloc(3*(nf+1), sizeof(double));  //Moon  position in EM coordinates
    csys->ps  = (double*) calloc(3*(nf+1), sizeof(double));  //Sun   position in NC coordinates
    csys->pe  = (double*) calloc(3*(nf+1), sizeof(double));  //Earth position in NC coordinates
    csys->pm  = (double*) calloc(3*(nf+1), sizeof(double));  //Moon  position in NC coordinates

    //--------------------------------------
    // Retrieving the arrays of coefficients
    //--------------------------------------
    if(isNew) //the QBCP is new, no coefficients exist
    {
        cout << "init_CSYS. The considered QBCP is new:" << endl;
        cout << "Integration coefficient are not retrieved from existing files." << endl;
    }
    else    //the coefficients can be retrieved
    {
        //cout << "init_CSYS. Integration coefficient are retrieved from existing files." << endl;
        //The flag compType is here to tell if we want the coefficients
        //computed from the FFT or from direct computation
        double compType;

        //Switch between the models
        switch(model)
        {
        case M_QBCP:
        case M_BCP:
        case M_ERTBP:
        {
            compType = 1; //from FFTs
            //-------------------------
            // Default set of vector field coefficients
            //-------------------------
            coefRetrieving(csys->F_COEF+"alpha", csys->coeffs, nf, 0, compType, coefNumber);
            //-------------------------
            // primaries position
            //-------------------------
            coefRetrieving(csys->F_COEF+"Ps", csys->Ps, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"Pe", csys->Pe, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"Pm", csys->Pm, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"ps", csys->ps, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"pe", csys->pe, nf, 0, compType, 3);
            coefRetrieving(csys->F_COEF+"pm", csys->pm, nf, 0, compType, 3);
            break;

        }

        case M_RTBP:
        {
            //cout << "init_CSYS. The use of the RTBP has been detected." << endl;
            //----------------------------------------
            // Default set of vector field coefficients
            // A REVOIR: METTRE UN SWITCH GEANT POUR TOUS LES CAS
            // OU: LES METTRE DANS DES FICHIERS TXT via une routine type bcp
            //----------------------------------------
            csys->coeffs[0] = 1.0;
            csys->coeffs[1] = 0.0;
            csys->coeffs[2] = 1.0;
            csys->coeffs[3] = 0.0;
            csys->coeffs[4] = 0.0;
            csys->coeffs[5] = 1.0;
            //Sun position
            csys->coeffs[6] = 0.0;
            csys->coeffs[7] = 0.0;
            //Earth position
            csys->coeffs[8] = csys->mu;
            csys->coeffs[9] = 0.0;
            //Moon position
            csys->coeffs[10] = csys->mu-1.0;
            csys->coeffs[11] = 0.0;
            //NC coeffs
            csys->coeffs[12] = -csys->c1;     //alpha13 = alpha[12] = -c1
            csys->coeffs[13] = 0.0;           //alpha14 = alpha[13] = 0


            //----------------------------------------
            // Earth, Moon, and Sun €€TODO: a revoir!
            // Mettre dans un fichier txt...
            //----------------------------------------
            csys->Pe[0] = csys->mu;
            csys->Pe[1] = 0.0;
            csys->Pe[2] = 0.0;

            csys->Pm[0] = csys->mu-1.0;
            csys->Pm[1] = 0.0;
            csys->Pm[2] = 0.0;

            csys->Ps[0] = 0.0;
            csys->Ps[1] = 0.0;
            csys->Ps[2] = 0.0;

            EMtoNC_prim(csys->Pe, csys->pe, csys->c1, csys->gamma);
            EMtoNC_prim(csys->Pm, csys->pm, csys->c1, csys->gamma);
            EMtoNC_prim(csys->Ps, csys->ps, csys->c1, csys->gamma);

            break;
        }

        default:
            cout << "init_QBCP. Warning: unknown model." << endl;
        }

    }


    //--------------------------------------
    // Storing the solutions of the QBTBP
    //--------------------------------------
    csys->zt = Ofsc(nf);
    csys->Zt = Ofsc(nf);
    csys->ztdot = Ofsc(nf);
    csys->Ztdot = Ofsc(nf);
    switch(model)
    {
    case M_RTBP:
    case M_BCP:
    {
        //Perfect circles
        csys->zt.setCoef(1.0+0.0*I, 0);
        csys->Zt.setCoef(1.0+0.0*I, 0);
        break;
    }


    case M_QBCP:
    default:
    {
        //Taken from files
        string filename = "data/qbtbp/";
        readOFS_txt(csys->zt, filename+"bjc");
        readOFS_txt(csys->Zt, filename+"cjc");
        //Derivatives
        csys->ztdot.dot(csys->zt, csys->us.n);
        csys->Ztdot.dot(csys->Zt, csys->us.n);
        break;
    }
    }
}

/**
* \brief Initialize the Quasi-Bicircular Four-Body Problem in the form of a QBCP structure.
* \param qbcp pointer on the QBCP structure to init.
* \param n1 name of the first primary (e.g Sun)
* \param n2 name of the second primary (e.g. Earth)
* \param n3 name of the third primary  (e.g. Moon)
*
* NEEDS TO BE MODIFIED TO GET RID OF THE HARD CODED VALUES
**/
void init_QBCP(QBCP *qbcp, int BIG, int MEDIUM, int SMALL, int model)
{
    //E.g. Earth-Moon init
    init_CR3BP(&qbcp->cr3bp1, MEDIUM, SMALL);
    //E.g. Sun-Earth+Moon init
    //WARNING: in the case of the Sun-Earth-Moon system, we need to set EARTH_AND_MOON, instead of EARTH alone.
    //In fact, ALL Four-body configurations should be built that way, so
    //keep in mind that if the configuration is not SUN/EARTH/MOON, the result will be biased.
    if(MEDIUM ==  EARTH && SMALL == MOON) init_CR3BP(&qbcp->cr3bp2, BIG, EARTH_AND_MOON);
    else init_CR3BP(&qbcp->cr3bp2, BIG, MEDIUM);
}

/**
 *  \brief Initializes two QBCP_L with two different models for continuation process.
 **/
void init_QBCP_I(QBCP_I *model,
                 QBCP_L *model1,
                 QBCP_L *model2,
                 int n1, int n2, int n3,
                 int isNormalized, int li_EM, int li_SEM,
                 int isNew,
                 int mod1,
                 int mod2,
                 int frwk,
                 int pmStyle)
{
    //Initialize the models
    QBCP qbp1, qbp2;
    init_QBCP(&qbp1, n1, n2, n3, mod1);
    init_QBCP(&qbp2, n1, n2, n3, mod2);

    //Initialize the models around the given li point
    init_QBCP_L(model1, &qbp1, isNormalized, li_EM, li_SEM, isNew, mod1, frwk, pmStyle, MAN_CENTER, MAN_CENTER);
    init_QBCP_L(model2, &qbp2, isNormalized, li_EM, li_SEM, isNew, mod2, frwk, pmStyle, MAN_CENTER, MAN_CENTER);

    //Store in model
    model->model1 = *model1;
    model->model2 = *model2;
    model->epsilon = 0.0;
}


/**
 * \brief Initializes a libration point.
 * \param libp a pointer towards the LibrationPoint structure to init.
 * \param cr3bp a CR3BP structure that contains useful coefficients.
 * \param number the indix of the libration point to init.
 **/
void init_libp(LibrationPoint *libp, CR3BP cr3bp, int number)
{

    double gamma_i;

    //Number
    libp->number = number;

    switch(number)
    {
    case 1:
        //Gamma
        gamma_i = cr3bp.rh - 1.0/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);                    //initial guess
        gamma_i = rtnewt(polynomialLi, gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);   //newton-raphson method
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu - gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;

        //printf("l1 = (%5.10f, 0, 0)\n", libp->position[0]);
        break;

    case 2:
        //Gamma
        gamma_i = cr3bp.rh + 1.0/3.0*pow(cr3bp.rh,2.0)- 1/9*pow(cr3bp.rh,3);
        gamma_i = rtnewt(polynomialLi, gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);
        libp->gamma_i = gamma_i;

        //Position
        libp->position[0] = 1 - cr3bp.mu + gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;

        //printf("l2 = (%5.10f, 0, 0)\n", libp->position[0]);
        break;

    case 3:
        //Gamma
        gamma_i = 7/12.0*cr3bp.mu + pow(237,2.0)/pow(12,4.0)*pow(cr3bp.mu,3.0);
        gamma_i = rtnewt(polynomialLi, gamma_i, LIBRATION_POINT_PRECISION, cr3bp.mu, number);
        libp->gamma_i = 1-gamma_i;  //BEWARE: for L3, gamma3 = L3-M1 distance != L3-M2


        //Position
        libp->position[0] = - cr3bp.mu - libp->gamma_i;
        libp->position[1] = 0;
        libp->position[2] = 0;

        //printf("l3 = (%5.10f, 0, 0)\n", libp->position[0]);
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
* \param body a pointer on the Body structure to init.
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
        body->Rm  = 6371.00;         //[km]
        body->M   = 5.97219e24;       //[kg]
        body->GM  = 398600.440;      //[km^3.s^-2]

        //Orbital parameters
        body->a = 149.60e6;          //[km]
        body->T = 365.25636*days;     //[s]

        strcpy(body->name, "Earth");
        break;

    case MOON:

        //Physical parameters
        body->Req = 1737.5;       //[km]
        body->Rm = 1737.5;        //[km]
        body->M =  0.07345814120628661e24;    //[kg] TO BE CONSISTENT WITH HARD CODED VALUE OF mu(Earth-Moon)
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

//-----------------------------------------------------------------------------------------------------------------------------
//            Change Coord. System
//-----------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Change the default coordinate system
 **/
void changeDCS(QBCP_L &qbcp_l, int fwrk)
{
    qbcp_l.fwrk = fwrk;
    switch(fwrk)
    {
    case F_EM:
        qbcp_l.us = qbcp_l.us_em;
        qbcp_l.cs = qbcp_l.cs_em;
        qbcp_l.li = qbcp_l.li_EM;
        break;
    case F_SEM:
        qbcp_l.us = qbcp_l.us_sem;
        qbcp_l.cs = qbcp_l.cs_sem;
        qbcp_l.li = qbcp_l.li_SEM;
        break;
    case F_SE:
        qbcp_l.us = qbcp_l.us_sem;
        qbcp_l.cs = qbcp_l.cs_sem;
        qbcp_l.li = qbcp_l.li_SEM;
        break;
    default:
        cout << "changeDCS. Warning: unknown fwrk." << endl;
    }
}

/**
 *  \brief Change the default coordinate system and the libration point for this coordinate system
 **/
void changeLIDCS(QBCP_L &qbcp_l, int fwrk, int li)
{
    //Default settings
    qbcp_l.fwrk = fwrk;  //new default cs
    qbcp_l.li = li;      //new default libration point

    //Change the coord. system approprietly: "fwrk" around "li"
    switch(fwrk)
    {
    case F_EM:
        switch(li)
        {
        case 1:
            qbcp_l.cs_em  = qbcp_l.cs_em_l1;
            break;
        case 2:
            qbcp_l.cs_em  = qbcp_l.cs_em_l2;
            break;
        case 3:
            qbcp_l.cs_em  = qbcp_l.cs_em_l3;
            break;
        }
        qbcp_l.us = qbcp_l.us_em;
        qbcp_l.cs = qbcp_l.cs_em;
        break;
    case F_SEM:
        switch(li)
        {
        case 1:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l1;
            break;
        case 2:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l2;
            break;
        case 3:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l3;
            break;
        }
        qbcp_l.us = qbcp_l.us_sem;
        qbcp_l.cs = qbcp_l.cs_sem;
        break;
    case F_SE:
        switch(li)
        {
        case 1:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l1;
            break;
        case 2:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l2;
            break;
        case 3:
            qbcp_l.cs_sem  = qbcp_l.cs_sem_l3;
            break;
        }
        qbcp_l.us = qbcp_l.us_sem;
        qbcp_l.cs = qbcp_l.cs_sem;
        break;
    default:
        cout << "changeLIDCS. Warning: unknown fwrk." << endl;
    }
}


//-----------------------------------------------------------------------------------------------------------------------------
//            Subroutines
//-----------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Return the string corresponding to the libration point number provided (e.g. "L1" if li == 1).
 **/
string init_F_LI(int li)
{
    switch(li)
    {
    case 1:
        return "L1";
    case 2:
        return "L2";
    case 3:
        return "L3";
    case 4:
        return "L4";
    case 5:
        return "L5";
    default:
        cout << "init_F_LI. Warning: supplied libration number is greater than 5." << endl;
    }
    return "L1"; //never here
}

/**
 *  \brief Return the string corresponding to the model indix provided (e.g. "QBCP" if model == M_QBCP).
 **/
string init_F_MODEL(int model)
{
    switch(model)
    {
    case M_QBCP:
        return "QBCP";
    case M_BCP:
        return "BCP";
    case M_RTBP:
        return "RTBP";
    case M_ERTBP:
        return "ERTBP";
    default:
        cout << "init_F_MODEL. Warning: unknown model." << endl;
    }
    return "QBCP"; //never here
}

/**
 *  \brief Return the string corresponding to the framework (coord. syst.) indix provided (e.g. "EM" if fwrk == F_EM).
 **/
string init_F_FWRK(int fwrk)
{
    switch(fwrk)
    {
    case F_EM:
        return  "EM";
    case F_SEM:
        return  "SEM";
    case F_SE:
        return  "SE";
    default:
        cout << "init_F_FWRK. Warning: unknown model." << endl;
    }

    return "EM"; //never here
}

/**
 *  \brief Return the folder name corresponding to the prefix/model/framework/libration point number combination provided (e.g. "prefix/QBCP/EM/L1").
 **/
string init_F_FOLDER(string prefix, int model, int fwrk, int li)
{
    return prefix+"/"+init_F_MODEL(model)+"/"+init_F_FWRK(fwrk)+"/"+init_F_LI(li)+"/";
}

/**
 * \brief Retrieve a set of coefficients, given as Fourier series from a txt file.
 * \param filename the name of the txt file.
 * \param params a pointer toward the array to update.
 * \param nf the order of the Fourier series.
 * \param shift the indix from which to start the storage of the coefficients in params.
 * \param flag: if flag == 1, the coefficients computed via FFT are used. Otherwise, the expansions obtained through Fourier series algebraic manipulations are used.
 *
 *  Warning: As of now, FFT coefficients must be used for betas and deltas (see QBCP_L structure).
 **/
void coefRetrieving(string filename, double *params, int nf, int shift, int flag, int number)
{
    //Reading tools
    ifstream readStream;
    double cDouble1;
    int alphaNumber = 1;
    string ss = static_cast<ostringstream*>( &(ostringstream() << alphaNumber) )->str();
    for(int header = shift; header <= (nf+1)*(number-1)+shift; header+=(nf+1))
    {
        if(flag) readStream.open((filename+ss+"c_fft.txt").c_str());
        else readStream.open((filename+ss+"c.txt").c_str());
        for(int i=0; i<= nf; i++)
        {
            readStream >> cDouble1;  //current order
            readStream >> params[i+header];
        }
        readStream.close();
        alphaNumber++;
        ss = static_cast<ostringstream*>( &(ostringstream() << alphaNumber) )->str();
    }

    if(!flag) cout << "coefRetrieving: the FFT coefficients have not been used." << endl;
}

/**
 * \brief Compute the potential energy for the given state and CR3BP (mu).
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double energy(double y[], double mu)
{
    double r1 = sqrt( pow(y[0]+mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );
    double r2 = sqrt( pow(y[0]- 1 + mu,2.0) + pow(y[1],2.0) + pow(y[2],2.0) );

    return - ( 1.0/2.0*(pow(y[0],2) + pow(y[1],2)) + (1-mu)/r1 + mu/r2 + 1.0/2.0*mu*(1-mu) );
}

/**
 * \brief Compute the Jacobi constant for the given state and CR3BP (mu)
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double jacobi(double y[], double mu)
{
    return -2.0*energy(y, mu) - (pow(y[3],2.0)+pow(y[4],2.0)+pow(y[5],2.0));
}

/**
 * \brief Compute the coefficient cn for a given libration point (L1 or L2 for now)
 * \param qbcp_l a reference to the QBCP_L initialized around the selected libration point.
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 *
 * The L3 value is taken from Richardson, 1980:
 *
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( 1 - \mu  + \frac{\mu \gamma_j^{n+1}}{(1 + \gamma_j)^{n+1}} \right) \f$.
 **/
double cn(QBCP_L& qbcp_l, int n)
{
    double gamma = qbcp_l.cs.gamma;
    double mu;
    switch(qbcp_l.fwrk)
    {
        case F_EM:
            mu   = qbcp_l.us.mu_EM;
        break;
        case F_SEM:
            mu   = qbcp_l.us.mu_SEM;
        break;
        case F_SE:
            mu   = qbcp_l.us.mu_SE;
        break;
        default: //EM by default
            cout << "WARNING in cn(): unknown framework. EM by default." << endl;
            mu   = qbcp_l.us.mu_EM;
    }


    double res = 0.0;
    switch(qbcp_l.cs.li)
    {
    case 1:
        res =  pow(gamma,-3.0)*(pow(+1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0-gamma), n+1));
        break;
    case 2:
        res =  pow(gamma,-3.0)*(pow(-1.0, n)*mu + pow(-1.0, n)*(1-mu)*pow(gamma/(1.0+gamma), n+1));
        break;
    case 3:
        res =  pow(gamma,-3.0)*(1.0 - mu + mu*pow(gamma/1.0+gamma, n+1)); //convention by Richardson 1980
        break;
    default:
        cout << "cn. Warning: supplied Li number is out of scope. 0.0 is returned." << endl;
    }
    return res;
}

/**
* \brief Using the Newton-Raphson method, find the root of a function known to lie close to x1 The root
*  /rtnewt will be refined until its accuracy is known within ± xacc.
*  funcd is a user-supplied routine that returns both the function value and the first derivative of the
*  function at the point x.
**/
double rtnewt(void (*funcd)(double, int, double, double *, double *), double x1, double xacc, double mu, int number)
{
    void nrerror(char error_text[]);
    int j;
    double df,dx,f,rtn;

    rtn=x1;   //initial guess

    for (j=1; j<=50; j++)
    {
        (*funcd)(mu, number, rtn,&f,&df);
        dx=f/df;
        rtn -= dx;

        if (fabs(dx) < xacc) return rtn;  //Convergence
    }
    printf("WARNING: Maximum number of iterations exceeded in rtnewt");
    return 0.0;   //Never get here.
}

/**
* \brief Provides the function value and its first derivative for the newton-raphson method.
* f corresponds to the equation satisfied by the Li-m2 distance for the L1/L2 cases
* and by 1-(Li-m1 distance) for the L3 case
**/
void polynomialLi(double mu, int number, double y, double *f, double *df)
{
    switch(number)
    {

    case 1:
        *f =  pow(y,5.0)   - (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) +  2*mu*y - mu;
        *df = 5*pow(y,4.0) - 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    +  2*mu;
        break;

    case 2:
        *f =  pow(y,5.0)   + (3.0-mu)*pow(y,4.0) + (3-2*mu)*pow(y,3.0) - mu*pow(y,2.0) -  2*mu*y - mu;
        *df = 5*pow(y,4.0) + 4*(3.0-mu)*pow(y,3.0) + 3*(3-2*mu)*pow(y,2.0) - 2*mu*y    -  2*mu;
        break;

    case 3:
        //*f =  pow(y,5.0) + (2.0+mu)*pow(y,4.0) + (1+2*mu)*pow(y,3.0) + (1+mu)*pow(y,2.0) +  2*(1-mu)*y + 1-mu;
        *f= pow(y,5.0) + (7+mu)*pow(y,4.0) + (19+6*mu)*pow(y,3.0) -(24+13*mu)*pow(y,2.0) +  (12+14*mu)*y -7*mu;
        *df= 5*pow(y,4.0) + 4*(7+mu)*pow(y,3.0) + 3*(19+6*mu)*pow(y,2.0) -2*(24+13*mu)*pow(y,1.0) +  (12+14*mu);
        //*df = 5*pow(y,4.0) + 4*(2.0+mu)*pow(y,3.0) + 3*(1+2*mu)*pow(y,2.0) + 2*(1+mu)*y +  2*(1-mu);
        break;
    }
}
