#include "Ofsc.h"

int OFS_ORDER;

/**
 * \file ofs.cpp
 * \brief Fourier series template class (src)
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Ofsc.
 */
Ofsc::Ofsc()
{
    order = OFS_ORDER;
    coef = new cdouble[2*order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor with a given order.
 */
Ofsc::Ofsc(const int newOrder)
{
    order = newOrder;
    coef = new cdouble[2*order+1];
    this->zero(); //all coefficients to zero
}

/**
 *  \brief Constructor from a given Ofsc object.
 */
Ofsc::Ofsc(Ofsc const& b)
{
    order = b.order;
    coef = new cdouble[2*order+1];
    for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order]; //this->setCoef(b.getCoef(i), i);
}

//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief  Copy from a given Ofs object (only the coefficients).
 */
Ofsc& Ofsc::ccopy(Ofsc const& b)
{
    if(order != b.order)
    {
        cout << "Erreur in Ofsc::ccopy: orders do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {
        for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order];
        return *this;
    }
}


//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Ofsc.
 */
Ofsc::~Ofsc()
{
    if(coef != NULL) delete[] coef;
    coef = 0;
}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets a coefficient at a given position in the serie.
 */
void Ofsc::setCoef(cdouble const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] = value;
    else cout << "Error in Ofsc::setCoef: position is out of scope. No coefficient is set." << endl;
}


/**
 *  \brief Adds a coefficient at a given position in the serie.
 */
void Ofsc::addCoef(cdouble const& value, int const& pos)
{
    if(fabs(pos) <= order) coef[pos+order] += value;
    else cout << "Error in Ofsc::addCoef: position is out of scope\n" << endl;
}


//---------------------------------------------------------------------------
//Operators (+=, -=, ...)
//---------------------------------------------------------------------------
/**
 *  \brief  An operator. Constructor from a given Ofsc object (only the coefficients).
 */
Ofsc& Ofsc::operator = (Ofsc const& b)
{
    if(this != &b)
    {
        order = b.order;
        if(coef != NULL) delete coef;
        coef = new cdouble[2*order+1];
        for(int i = -order ; i<= order; i++) coef[i+order] = b.coef[i+order]; //this->setCoef(b.getCoef(i), i);
    }
    return *this; //same object if returned
}

/**
 *  \brief  An operator. Adds all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.order != order.
 */
Ofsc& Ofsc::operator += (Ofsc const& b)
{
    if(b.order > order) //if b.order > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        cdouble temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new cdouble[2*b.order+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+b.order] = temp[i+order];
        order = b.order;
    }

    //Adding the coefficients
    for(int i = -b.order ; i <= b.order ; i++) coef[i+order] += b.coef[i+b.order];
    return *this;
}

/**
 *  \brief  An operator. Subtracts all coefficients term by term  from a given Ofs object.
 *
 *   Allows b.order != order.
 */
Ofsc& Ofsc::operator -= (Ofsc const& b)
{
    if(b.order > order) //if b.order > order, a new array of coefficients must be set
    {
        //Copy coef into temporary array
        cdouble temp[2*order+1];
        for(int i = 0 ; i< 2*order + 1; i++) temp[i] = coef[i];
        //Recreate a good array
        delete coef;
        coef = new cdouble[2*b.order+1];
        //Store the coefficients again
        for(int i = -order ; i<= order; i++) coef[i+b.order] = temp[i+order];
        order = b.order;
    }

    //Adding the coefficients
    for(int i = -b.order ; i <= b.order ; i++) coef[i+order] -= b.coef[i+b.order];
    return *this;
}

/**
 *  \brief  An operator. Multiplies all coefficients by a given \c cdouble coefficient.
 */
Ofsc& Ofsc::operator *= (cdouble const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] *= c;
    return *this;
}

/**
 *  \brief  An operator. Divides all coefficients by a given \c cdouble coefficient.
 */
Ofsc& Ofsc::operator /= (cdouble const& c)
{
    for(int i=0; i<2*order+1; i++) coef[i] /= c;
    return *this;
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the serie.
 */
int Ofsc::getOrder() const
{
    return order;
}

/**
 *  \brief  Gets the pointer address of the Ofsc object
 */
Ofsc* Ofsc::getAddress() const
{
    return (Ofsc*) this;
}

/**
 *  \brief  Gets the coefficient at a given position. cdouble case
 */
cdouble Ofsc::getCoef(int pos) const
{
    if(fabs(pos) <= order)  return coef[pos+order];
    else
    {
        cout << "Warning in Ofsc::getCoef: position " << pos << " is out of scope. 0.0 is returned by default." << endl;
        return  0.0+0.0*I;
    }
}


//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------

/**
 *  \brief  Sets all coefficients to zero. cdouble case
 */
void Ofsc::zero()
{
    for(int i = -order ; i<= order; i++) coef[i+order] = 0.0+0.0*I;
}

/**
 *  \brief  Is the Ofsc object equal to zero, at order ofs_order?
 */
bool Ofsc::isnull(const int ofs_order) const
{
    for(int i = -min(ofs_order, order); i <= min(ofs_order, order); i++)
    {
        if(cabs(getCoef(i)) != 0.0) return false;
    }
    return true;
}


//---------------------------------------------------------------------------
// Functions (evaluate)
//---------------------------------------------------------------------------
/**
 *  \brief  Evaluates the Ofsc object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
 cdouble Ofsc::fevaluate(double cR[], double sR[], int eff_order)
{
    cdouble result = 0+0.0*I;
    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];
    //Return result
    return result;
}

/**
 *  \brief  Evaluates the Ofsc object at angle theta (theta = nt) and at a certain order eff_order.
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
 cdouble Ofsc::evaluate(double const& theta, int eff_order)
{
    cdouble result = 0+0.0*I;
    double cR[eff_order];
    double sR[eff_order];

    cR[0] = cos(theta);
    sR[0] = sin(theta);
    for(int i = 1; i< eff_order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=eff_order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];//getCoef(0);


    return result;
}

/**
 *  \brief  Evaluates the Ofsc object at angle theta (theta = nt).
 *
 *  The sum is given in the form: \f$ c_0 + \sum \limits_{k = 1}^{J} \cos(k\theta) (c_k + c_{-k}) + i\sin(k\theta) (c_k - c_{-k}) \f$
 *  in order to avoid any additionnal roundoff errors by using cos and sin functions as little as possible.
 */
 cdouble Ofsc::evaluate(double const& theta)
{
    cdouble result = 0+0.0*I;
    double cR[order];
    double sR[order];

    cR[0] = cos(theta);
    sR[0] = sin(theta);
    for(int i = 1; i< order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    //Order -k and +k are calculated together to avoid any additionnal roundoff errors
    for(int k=order; k>=1; k--) result+= cR[k-1]*(coef[k+order] + coef[-k+order]) + I*sR[k-1]*(coef[k+order] - coef[-k+order]);
    //Order 0
    result+= coef[0+order];//getCoef(0);

    return result;
}

//---------------------------------------------------------------------------
// Functions (operations)
//---------------------------------------------------------------------------
/**
 *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ at a certain order eff_order.
 *
 *  Note: can be used in place.
 */
void Ofsc::ofs_smult(Ofsc const& a, cdouble c, int eff_order)
{
    //Sum
    for(int i = -eff_order; i <= eff_order; i++)
    {
        addCoef(c*a.getCoef(i), i);
    }
}


/**
 *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$.

    Notes:
    1. \c this \f$ += a \times b \f$. with new order = max(order, a.order, b.order): the sum is truncated at max(a.getOrder(),b.getOrder()).
    2. WARNING: Need improvement: for n <= a.getOrder(), b.getOrder(), some products are out of scope in: this->addCoef(a.getCoef(p)*b.getCoef(n-p), n). The getCoef function set these coefficients to zero, which guarantees the good result. However, unecessary product are made. psup and pinf must be redefined.
    3. Works fine when a.order = b.order which is the default case.
 */
void Ofsc::ofs_sprod(Ofsc const& a, Ofsc const& b)
{
    int psup, pinf;
    //Product
    for(int n=-order ; n<= order; n++)
    {
        psup = min(n+order,  order);
        pinf = max(n-order, -order);
        for(int p=pinf; p<= psup; p++) coef[n+order] += a.coef[p+order]*b.coef[n-p+order];
    }
}

/**
 *  \brief  An operation. Sets the sum-product: \c this \f$  = m_a a + m_b b \f$.
 *
 *  Note: can be used in place.
 */
void Ofsc::ofs_fsum(Ofsc const& a, cdouble const& ma, Ofsc const& b, cdouble const& mb)
{
    if(order != a.order ||  order != b.order)
    {
        cout << "Error using fsum: the order does not match. Initial Ofsc is returned" << endl;
    }
    else
    {
        for(int i=-order; i<=order; i++)
        {
            coef[i+order] = ma*a.coef[i+order] + mb*b.coef[i+order];//setCoef(ma*a.getCoef(i)+mb*b.getCoef(i), i);
        }
    }
}


//---------------------------------------------------------------------------
//Print
//---------------------------------------------------------------------------
/**
 *  \brief  A stream operator
 */
std::ostream& operator << (std::ostream& stream, Ofsc const& ofs)
{
    //Coefficients
    for(int i = 0 ; i< 2*ofs.order + 1; i++)
    {
        stream << setw(3) << setiosflags(ios::right) << std::showpos << i-ofs.order << "   " <<  setiosflags(ios::scientific) << setprecision(15) << creal(ofs.coef[i]) << "  " << cimag(ofs.coef[i]) << endl;

    }
    return stream;
}


//---------------------------------------------------------------------------
// Reading an OFS from a text file
//---------------------------------------------------------------------------
/**
 * \fn void inline readOFS_txt(Ofsc& xFFT, string filename, int fftN)
 * \brief Reading an Ofsc object from a text file.
 */
void readOFS_txt(Ofsc& xFFT, string filename)
{
    //Init
    ifstream readStream;
    double ct, cr, ci;
    int fftN = xFFT.getOrder();

    //Reading
    readStream.open((filename+".txt").c_str());
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.setCoef(cr+I*ci, i);
    }
    readStream.close();
}

//---------------------------------------------------------------------------
//Derivation
//---------------------------------------------------------------------------
/**
 *  \brief  An operation. Set the time derivative of object \c a with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
 */
void Ofsc::dot(Ofsc const& a, double const& n)
{
    //d(a_k)/dt = k*n*I*a_k
    for(int k=-order; k<= order; k++) coef[k+order] = k*n*I*a.coef[k+order]; //this->setCoef(k*n*I*a.getCoef(k), k);
}
