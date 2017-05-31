#ifndef OFS_H
#define OFS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include "parameters.h"
#include <math.h>

using namespace std;

/**
 *\class Ofsc
 *  \brief Fourier series template class.
 *
 *  The Ofsc class computes the operations of fourier series. Note that the array containing the coefficients
 *  is indexed from 0 to 2*order of the Fourier series, with the positions 0 to order-1 containing the terms
 *  of order -order to -1, and the positions order to 2*order containing the terms of order 0 to order.
 *  Ideally, the getters and the setters should be the ONLY routines to take into account this shift.
 */
class Ofsc
{
private:
    int order;  /// order of the expansion (from -order to +order) !
    cdouble *coef;    /// array of coefficients

public:
    //------------------
    //Create
    //------------------
    /**
     *  \brief Default constructor of the class Ofsc.
     */
    Ofsc();

    /**
     *  \brief Constructor with a given order.
     *  \param order_: order of the serie
     */
    Ofsc(const int order_);

    /**
     *  \brief Constructor from a given Ofsc object.
     *  \param ofs_:  a reference to the Ofsc object to copy in the new object
     */
    Ofsc(Ofsc const& ofs_);

    //------------------
    //Delete
    //------------------
    /**
     *  \brief Default destructor of the class Ofsc.
     */
    ~Ofsc();

    //------------------
    //Copy
    //------------------
    /**
     *  \brief  Copy from a given Ofs object (only the coefficients).
     *  \param  ofs_: a reference to the Ofs object to copy
     *  \return a reference to the current object
     */
    Ofsc& ccopy(Ofsc const& ofs_);

    //------------------
    //Setters
    //------------------
    /**
     *  \brief Sets a coefficient at a given position in the serie.
     *  \param value: the value to set
     *  \param pos: the position (indix in the serie) to modify
     */
    void setCoef(cdouble  const&  value, int const& pos);

    /**
     *  \brief Adds a coefficient at a given position in the serie.
     *  \param value: the value to add
     *  \param pos: the position (indix in the serie) to modify
     */
    void addCoef(cdouble  const&  value, int const& pos);


    //------------------
    //Getters
    //------------------
    /**
     *  \brief  Gets the order of the serie.
     *  \return the order of the serie as an \c int
     */
    int getOrder() const;

    /**
     *  \brief  Gets the pointer address of the Ofsc object
     *  \return A pointer to the Ofsc object
     */
    Ofsc* getAddress() const;

    /**
     *  \brief  Gets the coefficient at a given position.
     *  \param  pos: the position to get
     *  \return the coefficient of type \c cdouble at the position \c pos
     */
    cdouble getCoef(int pos) const;

    //------------------
    //Operators
    //------------------
    /**
     *  \brief  An operator. Constructor from a given Ofsc object (only the coefficients).
     *  \param  ofs_: a reference to the Ofsc object to copy
     *  \return a reference to the current object
     */
    Ofsc& operator  = (Ofsc const& ofs_);

    /**
     *  \brief  An operator. Adds all coefficients term by term  from a given Ofs object.
     *  \param  ofs_: the Ofs object to add
     *  \return a reference to the current object
     */
    Ofsc& operator += (Ofsc const& ofs_);

    /**
     *  \brief  An operator. Subtracts all coefficients term by term  from a given Ofs object.
     *  \param  ofs_: the Ofs object to subtract
     *  \return a reference to the current object
     */
    Ofsc& operator -= (Ofsc const& ofs_);

    /**
     *  \brief  An operator. Multiplies all coefficients by a given \c cdouble coefficient.
     *  \param  c: a reference to the multiplicating coefficient.
     *  \return a reference to the current object
     */
    Ofsc& operator *= (cdouble const& c);

    /**
     *  \brief  An operator. Divides all coefficients by a given \c cdouble coefficient.
     *  \param  c : a reference to the dividing coefficient.
     *  \return a reference to the current object
     */
    Ofsc& operator /= (cdouble const& c);

    //------------------
    //Operations
    //------------------
    /**
     *  \brief  An operation. Adds the product: \c this \f$  += m a \f$ at order eff_order.
     *  \param  a: a reference to an Ofsc object
     *  \param  c: reference to an cdouble coefficient
     *  \param  eff_order: the effective order of the operation
     *  \return a reference to the current object
     */
    void ofs_smult(Ofsc const& a, cdouble c, int eff_order);

    /**
     *  \brief  An operation. Set the time derivative of object \c a with pulsation \f$ \omega = n \f$, so that \c this \f$ = \dot{a} \f$.
     *  \param  a: a reference to an Ofsc object
     *  \param  n: reference to the pulsation
     *  \return a reference to the current object
     */
    void dot(Ofsc const& a, double const& n);

    /**
     *  \brief  An operation. Adds the product: \c this \f$ += a \times b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  b: reference to an Ofs object
     *  \return a reference to the current object
     */
    void ofs_sprod(Ofsc const& a, Ofsc const& b);

    /**
     *  \brief  An operation. Sets the sum-product: \c this \f$  = m_a a + m_b b \f$.
     *  \param  a: a reference to an Ofs object
     *  \param  ma: reference to an cdouble coefficient
     *  \param  b: a reference to an Ofs object
     *  \param  mb: reference to an cdouble coefficient
     *  \return a reference to the current object
     */
    void ofs_fsum(Ofsc const& a, cdouble const& ma, Ofsc const& b, cdouble const& mb);


    //------------------
    //Zeroing
    //------------------
    /**
     *  \brief  Sets all coefficients to zero.
     */
    void zero();

    /**
     *  \brief  Is the Ofsc object equal to zero, at order ofs_order?
     */
    bool isnull(const int ofs_order) const;

    //------------------
    //Evaluate
    //------------------
    /**
     *  \brief  Evaluates the Ofsc object at time t.
     *  \param  t: a reference to a \c double t
     *  \return the evaluation \f$ F_s(t) \f$ at time t.
     */
    cdouble evaluate(double const& t);

    /**
     *  \brief  Evaluates the Ofsc object at time t.
     *  \param  t: a reference to a \c double t
     *  \param  eff_order: the effective order of the operation
     *  \return the evaluation \f$ F_s(t) \f$ at time t.
     */
    cdouble evaluate(double const& t, int eff_order);
    cdouble fevaluate(double cR[], double sR[], int eff_order);

    //------------------
    //Print
    //------------------
    /**
     *  \brief  A stream operator
     */
    friend std::ostream& operator << (std::ostream& stream, Ofsc const& ofs);
};

//------------------
//Read
//------------------
/**
 * \fn void inline readOFS_txt(Ofsc& xFFT, string filename, int fftN)
 * \brief Reading an Ofsc object from a text file.
 * \param xFFT: a reference to the Ofsc object to update.
 * \param filename: a \c string containing the path to the txt file (without the suffix ".txt").
 */
void readOFS_txt(Ofsc& xFFT, string filename);

#endif // OFS_H
