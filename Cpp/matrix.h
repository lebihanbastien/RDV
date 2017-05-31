#ifndef MATRIX_H
#define MATRIX_H


#include "Oftsc.h"
#include "parameters.h"

/**
 * \file matrix.h
 * \brief Some extension of the vector class of C++ for matrix manipulation
 * \author BLB
 * \date May 2015
 * \version 1.0
 */


using namespace std;

//--------------------------------------------------------------------------
// matrix<T> class
//--------------------------------------------------------------------------
template<typename T>
class matrix;



template <typename T>
class matrix
{
private:
    int size1;      //rows
    int size2;      //columns
    vector<T> coef; //coefficients

public:
    //--------------------------------------------
    //Create
    //--------------------------------------------
    matrix<T>();
    matrix<T>(const int size1_, const int size2_);

    //--------------------------------------------
    //Copy
    //--------------------------------------------
    matrix<T>(matrix const& ofs_);
    matrix<T>& ccopy(matrix<T> const& b);
    matrix<T>& lcopy(matrix<T> const& b);

    //--------------------------------------------
    //Delete
    //--------------------------------------------
    ~matrix<T>();

    //--------------------------------------------
    //Setters
    //--------------------------------------------
    void setCoef(T const & value, int i, int j);
    void addCoef(T const & value, int i, int j);
    template <typename U> void setCoef(U const & value, int i, int j);
    void setCoef(cdouble const & value, int i, int j);

    //--------------------------------------------
    //Getters
    //--------------------------------------------
    T getCoef(int i, int j) const;
    T* getCA(int i, int j) const;
    int getSize(int num) const;

    //--------------------------------------------
    //Operators
    //--------------------------------------------
    matrix<T>& operator  = (matrix<T> const& b);

    //--------------------------------------------
    //Operations
    //--------------------------------------------
    void der(T const &a, int ni, int i, int j);
    void der(T const &a, int ni, int i, int j, int k);
    void dot(T const &a, double n, int i, int j);
    void dot(matrix<T> const &a, double n);
    void zero();
    void tfts_der(T const &a, int ni, int i, int j, int k);


    T operator [](int i) const    {return coef[i];}
    T & operator [](int i) {return coef[i];}
    T operator ()(int i, int j) const    {return coef[size2*i+j];}
    T & operator ()(int i, int j) {return coef[size2*i+j];}

};


//---------------------------------------------------------------------------
//Functions used with T = Ofs<U>
//---------------------------------------------------------------------------
void smvprod_ofs(matrix<Ofsc> const& a, vector<Ofsc> const& vIn, vector<Ofsc>& vOut);

//---------------------------------------------------------------------------
//Include the implementation .tpp
//---------------------------------------------------------------------------
#include "matrix.tpp"

#endif // MATRIX_H
