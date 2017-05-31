/**
 * \file matrix.tpp
 * \brief Some extension of the vector class of C++ for matrix manipulation
 * \author BLB
 * \date May 2015
 * \version 1.0
 */



//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default creator for matrix class. Never use.
 **/
template <typename T> matrix<T>::matrix()
{
    size1 = 0;//NV;
    size2 = 0;//REDUCED_NV;
    coef = *(new vector<T>());
}

/**
 *  \brief Default creator for matrix class, with size size1 x size2.
 **/
template <typename T> matrix<T>::matrix(const int size1_, const int size2_)
{
    size1 = size1_;
    size2 = size2_;
    coef = *(new vector<T>(size1*size2));
}


//---------------------------------------------------------------------------
//Copy
//---------------------------------------------------------------------------
/**
 *  \brief Create a matrix equal to the matrix b. Needs the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>::matrix(matrix const& b)
{
    size1 = b.size1;
    size2 = b.size2;
    coef  = *(new vector<T>(size1*size2));
    for(int i = 0 ; i< size1*size2; i++) coef[i].ccopy(b.coef[i]);
}

/**
 *  \brief Create a matrix equal to the matrix b. Needs the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>& matrix<T>::operator = (matrix<T> const& b)
{
    if(this != &b)
    {
        size1 = b.size1;
        size2 = b.size2;
        coef  = *(new vector<T>(size1*size2));
        for(int i = 0 ; i < size1*size2; i++) coef[i].ccopy(b.coef[i]);
    }
    return *this; //same object if returned
}

/**
 *  \brief Copy the matrix b in this. Needs the routine:
 *          - ccopy
 **/
template <typename T> matrix<T>& matrix<T>::ccopy(matrix<T> const& b)
{
    if(size1 != b.size1 || size2 != b.size2)
    {
        cout << "Erreur in ccopy for matrix: sizes do not match. Nothing is done." << endl;
        return *this;
    }
    else
    {

    for(int i = 0 ; i < size1*size2; i++) coef[i].ccopy(b.coef[i]);
    return *this; //same object if returned
    }
}

/**
 *  \brief Link the matrix b with this. Needs the routine:
 *          - lcopy
 **/
template <typename T> matrix<T>& matrix<T>::lcopy(matrix<T> const& b)
{
    size1 = b.size1;
    size2 = b.size2;
    for(int i = 0 ; i < size1*size2; i++) coef[i].lcopy(b.coef[i]);
    return *this;
}

//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Delete function, empty.
 *         Do we need to implement something here? Probably not, not pointer used.
 **/
template <typename T> matrix<T>::~matrix<T>()
{
}


//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief Gets the coefficient at the position (i,j)
 **/
template <typename T> T matrix<T>::getCoef(int i, int j) const
{
    if( i >= size1 || j >= size2)
    {
        cout << "Error in matrix<T>::getCoef: indices are out of scope. First coefficient is returned." << endl;
        return coef[0];
    }
    else return coef[i*size2 + j];
}

/**
 *  \brief Gets the address of the coefficient at the position (i,j)
 **/
template <typename T> T* matrix<T>::getCA(int i, int j) const
{
    if( i >= size1 || j >=size2)
    {
        cout << "Error in matrix<T>::getCoef: indices are out of scope. First coefficient is returned." << endl;
        return coef[0].getAddress();
    }
    else return coef[i*size2 + j].getAddress();
}

/**
 *  \brief Gets the size (either size1 or size2) of the matrix.
 **/
template <typename T> int matrix<T>::getSize(int num) const
{
    if(num == 1) return size1;
    else if(num == 2) return size2;
    else cout << "Error in matrix<T>::getSize: required number must be 1 or 2." << endl; return 0;

}

//---------------------------------------------------------------------------
//Setters
//---------------------------------------------------------------------------
/**
 *  \brief Sets the coefficient T at the position (i,j). Requires ccopy.
 **/
template <typename T> void matrix<T>::setCoef(T const & value, int i, int j)
{
    this->getCA(i,j)->ccopy(value);
}

/**
 *  \brief Sets the subcoefficient value at the order zero of the coefficient (i,j). Requires setCoef.
 **/
template <typename T> template <typename U> void matrix<T>::setCoef(U const & value, int i, int j)
{
    this->getCA(i,j)->setCoef(value, (int const) 0);
}

/**
 *  \brief Adds the coefficient T at the position (i,j). Requires smult.
 **/
template <typename T> void matrix<T>::addCoef(T const & value, int i, int j)
{
   this->getCoef(i,j).smult(value, 1.0);
}

/**
 *  \brief Zeroing of the matrix. Requires zero().
 **/
template <typename T> void matrix<T>::zero()
{
    for(int i = 0 ; i < size1*size2; i++) coef[i].zero();
}

//---------------------------------------------------------------------------
//Functions used with T = Ofs<U>
//---------------------------------------------------------------------------
/**
 *  \brief Matrix-vector product: vOut += a x vIn. Used with T = Ofsc. Requires sprod.
 **/
inline void smvprod_ofs(matrix<Ofsc> const& a, vector<Ofsc> const& vIn, vector<Ofsc>& vOut)
{
    if((unsigned int)  a.getSize(2) != vIn.size() || (unsigned int) a.getSize(1) != vOut.size() )
    {
        cout << "Error in smvprod_ofs (matrix.tpp): lengths do not match. Nothing is done." << endl;
    }
    else
    {
        for(int i =0; i < a.getSize(1) ; i++)
        {
            for(int j =0; j < a.getSize(2); j++)
            {
                ( (Ofsc*)&vOut[i])->ofs_sprod(a.getCoef(i,j), vIn[j]);
            }
        }
    }
}
