#ifndef OFTS_H
#define OFTS_H

#include <fstream>
#include <vector>
#include <sstream>

#include "Ofsc.h"
#include "FTA.h"

class Oftsc
{
    private:

    int nv;            ///number of variables
    int order;         ///order

    int cnv;           ///number of variables of the coefficients
    int corder;        ///order of the coefficients

    Ofsc **coefs;       ///matrix of coefficients

public:

    //------------------
    //Create
    //------------------
    /**
     *  \brief Default constructor of the class Oftsc.
     */
    Oftsc();
    /**
     *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
     *  \param newNv: number of variables of the serie
     *  \param newOrder: order of the serie
     *  \param newCnv: number of variables of the coefficients
     *  \param newOrder: order of the coefficients
     */
    Oftsc(int newNv, int newOrder, int newCnv, int newCorder);
    /**
     *  \brief Constructor from a given Oftsc object (without any link).
     *  \param b:  a reference to the Oftsc object to copy in the new object
     */
    Oftsc(Oftsc const& b);


    //------------------
    //Getters
    //------------------
    /**
     *  \brief  Gets the order of the serie.
     *  \return the order of the serie as an \c int
     */
    int getOrder() const;

    /**
     *  \brief  Gets the number of variables of serie.
     *  \return the number of variables of the serie as an \c int
     */
    int getNV() const;

    /**
     *  \brief  Gets the adress of the coefficient at order \c ord and position \c pos
     *  \return a pointer to the desired coefficient
     */
    Ofsc* getCoef(int const& ord, int const& pos) const;

    //------------------
    //Delete
    //------------------
    /**
     *  \brief Default destructor of the class Oftsc. WARNING: potential memory leak here, through the terms of type Oftsh.
     */
    ~Oftsc();

    //------------------
    //Zeroing
    //------------------
    /**
     *  \brief  Sets all coefficients to zero.
     */
    void zero();

    //------------------
    //Evaluate
    //------------------
    void evaluate(cdouble X[], Ofsc& z, int const& m, int const& ofs_order);
    void evaluate(cdouble X[], Ofsc* z, int const& m);

    //------------------
    //Friendly streaming
    //------------------
    friend std::ostream& operator << (std::ostream& stream, Oftsc const& ofts);
};


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Reading & writing
//
//---------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------
// Text format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in txt format.
 **/
void readOFS_txt(Ofsc &xFFT, ifstream &readStream, int fftN);
/**
 * \brief Reads a given \c Oftsc  object, in txt format.
 **/
int  readOFTS_txt(Oftsc &x, string filename, int fftN);
/**
 * \brief Reads a given vector W of type \c Oftsc  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
void readVOFTS_txt(vector<Oftsc >  &W, string filename, int fftN);


//----------------------------------------------
// Binary format, read
//----------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in bin format.
 **/
void readOFS_bin(Ofsc  &xFFT, fstream &myfile);

/**
 * \brief Reads a given \c Oftsc  object, in bin format.
 **/
int readOFTS_bin(Oftsc &W, string filename, int fftN);

/**
 * \brief Reads a given vector W of type \c Oftsc  in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
void readVOFTS_bin(vector<Oftsc >  &W, string filename, int fftN);

//----------------------------------------------
// Text format, write
//----------------------------------------------
/**
 * \brief Writes a given vector W of type \c Oftsc  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
void  writeVOFTS_txt(vector<Oftsc > &W, string filename);

//----------------------------------------------
// Binary format, write
//----------------------------------------------
/**
 * \brief Writes a given \c Ofsc  object within a \c Oftsc, in bin format.
 **/
void  writeOFS_bin(Ofsc &xFFT, fstream &myfile);
/**
 * \brief Writes a given \c Oftsc  object, in bin format.
 **/
void  writeOFTS_bin(Oftsc &W, string filename);
/**
 * \brief Writes a given vector W of type \c Oftsc  in a binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
void  writeVOFTS_bin(vector<Oftsc > &W, string filename);



#endif // OFTS_H
