#include "Oftsc.h"

int OFTS_ORDER;
int REDUCED_NV;

/**
 * \file ofts.cpp
 * \brief Fourier-Taylor series template class (src)
 * \author BLB
 * \date 2016
 * \version 1.0
 */

//---------------------------------------------------------------------------
//Create
//---------------------------------------------------------------------------
/**
 *  \brief Default constructor of the class Oftsc.
 */
Oftsc::Oftsc()
{
   int m;
   //Update the different orders and number of variables
   nv     = REDUCED_NV;
   order  = OFTS_ORDER;
   cnv    = OFS_NV;
   corder = OFS_ORDER;

   //Allocate the matrix of Ofsc coefficients
   coefs = (Ofsc**) new Ofsc*[binomial(nv+order,nv)]; //coefs = (Ofsc**) calloc(binomial(nv+order,nv), sizeof(Ofsc*));
   for (int k=0; k<=order; k++)
   {
      //Allocate each columns of coefs
      m = FTA::nmon(nv, k);  //number of monomials of order k
      coefs[k] = (Ofsc*) new Ofsc[m];  //coefs[k]=(Ofsc*) calloc(m, sizeof(Ofsc));
      //Allocate each coefficient
      //for(int p = 0; p < m; p++) coefs[k][p] = Ofsc(corder);
   }
}

/**
 *  \brief Constructor with given order/number of variables both for the Fourier-Taylor series and the coefficients.
 *  \param newNv: number of variables of the serie
 *  \param newOrder: order of the serie
 *  \param newCnv: number of variables of the coefficients
 *  \param newOrder: order of the coefficients
 */
Oftsc::Oftsc(int newNv, int newOrder, int newCnv, int newCorder)
{
   int m;
   //Update the different orders and number of variables
   nv     = newNv;
   order  = newOrder;
   cnv    = newCnv;
   corder = newCorder;

   //Allocate the matrix of Ofsc coefficients
   coefs = (Ofsc**) new Ofsc*[binomial(nv+order,nv)]; //coefs = (Ofsc**) calloc(binomial(nv+order,nv), sizeof(Ofsc*));
   for (int k=0; k<=order; k++)
   {
      //Allocate each columns of coefs
      m = FTA::nmon(nv, k);  //number of monomials of order k
      coefs[k] = (Ofsc*) new Ofsc[m];  //coefs[k]=(Ofsc*) calloc(m, sizeof(Ofsc));
      //Allocate each coefficient
      //for(int p = 0; p < m; p++) coefs[k][p] = Ofsc(corder);
   }

}


/**
 *  \brief Constructor from a given Oftsc object (without any link).
 */
Oftsc::Oftsc(Oftsc const& b)
{
    int m;
    //----------------------------------
    //Same nv/order
    //----------------------------------
    nv = b.nv;
    order = b.order;
    cnv = b.cnv;
    corder = b.corder;

    //----------------------------------
    //Copy of all the coefficients at every order in new array
    //----------------------------------
    //Allocate the matrix of Ofsc coefficients
    coefs = (Ofsc**) new Ofsc*[binomial(nv+order,nv)];
    for (int k=0; k<=order; k++)
    {
        //Allocate each columns of coefs
        m = FTA::nmon(nv, k);  //number of monomials of order k
        coefs[k] = (Ofsc*) new Ofsc[m];
        //Allocate each coefficient
        //for(int p = 0; p < m; p++) coefs[k][p] = Ofsc(corder);
    }

    //Copy the coefficients of b
    for(int k=0; k<= order; k++)
    {
        for (int p=0; p< FTA::nmon(nv, k); p++)  coefs[k][p] = b.coefs[k][p];
    }
}

//---------------------------------------------------------------------------
//Delete
//---------------------------------------------------------------------------
/**
 *  \brief Default destructor of the class Oftsc. WARNING: potential memory leak here, through the terms of type Oftsh.
 *
 * Certainly a problem here: the line thats delete the coefficient is commented, because it leads to an error when programs end.
 * May lead to memory leak if the objects are created "on the fly", which may be the case in some inner routines like smprod_t.
 */
Oftsc::~Oftsc()
{
   for (int k=0; k<=order; k++)
   {
      if(coefs[k] != NULL) delete[] coefs[k];
   }
   delete[] coefs;
}

//---------------------------------------------------------------------------
//Getters
//---------------------------------------------------------------------------
/**
 *  \brief  Gets the order of the serie.
 */
int Oftsc::getOrder() const
{
    return order;
}

/**
 *  \brief  Gets the number of variables of serie.
 */
int Oftsc::getNV() const
{
    return nv;
}

/**
 *  \brief  Gets the adress of the coefficient at order \c ord and position \c pos
 */
 Ofsc* Oftsc::getCoef(int const& ord, int const& pos) const
{
    if(ord > order || pos >= FTA::nmon(nv, ord))
    {
        cout << "Error in getCoef: out of range. First term is returned" << endl;
        cout << "Requested order: " << ord << ", Maximum allowed: " <<  order << endl;
        cout << "Requested pos:   " << pos << ", Maximum allowed: " <<  FTA::nmon(nv, ord) << endl;
        return coefs[0];
    }
    else return &coefs[ord][pos];
}



//---------------------------------------------------------------------------
//Zeroing
//---------------------------------------------------------------------------
/**
 *  \brief  Sets all coefficients to zero.
 */
void Oftsc::zero()
{
    for(int nrc=0; nrc<= order; nrc++)
        for (int i=0; i< FTA::nmon(nv, nrc); i++)
            coefs[nrc][i].zero();
}


//Evaluate up to order m
void Oftsc::evaluate(cdouble X[], Ofsc& z, int const& m, int const& ofs_order)
{
    //Parameters
    int *kv = (int*) calloc(nv, sizeof(int));
    cdouble aux, bux;

    //Zeroing the result
    z.zero();

    //For each order in [[k, 0]]
    for(int k = m; k >= 0 ; k--)
    {
        //kv = (k 0 0 0 ...)
        kv[0] = k;
        for(int i=1; i<nv; i++) kv[i] = 0;

        //Loop on all monomials of degree k
        for (int i=0; i< FTA::nmon(nv, k); i++)
        {
            //Evaluate one coefficient
            if(!coefs[k][i].isnull(1))
            {
                //z += X[ii]^kv[ii]*coef(i)
                bux = 1.0+0.0*I;
                for(int ii = 0; ii < nv; ii++)
                {
                    aux = 1.0+0.0*I;
                    if(kv[ii] != 0.0) for(int j = 1; j <= kv[ii]; j++) aux*= X[ii];
                    bux *= aux;
                }
                z.ofs_smult(coefs[k][i], bux, ofs_order);
            }

            //update the exponents
            if(i< FTA::nmon(nv, k)-1)  FTA::prxkt(kv, nv);
        }

    }

    free(kv);
}

//---------------------------------------------------------------------------
//Evaluate
//---------------------------------------------------------------------------

std::ostream& operator << (std::ostream& stream, Oftsc const& ofts)
{
    int k[ofts.nv];

    for(int nrc=0; nrc<= ofts.order; nrc++)
    {
        //Update k
        k[0] = nrc;
        for(int i=1; i<ofts.nv; i++) k[i] = 0;
        //Current homogeneous polynomial
        for (int i=0; i< FTA::nmon(ofts.nv, nrc); i++)
        {
            //Print the current exponents
            for(int j=0; j<ofts.nv; j++) stream <<   setiosflags(ios::right) <<  k[j] << " ";
            stream << endl;
            //Print the Fourier series
            stream << std::showpos << setiosflags(ios::scientific)  << setprecision(15)  <<  ofts.coefs[nrc][i] << std::noshowpos << endl;
            //Update the exponents
            if(i< FTA::nmon(ofts.nv, nrc)-1)  FTA::prxkt(k, ofts.nv);
        }
    }
    return stream;
}

//---------------------------------------------------------------------------
// Text format, read
//---------------------------------------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in txt format.
 **/
void readOFS_txt(Ofsc &xFFT, ifstream &readStream, int fftN)
{
    //Init
    double ct, cr, ci;
    //Reading
    for(int i = -fftN; i<=fftN; i++)
    {
        readStream >> ct;  //current order
        readStream >> cr;  //real part
        readStream >> ci;  //imag part
        xFFT.setCoef(cr+I*ci, i);
    }
}

/**
 * \brief Reads a given \c Oftsc  object, in txt format.
 **/
int readOFTS_txt(Oftsc &x, string filename, int fftN)
{
    //Init
    ifstream readStream;
    string ct;
    //Reading
    readStream.open((filename).c_str(), ios::in);

    //Check that the opening went well
    if (!readStream.is_open())
    {
        cout << "readOFTS_txt. Cannot open file " << filename << endl;
        cout << "Check the text data exist." << endl;
        return -1;
    }
    else
    {
        for(int k = 0 ; k <= x.getOrder() ; k++)
        {
            for(int p = 0; p < FTA::nmon(x.getNV(), k); p++)
            {
                //Current kv
                getline(readStream, ct);
                //Reading the coefficient
                readOFS_txt(*x.getCoef(k,p), readStream, fftN);
                getline(readStream, ct);
                getline(readStream, ct);
            }
        }
        readStream.close();
        return 0;
    }
}

/**
 * \brief Reads a given vector W of type \c Oftsc  in a txt files of the form "filename+i.txt", with i = 0, length(W)-1.
 **/
void readVOFTS_txt(vector<Oftsc> &W, string filename, int fftN)
{
    string ss1;
    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
        readOFTS_txt(W[i], (filename+"["+ss1+"].txt"), fftN);
    }
}


//---------------------------------------------------------------------------
// Binary format, read
//---------------------------------------------------------------------------
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in bin format.
 **/
void readOFS_bin(Ofsc  &xFFT, fstream &myfile)
{
    int fftN = xFFT.getOrder();
    double cr, ci;
    //Writing
    for(int i = -fftN; i<=fftN; i++)
    {
        //Real part
        myfile.read((char*)&cr, sizeof(double));
        //Imag part
        myfile.read((char*)&ci, sizeof(double));
        //Put in current position
        xFFT.setCoef(cr+I*ci, i);
    }
}

/**
 * \brief Reads a given \c Oftsc  object, in bin format.
 **/
int readOFTS_bin(Oftsc &W, string filename, int fftN)
{
    //Init
    fstream myfile;

    //Open the stream
    myfile.open((filename).c_str(), ios::binary | ios::in);

    //Check that the opening went well
    if (!myfile.is_open())
    {
        cout << "readOFTS_bin. Cannot open file " << filename << endl;
        cout << "readOFTS_bin. Check the binary data exist." << endl;
        return -1;
    }
    else
    {
        //Loop on order
        for(int nrc=0; nrc<= W.getOrder(); nrc++)
        {
            //Current homogeneous polynomial
            for (int i=0; i< FTA::nmon(W.getNV(), nrc); i++)
            {
                //Read each Ofsc coefficient
                readOFS_bin(*W.getCoef(nrc,i), myfile);
            }
        }
        myfile.close();
        return 0;
    }
}

/**
 * \brief Reads a given vector W of type \c Oftsc  in binary files of the form "filename+i.bin", with i = 0, length(W)-1.
 **/
void readVOFTS_bin(vector<Oftsc >  &W, string filename, int fftN)
{
    string ss1;
    int status, global_status;
    global_status = 0;

    //Loop on all coefficients
    for(unsigned int i = 0; i < W.size(); i++)
    {
        ss1 = static_cast<ostringstream*>( &(ostringstream() << i) )->str();

        //Read binary format
        status = readOFTS_bin(W[i], (filename+"["+ss1+"].bin"), fftN);

        //Try txt format if failure
        if(status != 0)
        {
            cout << "readVOFTS_bin. Last reading went wrong. Trying to find data in txt format..." << endl;
            status = readOFTS_txt(W[i], (filename+"["+ss1+"].txt"), fftN);
            if(status != 0)
            {
                cout << "readVOFTS_bin. Txt format also went wrong. Check data manually." << endl;
                global_status--;
            }
            else
            {
                cout << "readVOFTS_bin. Success with txt format." << endl;
                global_status++;
            }
        }
    }
}
