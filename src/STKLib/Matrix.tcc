
/** @file Matrix.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */


#ifndef STK_Matrix_tcc
#define STK_Matrix_tcc

#pragma GCC system_header

#include "common.h"
#include <cstdlib>
#include <malloc.h>
//#include <cblas.h>
#include <cassert>
#include <cmath>
extern "C"{
  #include <cblas.h>
}


static union
{
  double d;
  struct
  {
    int j;
    int i;
  } n;
} qn_d2i;

#define QN_EXP_A (1048576/M_LN2)
#define QN_EXP_C 60801
//#define EXP2(y) (qn_d2i.n.j = (int) (QN_EXP_A*(y)) + (1072693248 - QN_EXP_C), qn_d2i.d)
#define FAST_EXP(y) (qn_d2i.n.i = (int) (QN_EXP_A*(y)) + (1072693248 - QN_EXP_C), qn_d2i.d)

/*
void sigmoid_vec(float *in, float *out, int size)
{
  while(size--) *out++ = 1.0/(1.0 + FAST_EXP(-*in++));
}
*/
/*
float i_max_double(float *a, int len) 
{
  int i;
  float max;
  max = a[0];
  for (i=1; i<len; i++) {
    if (a[i] > max) {
      max = a[i];
    }
  }
  return max;
}


/*
void softmax_vec(float *in, float *out, int size)
{
  int i;
  float maxa,sum;
  // first find the max
  maxa = i_max_double (in, size);
  // normalize, exp and get the sum
  sum = 0.0;
  for (i=0; i<size; i++) {
    out[i] = FAST_EXP(in[i] - maxa);
    sum += out[i];
  }
  // now normalize bu the sum
  for (i=0; i<size; i++) {
    out[i] /= sum;
  }
}
*/
namespace STK
{
//******************************************************************************
  // Basic constructor
  template<typename _ElemT>
  Matrix<_ElemT>::
  Matrix<_ElemT> (const size_t       r,
                  const size_t       c,
                  const StorageType  st)
  {
    Init(r, c, st);
  }


//******************************************************************************    
  template<typename _ElemT>
  void
  Matrix<_ElemT>::
  Init(const size_t r,
       const size_t c,
       const StorageType st)
  {
    // initialize some helping vars
    size_t  rows;
    size_t  cols;
    size_t  skip;
    size_t  real_cols;
    size_t  size;
    void *  data;

    // at this moment, nothing is clear
    mStorageType = STORAGE_UNDEFINED;

    // if we store transposed, we swap the rows and cols
    // we assume that the user does not specify STORAGE_UNDEFINED
    if (st == STORAGE_TRANSPOSED)
    { rows = c; cols = r; }
    else
    { rows = r; cols = c; }


    // compute the size of skip and real cols
    skip      = ((16 / sizeof(_ElemT)) - cols % (16 / sizeof(_ElemT))) % (16 / sizeof(_ElemT));
    real_cols = cols + skip;
    size      = rows * real_cols * sizeof(_ElemT);


    // allocate the memory and set the right dimensions and parameters
    if (!posix_memalign(& data, 16, size))
    {
      this->mpData        = static_cast<_ElemT *> (data);

      this->mStorageType = st;
      this->mainMRows      = rows;
      this->mMCols      = cols;
      this->mMRealCols  = real_cols;
      this->mMSkip      = skip;
      this->mMSize      = size;

      this->mTRows      = r;
      this->mTCols      = c;

      // set all bytes to 0
      memset(this->mpData, this->mMSize, 0);
    }
    else
    {
      throw std::bad_alloc();
    }
  } //
  

//******************************************************************************
  // Copy constructor
  template<typename _ElemT>
  Matrix<_ElemT>::
  Matrix<_ElemT> (const ThisType & t)
  {
    // we need to be sure that the storage type is defined
    if (t.st != STORAGE_UNDEFINED)
    {
      void * data;

      // first allocate the memory
      // this->mpData = new char[t.M_size];
      // allocate the memory and set the right dimensions and parameters
      if (!posix_memalign(& data, 16, t.mMSize))
      {
        this->mpData  = static_cast<_ElemT *> (data);

        // copy the memory block
        memcpy(this->mpData, t.mpData, t.mMSize);

        // set the parameters
        this->mStorageType = t.mStorageType;
        this->mMRows       = t.mMRows;
        this->mMCols       = t.mMCols;
        this->mMRealCols   = t.mMRealCols;
        this->mMSkip       = t.mMSkip;
        this->mMSize       = t.mMSize;
      }
      else
      {
        // throw bad_alloc exception if failure
        throw std::bad_alloc();
      }
    }
  }


//******************************************************************************
  template<typename _ElemT>
  Matrix<_ElemT> &
  Matrix<_ElemT>::
  operator += (const Matrix<_ElemT> & a)
  {
#ifdef CHECKSIZE
    // the sizes must match, else throw an exception
    if (a.mMSize != this->mMSize)
      throw std::logic_error("Sizes of matrices must be identical");
#endif

    // return a reference to this
    return *this;
  }


//******************************************************************************
  // The destructor
  template<typename _ElemT>
  Matrix<_ElemT>::
  ~Matrix<_ElemT> ()
  {
    // we need to free the data block if it was defined
    if (mStorageType != STORAGE_UNDEFINED)
    {
      free(mpData);
    }
  }

  
//******************************************************************************
  template<typename _ElemT>
  _ElemT &
  Matrix<_ElemT>::
  operator ()(const size_t r, const size_t c)
  {
    assert(mStorageType == STORAGE_REGULAR || mStorageType == STORAGE_TRANSPOSED);
    if (mStorageType == STORAGE_REGULAR)
    {
      return *(mpData + r * mMRealCols + c);
    }
    else
    {
      return *(mpData + c * mMRealCols + r);      
    }
    
  }

    
//******************************************************************************
  template<typename _ElemT>
  Matrix<_ElemT> &
  Matrix<_ElemT>::
  AddMatMult(ThisType & a, ThisType & b)
  { 
    if(!(a.Cols() == b.Rows() && this->Rows() == a.Rows() &&  this->Cols() == b.Cols()))
      STK::Error("Matrix multiply: bad matrix sizes (%d %d)*(%d %d) -> (%d %d)", a.Rows(), a.Cols(), b.Rows(), b.Cols(), this->Rows(), this->Cols());
    
    // :KLUDGE: Dirty for :o)        
    assert(a.mStorageType == STORAGE_REGULAR);
    for(int r = 0; r < this->Rows() ; r++){
      for(int c = 0; c < this->Cols(); c++){
        // for every out matrix cell
        for(int index = 0; index < a.Cols(); index++){
          (*this)(r, c) += a(r, index) * b(index, c);
          //std::cout << a(r, index);
        }
      }
    }
    
    STK::Warning("Could be slow, not well implemented");
    
    return *this;
  }; // AddMatMult(const ThisType & a, const ThisType & b)
  

//******************************************************************************
/*
  template<>
  inline
  Matrix<float> &
  Matrix<float>::
  AddMatMult(Matrix<float> & a, Matrix<float> & b)
  { 
    assert(a.Cols() == b.Rows());
    assert(this->Rows() == a.Rows());
    assert(this->Cols() == b.Cols());
    if(b.Storage() == STORAGE_TRANSPOSED){
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mMRealCols, b.mpData, b.mMRealCols, 1.0f, this->mpData, this->mMRealCols);
    }
    else{
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mMRealCols, b.mpData, b.mMRealCols, 1.0f, this->mpData, this->mMRealCols);
    }
    return *this;
  }; // AddMatMult(const ThisType & a, const ThisType & b)
*/

//******************************************************************************
  template<typename _ElemT>
  Matrix<_ElemT> &
  Matrix<_ElemT>::
  FastRowSigmoid()
  {
    STK::Error("Sigmoid only implemented for float");
    return *this;
  }
  

//******************************************************************************
  /*
  template<>
  inline
  Matrix<float> &
  Matrix<float>::
  FastRowSigmoid()
  {
    for(size_t row = 0; row < this->Rows(); row++)
    {
      sigmoid_vec(this->Row(row), this->Row(row), this->Cols());
    }
    return *this;
  }
  */
//******************************************************************************
  template<typename _ElemT>
  Matrix<_ElemT> &
  Matrix<_ElemT>::
  FastRowSoftmax()
  {
    STK::Error("Softmax only implemented for float");
    return *this;
  }
  
//******************************************************************************
/*
  template<>
  inline
  Matrix<float> &
  Matrix<float>::
  FastRowSoftmax()
  {
    for(size_t row = 0; row < this->Rows(); row++){
      softmax_vec(this->Row(row), this->Row(row), this->Cols());
    }
    return *this;
  }
*/
  template<typename _ElemT>
  Matrix<_ElemT> &
  Matrix<_ElemT>::
  Transpose()
  {
    if(this->mStorageType == STORAGE_REGULAR) 
      this->mStorageType = STORAGE_TRANSPOSED;
    else
      this->mStorageType = STORAGE_REGULAR;  
      
    int pom;
    
   /* pom = this->mMRows;
    this->mMRows = mMCols;
    this->mMCols = pom;*/

    pom = this->mTRows;
    this->mTRows = mTCols;
    this->mTCols = pom;  
    
    return *this;
  }
//******************************************************************************

//******************************************************************************
  // Copy constructor
  template<typename _ElemT>
    WindowMatrix<_ElemT>::
    WindowMatrix<_ElemT> (const ThisType & t) :
      Matrix<_ElemT> (t)
    {
      // we need to be sure that the storage type is defined
      if (t.st != STORAGE_UNDEFINED)
      {
        // the data should be copied by the parent constructor
        void * data;

        // first allocate the memory
        // this->mpData = new char[t.M_size];
        // allocate the memory and set the right dimensions and parameters
        if (!posix_memalign(& data, 16, t.mMSize))
        {
          Matrix<_ElemT>::mpData  = static_cast<_ElemT *> (data);

          // copy the memory block
          memcpy(this->mpData, t.mpData, t.mMSize);

          // set the parameters
          this->mStorageType = t.mStorageType;
          this->mMRows       = t.mMRows;
          this->mMCols       = t.mMCols;
          this->mMRealCols   = t.mMRealCols;
          this->mMSkip       = t.mMSkip;
          this->mMSize       = t.mMSize;
        }
        else
        {
          // throw bad_alloc exception if failure
          throw std::bad_alloc();
        }
      }
    }


//******************************************************************************
  template<typename _ElemT>
    void
    WindowMatrix<_ElemT>::
    SetSize(const size_t rowOff,
            const size_t rows)
    {
      if (Matrix<_ElemT>::mStorageType == STORAGE_REGULAR)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData = mpOrigData + mOrigMRealCols * rowOff;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = mOrigTCols;
        Matrix<_ElemT>::mMRows = rows;
        Matrix<_ElemT>::mMCols = mOrigMCols;
        Matrix<_ElemT>::mMRealCols = mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = mOrigMSkip;
      }
      else if (Matrix<_ElemT>::mStorageType == STORAGE_TRANSPOSED)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData = mpOrigData + rowOff;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = mOrigTCols;
        Matrix<_ElemT>::mMRows = mOrigMCols;
        Matrix<_ElemT>::mMCols = rows;
        Matrix<_ElemT>::mMRealCols = mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = mOrigMRealCols - rows;
      }
    }


//******************************************************************************
  template<typename _ElemT>
    void
    WindowMatrix<_ElemT>::
    SetSize(const size_t rowOff,
            const size_t rows,
            const size_t colOff,
            const size_t cols)
    {
      if (Matrix<_ElemT>::mStorageType == STORAGE_REGULAR)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData = mpOrigData + mOrigMRealCols * rowOff + cols;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = cols;
        Matrix<_ElemT>::mMRows = rows;
        Matrix<_ElemT>::mMCols = cols;
        Matrix<_ElemT>::mMRealCols = this->mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = this->mOrigMRealCols - cols;
      }
      else if (Matrix<_ElemT>::mStorageType == STORAGE_TRANSPOSED)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData = mpOrigData + rowOff;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = mOrigTCols;
        Matrix<_ElemT>::mMRows = mOrigMCols;
        Matrix<_ElemT>::mMCols = rows;
        Matrix<_ElemT>::mMRealCols = mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = mOrigMRealCols - rows;
      }

    }


  template<typename _ElemT>
    std::ostream &
    operator << (std::ostream & rOut, Matrix<_ElemT> & rM)
    {
      rOut << "Rows: " << rM.Rows() << "\n";
      rOut << "Cols: " << rM.Cols() << "\n";
      if(rM.Storage() == STORAGE_TRANSPOSED) rOut << "Transposed \n ";
      else rOut << "Regular \n ";
      
      size_t    rcount;
      size_t    ccount;

      _ElemT  * data    = rM.mpData;

      // go through all rows
      for (rcount = 0; rcount < rM.Rows(); rcount++)
      {
        // go through all columns
        for (ccount = 0; ccount < rM.Cols(); ccount++)
        {
          rOut << static_cast<_ElemT>(/*(*data)*/ rM(rcount, ccount)) << " ";
          /*data ++;*/
        }

        //data += rM.mMSkip;

        rOut << std::endl;
      }
      return rOut;
    }

}// namespace STK

// #define STK_Matrix_tcc
#endif 
