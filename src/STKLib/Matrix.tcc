#include <cstdlib>
#include <malloc.h>
#include <cblas.h>

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
      this->mMRows      = rows;
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
      // this->data = new char[t.M_size];
      // allocate the memory and set the right dimensions and parameters
      if (!posix_memalign(& data, 16, t.mMSize))
      {
        this->data  = static_cast<_ElemT *> (data);

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
  _ElemT *
  Matrix<_ElemT>::
  operator ()(const size_t r, const size_t c)
  {
    if (mStorageType == STORAGE_REGULAR)
    {
      return (mpData + r * mMRealCols + c);
    }
    else if (mStorageType == STORAGE_TRANSPOSED)
    {
      return (mpData + c * mMRealCols + r);
    }
    else
      return NULL;
  }



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
        // this->data = new char[t.M_size];
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
    operator << (std::ostream & rOut, const Matrix<_ElemT> & rM)
    {
      size_t    rcount;
      size_t    ccount;

      _ElemT  * data    = rM.data;

      // go through all rows
      for (rcount = 0; rcount < rM.rows(); rcount++)
      {
        // go through all columns
        for (ccount = 0; ccount < rM.cols(); ccount++)
        {
          rOut << static_cast<_ElemT> (*data) << " ";
          data ++;
        }

        data += rM.M_skip;

        rOut << std::endl;
      }
      return rOut;
    }

}// namespace STK
