#include <cstdlib>
#include <malloc.h>
#include <cblas.h>

namespace STK
{
//******************************************************************************
  // Basic constructor
  template<typename _ElemT>
  Matrix<_ElemT>::
  Matrix<_ElemT> (const size_t       _r,
                  const size_t       _c,
                  const storage_type _st)
  {
    // initialize some helping vars
    size_t  _rows;
    size_t  _cols;
    size_t  _skip;
    size_t  _realcols;
    size_t  _size;
    void * _data;

    // at this moment, nothing is clear
    mStorageType = STORAGE_UNDEFINED;

    // if we store transposed, we swap the rows and cols
    // we assume that the user does not specify STORAGE_UNDEFINED
    if (_st == STORAGE_TRANSPOSED)
    { _rows = _c; _cols = _r; }
    else
    { _rows = _r; _cols = _c; }


    // compute the size of skip and real cols
    _skip     = ((16 / sizeof(_ElemT)) - _cols % (16 / sizeof(_ElemT))) % (16 / sizeof(_ElemT));
    _realcols = _cols + _skip;
    _size     = _rows * _realcols * sizeof(_ElemT);


    // allocate the memory and set the right dimensions and parameters
    if (!posix_memalign(&_data, 16, _size))
    {
      this->data        = static_cast<_ElemT *> (_data);

      this->mStorageType = _st;
      this->mMRows      = _rows;
      this->mMCols      = _cols;
      this->mMRealCols  = _realcols;
      this->mMSkip      = _skip;
      this->mMSize      = _size;

      this->mTRows      = _r;
      this->mTCols      = _c;

      // set all bytes to 0
      memset(this->mpData, this->mMSize, 0);
    }
    else
    {
      throw std::bad_alloc();
    }
  }


//******************************************************************************
  // Copy constructor
  template<typename _ElemT>
  Matrix<_ElemT>::
  Matrix<_ElemT> (const __this_type & t)
  {
    // we need to be sure that the storage type is defined
    if (t.st != STORAGE_UNDEFINED)
    {
      void * _data;

      // first allocate the memory
      // this->mpData = new char[t.mMSize];
      // allocate the memory and set the right dimensions and parameters
      if (!posix_memalign(&_data, 16, t.mMSize))
      {
        this->mpData  = static_cast<_ElemT *> (_data);

        // copy the memory block
        memcpy(this->mpData, t.mpData, t.mMSize);

        // set the parameters
        this->mStorageType = t.mStorageType;
        this->mMRows      = t.mMRows;
        this->mMCols      = t.mMCols;
        this->mMRealCols  = t.mMRealCols;
        this->mMSkip      = t.mMSkip;
        this->mMSize      = t.mMSize;
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
    // we need to free the mpData block if it was defined
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
    if (this->mStorageType == STORAGE_REGULAR)
    {
      return (this->mpData + r * this->mMRealCols + c);
    }
    else if (this->mStorageType == STORAGE_TRANSPOSED)
    {
      return (this->mpData + c * this->mMRealCols + r);
    }
    else
      return NULL;
  }



//******************************************************************************
  // Copy constructor
  template<typename _ElemT>
    WindowMatrix<_ElemT>::
    WindowMatrix<_ElemT> (const __this_type & t) :
      Matrix<_ElemT> (t)
    {
      // we need to be sure that the storage type is defined
      if (t.st != STORAGE_UNDEFINED)
      {
        // the data should be copied by the parent constructor
        void * _data;

        // first allocate the memory
        // this->mpData = new char[t.mMSize];
        // allocate the memory and set the right dimensions and parameters
        if (!posix_memalign(&_data, 16, t.mMSize))
        {
          this->mpData  = static_cast<_ElemT *> (_data);

          // copy the memory block
          memcpy(this->mpData, t.mpData, t.mMSize);

          // set the parameters
          this->mStorageType = t.mStorageType;
          this->mMRows      = t.mMRows;
          this->mMCols      = t.mMCols;
          this->mMRealCols  = t.mMRealCols;
          this->mMSkip      = t.mMSkip;
          this->mMSize      = t.mMSize;
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
    setSize(const size_t rowOff,
            const size_t rows)
    {
      if (Matrix<_ElemT>::mStorageType == STORAGE_REGULAR)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData   = this->orig_data + this->mOrigMRealCols * rowOff;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = this->mOrigTCols;
        Matrix<_ElemT>::mMRows = rows;
        Matrix<_ElemT>::mMCols = this->mOrigMCols;
        Matrix<_ElemT>::mMRealCols = this->mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = this->mOrigMSkip;
      }
      else if (Matrix<_ElemT>::mStorageType == STORAGE_TRANSPOSED)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData   = this->orig_data + rowOff;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = this->mOrigTCols;
        Matrix<_ElemT>::mMRows = this->mOrigMCols;
        Matrix<_ElemT>::mMCols = rows;
        Matrix<_ElemT>::mMRealCols = this->mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = this->mOrigMRealCols - rows;
      }
    }


//******************************************************************************
  template<typename _ElemT>
    void
    WindowMatrix<_ElemT>::
    setSize(const size_t rowOff,
            const size_t rows,
            const size_t colOff,
            const size_t cols)
    {
      if (Matrix<_ElemT>::mStorageType == STORAGE_REGULAR)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData   = this->orig_data +
                                          this->orig_realCols * rowOff +
                                          cols;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = cols;
        Matrix<_ElemT>::mMRows = rows;
        Matrix<_ElemT>::mMCols = cols;
        Matrix<_ElemT>::mMRealCols = this->mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = this->orig_realCols - cols;
      }
      else if (Matrix<_ElemT>::mStorageType == STORAGE_TRANSPOSED)
      {
        // point to the begining of window
        Matrix<_ElemT>::mpData   = this->orig_data + rowOff;
        Matrix<_ElemT>::mTRows = rows;
        Matrix<_ElemT>::mTCols = this->mOrigTCols;
        Matrix<_ElemT>::mMRows = this->mOrigMCols;
        Matrix<_ElemT>::mMCols = rows;
        Matrix<_ElemT>::mMRealCols = this->mOrigMRealCols;
        Matrix<_ElemT>::mMSkip = this->orig_realCols - rows;
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

        data += rM.mMSkip;

        rOut << std::endl;
      }
      return rOut;
    }

}// namespace STK
