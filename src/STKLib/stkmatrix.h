#ifndef __STKMATRIX_H
#define __STKMATRIX_H

#include <cstddef>
//#include <cstdlib>
#include <stdlib.h>
#include <stdexcept>
#include <cblas.h>

#include <iostream>

#define CHECKSIZE


//  class matrix_error : public std::logic_error {};/
//  class matrix_sizes_error : public matrix_error {};

/// defines a storage type
typedef enum
{
  STORAGE_UNDEFINED = 0,
  STORAGE_REGULAR,
  STORAGE_TRANSPOSED
} storage_type;


  // declare the class so the header knows about it
  template<typename _ElemT>
    class Matrix;

   

  // we need to declare the friend << operator here
  template<typename _ElemT>
    std::ostream & operator << (std::ostream & out, const Matrix<_ElemT> & m);


  /**
   *  @brief Provides a matrix abstraction class
   *
   *  This class provides a way to work with matrices in STK.
   *  It encapsulates basic operations and memory optimizations.
   *
   */
  template<typename _ElemT>
    class Matrix
    {
    public:
      /// defines a type of this
      typedef Matrix<_ElemT>    ThisType;


    protected:
      /// keeps info about data layout in the memory
      storage_type mStorageType;


      //@{
      /// these atributes store the real matrix size as it is stored in memory
      /// including memalignment
      size_t  mMRows;       ///< Number of rows
      size_t  mMCols;       ///< Number of columns
      size_t  mMRealCols;   ///< true number of columns for the internal matrix.
                            ///< This number may differ from mMCols as memory
                            ///< alignment might be used
      size_t  mMSize;       ///< Total size of data block in bytes
      size_t  mMSkip;       ///< Bytes to skip (memalign...)

      size_t  mTRows;       ///< Real number of rows (available to the user)
      size_t  mTCols;       ///< Real number of columns
      //@}

      /// data memory area
      _ElemT * mpData;


    public:
      // Constructors

      /// Empty constructor
      Matrix<_ElemT> (): mStorageType(STORAGE_UNDEFINED) {}

      /// Copy constructor
      Matrix<_ElemT> (const ThisType & t);

      /// Basic constructor
      Matrix<_ElemT> (const size_t r,
                               const size_t c,
                               const storage_type st = STORAGE_REGULAR);

      /// Destructor
      ~Matrix<_ElemT> ();


      /// Returns number of rows in the matrix
      size_t
      rows() const
      {
        return mTRows;
      }

      /// Returns number of columns in the matrix
      size_t
      cols() const
      {
        return mTCols;
      }


      /// Returns the way the matrix is stored in memory
      storage_type
      storage() const
      {
        return mStorageType;
      }


      /**
       *  @brief Adds values of matrix a to this element by element
       *  @param a Matrix to be added to this
       *  @return Reference to this
       */
      ThisType &
      operator += (const ThisType & a);

      /**
       *  @brief Performs vector multiplication on a and b and and adds the
       *         result to this (elem by elem)
       */
      ThisType &
      AddVectorMult(const ThisType & a, const ThisType & b);


      /**
       *  @brief Gives access to the matrix memory area
       *  @return pointer to the first field
       */
      _ElemT *
      operator () () {return mpData;};

      /**
       *  @brief Gives access to a specified matrix row
       *  @return pointer to the first field of the row
       */
      _ElemT *
      row (const size_t r)
      {
        return mpData + (r * mMRealCols);
      }

      /**
       *  @brief Gives access to matrix elements (row, col)
       *  @return pointer to the desired field
       */
      _ElemT *
      operator () (const size_t r, const size_t c);


      friend std::ostream & 
      operator << <> (std::ostream & out, const ThisType &m);

      void printOut();

    }; // class Matrix



  /**
   *  @brief Provides a window matrix abstraction class
   *
   *  This class provides a way to work with matrix cutouts in STK.
   *  It encapsulates basic operations and memory optimizations.
   *
   */
  template<typename _ElemT>
    class WindowMatrix : public Matrix<_ElemT>
    {
    protected:
      /// points to the original begining of the data array
      /// The mpData atribute points now to the begining of the window
      _ElemT * orig_data;

      ///@{
      /// these atributes store the real matrix size as it is stored in memory
      /// including memalignment
      size_t  mOrigMRows;       ///< Number of rows
      size_t  mOrigMCols;       ///< Number of columns
      size_t  mOrigMRealCols;   ///< true number of columns for the internal matrix.
                                ///< This number may differ from mMCols as memory
                                ///< alignment might be used
      size_t  mOrigMSize;       ///< Total size of mpData block in bytes
      size_t  mOrigMSkip;       ///< Bytes to skip (memalign...)

      size_t  mOrigTRows;       ///< Original number of window rows
      size_t  mOrigTCols;       ///< Original number of window columns

      size_t  mTRowOffset;      ///< First row of the window
      size_t  mTColOffset;      ///< First column of the window
      ///@}


    public:
      /// defines a type of this
      typedef WindowMatrix<_ElemT>    ThisType;


      /// Empty constructor
      WindowMatrix() : Matrix<_ElemT>() {};

      /// Copy constructor
      WindowMatrix<_ElemT> (const ThisType & t);

      /// Basic constructor
      WindowMatrix<_ElemT> (const size_t r,
                                  const size_t c,
                                  const storage_type st = STORAGE_REGULAR):
        Matrix<_ElemT>(r, c, st), // create the base class
        mTRowOffset(0), mTColOffset(0),          // set the offset
        mOrigTRows    (Matrix<_ElemT>::T_rows),   // copy the original values
        mOrigTCols    (Matrix<_ElemT>::T_cols),
        mOrigMRows    (Matrix<_ElemT>::M_rows),
        mOrigMCols    (Matrix<_ElemT>::M_cols),
        mOrigMRealCols(Matrix<_ElemT>::M_realCols),
        mOrigMSize    (Matrix<_ElemT>::M_size),
        mOrigMSkip    (Matrix<_ElemT>::M_skip)
      {
        orig_data = Matrix<_ElemT>::mpData;
      }

      /// The destructor
      ~WindowMatrix<_ElemT>()
      {
        Matrix<_ElemT>::mpData = this->orig_data;
      }

      /// sets the size of the window (whole rows)
      void
      setSize(const size_t _ro,
              const size_t _r);


      /// sets the size of the window
      void
      setSize(const size_t _ro,
              const size_t _r,
              const size_t _co,
              const size_t _c);

      /**
       *  @brief Resets the window to the default (full) size
       */
      void
      reset ()
      {
        Matrix<_ElemT>::mpData    =  orig_data;
        Matrix<_ElemT>::mTRows    =  mOrigTRows; // copy the original values
        Matrix<_ElemT>::mTCols    =  mOrigTCols;
        Matrix<_ElemT>::mMRows    =  mOrigMRows;
        Matrix<_ElemT>::mMCols    =  mOrigMCols;
        Matrix<_ElemT>::mMRealCols=  mOrigMRealCols;
        Matrix<_ElemT>::mMSize    =  mOrigMSize;
        Matrix<_ElemT>::mMSkip    =  mOrigMSkip;

        mTRowOffset = 0;
        mTColOffset = 0;
      }
    };

// we need to include the implementation
#include "stkmatrix.tcc"

#endif //#ifndef __STKMATRIX_H
