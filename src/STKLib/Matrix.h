#ifndef STK_Matrix_h
#define STK_Matrix_h

#include "common.h"
#include "Error.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

#ifdef USE_BlAS
extern "C"{
  #include <cblas.h>
}
#endif


#define CHECKSIZE

namespace STK
{

  //  class matrix_error : public std::logic_error {};/
  //  class matrix_sizes_error : public matrix_error {};
  
  /// defines a storage type
  typedef enum
  {
    STORAGE_UNDEFINED = 0,
    STORAGE_REGULAR,
    STORAGE_TRANSPOSED
  } StorageType;
  
  
  // declare the class so the header knows about it
  template<typename _ElemT> class BasicMatrix;
  template<typename _ElemT> class Matrix;

  
  // we need to declare the friend << operator here
  template<typename _ElemT>
    std::ostream & operator << (std::ostream & out, Matrix<_ElemT> & m);

  
  /** **************************************************************************
   ** **************************************************************************
   *  @brief Provides a matrix abstraction class
   *
   *  This class provides a way to work with matrices in STK.
   *  It encapsulates basic operations and memory optimizations.
   *
   */
  template<typename _ElemT>
    class BasicMatrix
    {
    public:
      /// defines a type of this
      typedef Matrix<_ElemT>    ThisType;

      BasicMatrix<_ElemT> (): mStorageType(STORAGE_UNDEFINED) {}
      
      /**
       *  @brief Gives access to the matrix memory area
       *  @return pointer to the first field
       */
      _ElemT*
      operator () () {return mpData;};

      /**
       *  @brief Gives access to the matrix memory area
       *  @return pointer to the first field
       */
      _ElemT*
      Data() {return mpData;};
      
      
      /// Returns the way the matrix is stored in memory
      const StorageType
      Storage() const
      {
        return mStorageType;
      }
      
      
    protected:
      /// keeps info about data layout in the memory
      StorageType mStorageType;
      
      /// data memory area
      _ElemT*   mpData;

#ifdef STK_MEMALIGN_MANUAL
      /// data to be freed (in case of manual memalignment use, see common.h)
      _ElemT*   mpFreeData;
#endif
      
    }; // class BasicMatrix
    
  
  
  
  /** **************************************************************************
   ** **************************************************************************
   *  @brief Provides a matrix class
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
      StorageType mStorageType;


      //@{
      /// these atributes store the real matrix size as it is stored in memory
      /// including memalignment
      size_t  mMRows;       ///< Number of rows
      size_t  mMCols;       ///< Number of columns
      size_t  mMRealCols;   ///< true number of columns for the internal matrix.
                            ///< This number may differ from M_cols as memory
                            ///< alignment might be used
      size_t  mMSize;       ///< Total size of data block in bytes
      size_t  mMSkip;       ///< Bytes to skip (memalign...)

      size_t  mTRows;       ///< Real number of rows (available to the user)
      size_t  mTCols;       ///< Real number of columns
      //@}

      /// data memory area
      _ElemT *  mpData;

#ifdef STK_MEMALIGN_MANUAL
      /// data to be freed (in case of manual memalignment use, see common.h)
      _ElemT *  mpFreeData;
#endif

    public:
      // Constructors

      /// Empty constructor
      Matrix<_ElemT> (): mStorageType(STORAGE_UNDEFINED) {}

      /// Copy constructor
      Matrix<_ElemT> (const ThisType & t);

      /// Basic constructor
      Matrix<_ElemT> (const size_t r,
                      const size_t c,
                      const StorageType st = STORAGE_REGULAR);

      /// Destructor
      ~Matrix<_ElemT> ();


      /// Initializes matrix (if not done by constructor)
      void
      Init(const size_t r,
           const size_t c,
           const StorageType st = STORAGE_REGULAR);
      
      /// Returns number of rows in the matrix
      size_t
      Rows() const
      {
        return mTRows;
      }

      /// Returns number of columns in the matrix
      size_t
      Cols() const
      {
        return mTCols;
      }

      /// Returns vector (matrix) length
      size_t
      Length() const
      {
        return mTRows > mTCols ? mTRows : mTCols;
      }
      
      /// Returns the way the matrix is stored in memory
      StorageType
      Storage() const
      {
        return mStorageType;
      }
      
      size_t
      Size() const
      {
        return mMSize;
      }


      /**
       *  @brief Performs vector multiplication on a and b and and adds the
       *         result to this (elem by elem)
       */
      ThisType &
      AddMMMul(ThisType & a, ThisType & b);
      
      ThisType &
      AddMCMul(ThisType & a, _ElemT c);
      
      ThisType &
      RepMMSub(ThisType & a, ThisType & b);
      
      ThisType &
      RepMMTMul(ThisType & a, ThisType & b);
      
      ThisType &
      RepMTMMul(ThisType & a, ThisType & b);
      
      /**
       *  @brief Performs fast sigmoid on row vectors
       *         result to this (elem by elem)
       */
      ThisType &
      FastRowSigmoid();
      
      /**
       *  @brief Performs fast softmax on row vectors
       *         result to this (elem by elem)
       */
      ThisType &
      FastRowSoftmax();
      
      /**
       *  @brief Performs matrix transposition
       */
      ThisType &
      Transpose();

      ThisType &
      Clear();      
    
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
      Row (const size_t r)
      {
        return mpData + (r * mMRealCols);
      }

      /**
       *  @brief Gives access to matrix elements (row, col)
       *  @return pointer to the desired field
       */
      _ElemT &
      operator () (const size_t r, const size_t c);


      friend std::ostream & 
      operator << <> (
        std::ostream & out, 
        ThisType & m);

        
      void PrintOut(char *file);

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
      /// The data atribute points now to the begining of the window
      _ElemT * mpOrigData;

      //@{
      /// these atributes store the real matrix size as it is stored in memory
      /// including memalignment
      size_t  mOrigMRows;       ///< Number of rows
      size_t  mOrigMCols;       ///< Number of columns
      size_t  mOrigMRealCols;   ///< true number of columns for the internal matrix.
                                ///< This number may differ from M_cols as memory
                                ///< alignment might be used
      size_t  mOrigMSize;       ///< Total size of data block in bytes
      size_t  mOrigMSkip;       ///< Bytes to skip (memalign...)

      size_t  mOrigTRows;       ///< Original number of window rows
      size_t  mOrigTCols;       ///< Original number of window columns

      size_t  mTRowOff;         ///< First row of the window
      size_t  mTColOff;         ///< First column of the window
      //@}


    public:
      /// defines a type of this
      typedef WindowMatrix<_ElemT>    ThisType;


      /// Empty constructor
      WindowMatrix() : Matrix<_ElemT>() {};

      /// Copy constructor
      WindowMatrix<_ElemT> (const ThisType & rT);

      /// Basic constructor
      WindowMatrix<_ElemT> (const size_t r,
                            const size_t c,
                            const StorageType st = STORAGE_REGULAR):                            
        Matrix<_ElemT>(r, c, st), // create the base class
        mOrigMRows    (Matrix<_ElemT>::mMRows),
        mOrigMCols    (Matrix<_ElemT>::mMCols),
        mOrigMRealCols(Matrix<_ElemT>::mMRealCols),
        mOrigMSize    (Matrix<_ElemT>::mMSize),
        mOrigMSkip    (Matrix<_ElemT>::mMSkip),
        mOrigTRows    (Matrix<_ElemT>::mTRows),   // copy the original values
        mOrigTCols    (Matrix<_ElemT>::mTCols),
        mTRowOff(0), mTColOff(0)          // set the offset
      {
        mpOrigData = Matrix<_ElemT>::mpData;
      }
      
      
      /// The destructor
      ~WindowMatrix<_ElemT>()
      {
        Matrix<_ElemT>::mpData = this->mpOrigData;
      }

      /// sets the size of the window (whole rows)
      void
      SetSize(const size_t ro,
              const size_t r);


      /// sets the size of the window
      void
      SetSize(const size_t ro,
              const size_t r,
              const size_t co,
              const size_t c);

      /**
       *  @brief Resets the window to the default (full) size
       */
      void
      Reset ()
      {
        Matrix<_ElemT>::mpData    =  mpOrigData;
        Matrix<_ElemT>::mTRows    =  mOrigTRows; // copy the original values
        Matrix<_ElemT>::mTCols    =  mOrigTCols;
        Matrix<_ElemT>::mMRows    =  mOrigMRows;
        Matrix<_ElemT>::mMCols    =  mOrigMCols;
        Matrix<_ElemT>::mMRealCols=  mOrigMRealCols;
        Matrix<_ElemT>::mMSize    =  mOrigMSize;
        Matrix<_ElemT>::mMSkip    =  mOrigMSkip;

        mTRowOff = 0;
        mTColOff = 0;
      }
    };

} // namespace STK



//*****************************************************************************
//*****************************************************************************
// we need to include the implementation
#include "Matrix.tcc"
//*****************************************************************************
//*****************************************************************************


/******************************************************************************
 ******************************************************************************
 * The following section contains specialized template definitions
 * whose implementation is in Matrix.cc
 */
 
namespace STK
{
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSigmoid();
    
  template<>
    Matrix<float> &
    Matrix<float>::
    AddMMMul(Matrix<float> & a, Matrix<float> & b);
    
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSoftmax();
} // namespace STK
    
//#ifndef STK_Matrix_h
#endif 
