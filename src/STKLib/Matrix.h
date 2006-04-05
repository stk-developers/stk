#ifndef STK_Matrix_h
#define STK_Matrix_h

#include "common.h"
#include "Error.h"
#include "BasicVector.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

#ifdef USE_BLAS
extern "C"{
  #include <cblas.h>
  #include <clapack.h>
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
  template<typename _ElemT> class Matrix;

  
  // we need to declare the friend << operator here
  template<typename _ElemT>
    std::ostream & operator << (std::ostream & out, Matrix<_ElemT> & m);
  
  
  
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

      // Constructors

      /// Empty constructor
      Matrix<_ElemT> ():
        mMRows(0), mMCols(0), mStride(0), mpData(NULL)
#ifdef STK_MEMALIGN_MANUAL
        , mpFreeData(NULL)        
#endif      
      {}

      /// Copy constructor
      Matrix<_ElemT> (const ThisType & t);

      /// Basic constructor
      Matrix<_ElemT> (const size_t r,
                      const size_t c);

      /// Destructor
      ~Matrix<_ElemT> ();


      /// Initializes matrix (if not done by constructor)
      void
      Init(const size_t r,
           const size_t c);
      
      /**
       * @brief Dealocates the matrix from memory and resets the dimensions to (0, 0)
       */
      void
      Destroy();
           
      const bool
      IsInitialized() const
      { return mpData != NULL; }            
            
      /// Returns number of rows in the matrix
      const size_t
      Rows() const
      {
        return mMRows;
      }

      /// Returns number of columns in the matrix
      const size_t
      Cols() const
      {
        return mMCols;
      }

      /// Returns number of columns in the matrix memory
      const size_t
      Stride() const
      {
        return mStride;
      }
      
      
      /**
       *  @brief Gives access to a specified matrix row without range check
       *  @return Pointer to the const array
       */
      const _ElemT* const
      cpData () const
      {
        return mpData;
      }
      
      
      /**
       *  @brief Gives access to a specified matrix row without range check
       *  @return Pointer to the non-const data array
       */
      _ElemT*      
      pData () const
      {
        return mpData;
      }

      
      /// Returns size of matrix in memory
      const size_t
      MSize() const
      {
        return mMRows * mStride * sizeof(_ElemT);
      }
      
      
      //########################################################################
      //########################################################################
      // Math stuff
      /**
       * @brief Performs diagonal scaling
       * @param pDiagVector Array representing matrix diagonal
       * @return Refference to this
       */
      ThisType&
      DiagScale(_ElemT* pDiagVector);
      
      ThisType&
      DiagScale(BasicVector<_ElemT>& rDiagVector);
      
      /**
       *  @brief Performs vector multiplication on a and b and and adds the
       *         result to this (elem by elem)
       */
      ThisType &
      AddMMMul(ThisType & a, ThisType & b);
      
      ThisType &
      AddMMTMul(ThisType & a, ThisType & b);

      ThisType &
      AddCMMtMul(_ElemT c, ThisType & a, ThisType & b);
            
      ThisType &
      AddMCMul(ThisType & a, _ElemT c);
      
      ThisType &
      RepMMSub(ThisType & a, ThisType & b);
      
      ThisType &
      RepMMMul(ThisType & a, ThisType & b);
      
      ThisType &
      RepMTMMul(ThisType & a, ThisType & b);
      
      ThisType &
      AddCVVtMul(_ElemT c, BasicVector<_ElemT>& rA, BasicVector<_ElemT>& rB);
      
      ThisType &
      AddCVVtMul(_ElemT c, _ElemT* pA, _ElemT* pB);

      ThisType &
      AddCVVt(_ElemT c, BasicVector<_ElemT>& rA, _ElemT* pB);
      
      //########################################################################
      //########################################################################
      
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
       *  @brief Performs matrix inversion
       */
      ThisType &
      Invert();
      
      
      

      ThisType &
      Clear();
      
      /**
       *  @brief Gives access to a specified matrix row without range check
       *  @return pointer to the first field of the row
       */
      _ElemT*      
      operator []  (size_t i) const
      {
        return mpData + (i * mStride);
      }
      
      /**
       *  @brief Gives access to a specified matrix row with range check
       *  @return pointer to the first field of the row
       */
      _ElemT*
      Row (const size_t r)
      {
        if (0 <= r && r < mMRows)
        {
          return this->operator[] (r);
        }
        else
        {
          throw std::out_of_range("Matrix row out of range");
        }
      }

      
      /**
       *  @brief Gives access to matrix elements (row, col)
       *  @return pointer to the desired field
       */
      _ElemT&
      operator () (const size_t r, const size_t c)
      { return *(mpData + r * mStride + c); }


      friend std::ostream & 
      operator << <> (std::ostream & out, ThisType & m);

            
      void PrintOut(char *file);

    
    protected:
      
      //@{
      /// these atributes store the real matrix size as it is stored in memory
      /// including memalignment
      size_t    mMRows;       ///< Number of rows
      size_t    mMCols;       ///< Number of columns
      size_t    mStride;      ///< true number of columns for the internal matrix.
                              ///< This number may differ from M_cols as memory
                              ///< alignment might be used
      /// data memory area
      _ElemT*   mpData;

#ifdef STK_MEMALIGN_MANUAL
      /// data to be freed (in case of manual memalignment use, see common.h)
      _ElemT*   mpFreeData;
#endif
    }; // class Matrix



  /** **************************************************************************
   ** **************************************************************************
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
      size_t  mTRowOff;         ///< First row of the window
      size_t  mTColOff;         ///< First column of the window
      //@}


    public:
      /// defines a type of this
      typedef WindowMatrix<_ElemT>    ThisType;


      /// Empty constructor
      WindowMatrix() : Matrix<_ElemT>() {};

      /// Copy constructor
      WindowMatrix(const ThisType & rT);

      /// Basic constructor
      WindowMatrix(const size_t r, const size_t c):                            
        Matrix<_ElemT>(r, c), // create the base class
        mOrigMRows    (Matrix<_ElemT>::mMRows),
        mOrigMCols    (Matrix<_ElemT>::mMCols),
        mOrigMRealCols(Matrix<_ElemT>::mStride),
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
        Matrix<_ElemT>::mMRows    =  mOrigMRows;
        Matrix<_ElemT>::mMCols    =  mOrigMCols;
        Matrix<_ElemT>::mStride   =  mOrigMRealCols;

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
    DiagScale(BasicVector<float>& rDiagVector);
  
  
  template<>
    Matrix<float> &
    Matrix<float>::
    DiagScale(float* pDiagVector);
    

      
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSoftmax();
  
  
  template<>
    Matrix<float> &
    Matrix<float>::
    AddMMMul(Matrix<float> & a, Matrix<float>& b);
  
  template<>
    Matrix<double> &
    Matrix<double>::
    AddMMMul(Matrix<double> & a, Matrix<double>& b);

    
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCMMtMul(float c, Matrix<float>& a, Matrix<float>& b);

  template<>
    Matrix<double> &
    Matrix<double>::
    AddCMMtMul(double c, Matrix<double>& a, Matrix<double>& b);

      
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCVVtMul(float c, BasicVector<float>& rA, BasicVector<float>& rB);

  template<>
    Matrix<double> &
    Matrix<double>::
    AddCVVtMul(double c, BasicVector<double>& rA, BasicVector<double>& rB);
    
} // namespace STK
    
//#ifndef STK_Matrix_h
#endif 
