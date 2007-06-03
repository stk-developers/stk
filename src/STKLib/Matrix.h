#ifndef STK_Matrix_h
#define STK_Matrix_h

#include "common.h"
#include "Error.h"
#include "BasicVector.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

#ifdef HAVE_ATLAS
extern "C"{
  #include <cblas.h>
  #include <clapack.h>
}
#endif

//#define TRACE_MATRIX_OPERATIONS
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
  template<typename _ElemT> class MatrixRange;
  
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
      virtual
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
           
      /**
       * @brief Returns @c true if matrix is initialized
       */
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
      DiagScale(const _ElemT* pDiagVector);
      
      ThisType&
      DiagScale(const BasicVector<_ElemT>& rDiagVector);
      
      /**
       *  @brief Performs vector multiplication on a and b and and adds the
       *         result to this (elem by elem)
       */
      ThisType &
      AddCMtMMul(const _ElemT& c, const ThisType& a, const ThisType& b);
      
      ThisType &
      AddMMMul(const ThisType& a, const ThisType& b);
      
      ThisType &
      AddMMtMul(const ThisType& a, const ThisType& b);

      ThisType &
      AddCMMtMul(const _ElemT& c, const ThisType& a, const ThisType& b);
            
      ThisType &
      AddCMMul(const _ElemT& c, const ThisType& a);
      
      ThisType &
      RepMMSub(const ThisType& a, const ThisType& b);
      
      ThisType &
      RepMMMul(const ThisType& a, const ThisType& b);
      
      ThisType &
      RepMtMMul(const ThisType& a, const ThisType& b);
      
      ThisType &
      AddCVVtMul(const _ElemT& c, const BasicVector<_ElemT>& rA, 
                 const BasicVector<_ElemT>& rB);
      
      ThisType &
      AddCVVtMul(const _ElemT& c, const _ElemT* pA, const _ElemT* pB);

      ThisType &
      AddCVVt(const _ElemT& c, const BasicVector<_ElemT>& rA, const _ElemT* pB);
      
      ThisType &
      DivC(const _ElemT& c);
      
      
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
      inline _ElemT*      
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
      void ReadIn(char *file);

      /**
       * @brief Returns a matrix sub-range
       * @param ro Row offset
       * @param r  Rows in range
       * @param co Column offset
       * @param c  Coluns in range
       * See @c MatrixRange class for details
       */
      MatrixRange<_ElemT>
      Range(const size_t    ro, const size_t    r, 
            const size_t    co, const size_t    c)
      { return MatrixRange<_ElemT>(*this, ro, r, co, c); }

    
    protected:
      
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

    
    
    
  /** **************************************************************************
   ** **************************************************************************
   *  @brief Sub-matrix representation
   *
   *  This class provides a way to work with matrix cutouts in STK.
   *  
   *
   */
  template<typename _ElemT>
    class MatrixRange : public Matrix<_ElemT>
    {
    public:
      /// Constructor
      MatrixRange(const Matrix<_ElemT>& rT, 
                  const size_t    ro,
                  const size_t    r,
                  const size_t    co,
                  const size_t    c);

      /// The destructor
      virtual
      ~MatrixRange<_ElemT>()
      { 
#ifndef STK_MEMALIGN_MANUAL
        Matrix<_ElemT>::mpData = NULL;
#else
        Matrix<_ElemT>::mpFreeData = NULL;
#endif
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
    Matrix<double> &
    Matrix<double>::
    FastRowSigmoid();
    
  
  template<>
    Matrix<float> &
    Matrix<float>::
    DiagScale(const BasicVector<float>& rDiagVector);
  
  
  template<>
    Matrix<float> &
    Matrix<float>::
    DiagScale(const float* pDiagVector);
    

      
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSoftmax();
  
  template<>
    Matrix<double> &
    Matrix<double>::
    FastRowSoftmax();
  
  
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCMtMMul(const float& c, const Matrix<float>& a, 
                               const Matrix<float>& b);
  
  template<>
    Matrix<double> &
    Matrix<double>::
    AddCMtMMul(const double& c, const Matrix<double>& a, 
                                const Matrix<double>& b);

      
  template<>
    Matrix<float> &
    Matrix<float>::
    AddMMMul(const Matrix<float>& a, const Matrix<float>& b);
  
  template<>
    Matrix<double> &
    Matrix<double>::
    AddMMMul(const Matrix<double> & a, const Matrix<double>& b);

    
  template<>
    Matrix<float> &
    Matrix<float>::
    RepMMMul(const Matrix<float> & a, const Matrix<float>& b);
    
    
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCMMtMul(const float& c, const Matrix<float>& a, 
                               const Matrix<float>& b);

  template<>
    Matrix<double> &
    Matrix<double>::
    AddCMMtMul(const double& c, const Matrix<double>& a, 
                                const Matrix<double>& b);

      
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCVVtMul(const float& c, const BasicVector<float>& rA, 
                               const BasicVector<float>& rB);

  template<>
    Matrix<double> &
    Matrix<double>::
    AddCVVtMul(const double& c, const BasicVector<double>& rA, 
                                const BasicVector<double>& rB);
    
} // namespace STK
    
//#ifndef STK_Matrix_h
#endif 
