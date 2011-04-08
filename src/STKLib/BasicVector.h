//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STK_BasicVector_h
#define STK_BasicVector_h

#include "common.h"
#include "Error.h"

#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#ifdef HAVE_ATLAS
extern "C"{
  #include <cblas.h>
  #include <clapack.h>
}
#endif


namespace STK
{
  template<typename _ElemT> class BasicVector;
  template<typename _ElemT> class Matrix;
  
  // we need to declare the friend << operator here
  template<typename _ElemT>
    std::ostream & operator << (std::ostream & out, BasicVector<_ElemT> & m);
    
  
  /** **************************************************************************
   ** **************************************************************************
   *  @brief Provides a matrix abstraction class
   *
   *  This class provides a way to work with matrices in STK.
   *  It encapsulates basic operations and memory optimizations.
   *
   */
  template<typename _ElemT>
      class BasicVector
    {
    public:
      BasicVector(): mpData(NULL)
#ifdef STK_MEMALIGN_MANUAL
        ,mpFreeData(NULL)
#endif
        , mLength(0)
      {}
      
      /**
       * @brief Copy constructor
       * @param rV 
       */
      BasicVector(const BasicVector<_ElemT>& rV);
      
      BasicVector(const size_t s);
      
      BasicVector(const _ElemT* pData, const size_t s);
      
      ~BasicVector();
      
      void
      Init(size_t length);
      
      /**
       * @brief Dealocates the window from memory and resets the dimensions to (0)
       */
      void
      Destroy();
           
      /**
       * @brief Returns @c true if vector is initialized
       */
      const bool
      IsInitialized() const
      { return mpData != NULL; }            

      /**
       * @brief Sets all elements to 0
       */
      void 
      Clear();
      
      const size_t
      Length() const
      { return mLength; }
      
      /**
       * @brief Returns size of matrix in memory (in bytes)
       */
      const size_t
      MSize() const
      {
        return (mLength + (((16 / sizeof(_ElemT)) - mLength%(16 / sizeof(_ElemT))) 
                          % (16 / sizeof(_ElemT)))) * sizeof(_ElemT);
      }
      
      /**
       *  @brief Gives access to the vector memory area
       *  @return pointer to the first field
       */
      _ElemT*
      pData() const
      { return mpData; }

      /**
       *  @brief Gives access to the vector memory area
       *  @return pointer to the first field
       */
      const _ElemT* const 
      cpData() const
      { return mpData; }
      
      /**
       *  @brief Gives access to a specified vector element without range check
       *  @return pointer to the first field of the row
       */
      _ElemT&      
      operator [] (size_t i) const
      { return *(mpData + i);}
      
      
      //########################################################################
      //########################################################################
      BasicVector<_ElemT>&
      AddCVMul(const _ElemT c, const BasicVector<_ElemT>& rV);
      
      BasicVector<_ElemT>&
      AddCVMul(const _ElemT c, const _ElemT* pV);

      BasicVector<_ElemT>&
      AddCVVDotMul(const _ElemT c, const _ElemT* pA, const size_t nA, 
                                   const _ElemT* pB, const size_t nB);
      
      /// this = this + c * rM * rV 
      BasicVector<_ElemT>&
      AddCMVMul(const _ElemT c, const Matrix<_ElemT>& rM, 
                const BasicVector<_ElemT>& rV);

      /// this = this + c * rM * pV 
      BasicVector<_ElemT>&
      AddCMVMul(const _ElemT c, const Matrix<_ElemT>& rM, 
                const _ElemT* pV);
                
      /// this = this + c * rM * rV 
      BasicVector<_ElemT>&
      AddCMtVMul(const _ElemT c, const Matrix<_ElemT>& rM, 
                 const BasicVector<_ElemT>& rV);

      /// this = this + c * rM * rV 
      BasicVector<_ElemT>&
      AddCMtVMul(const _ElemT c, const Matrix<_ElemT>& rM, 
                 const _ElemT* pV);

      /// Adds a diagonal after Matrix-Matrix Multiplication
      BasicVector<_ElemT>&
      AddDiagCMMMul(const _ElemT c, const Matrix<_ElemT>& rMa, 
                    const Matrix<_ElemT>& rMb);
      
      _ElemT
      Dot(const BasicVector<_ElemT>& rV);
      
      _ElemT
      Dot(const _ElemT* pV);
      
      
      /// Performs a row stack of the matrix rMa
      BasicVector<_ElemT>&
      MatrixRowStack(const Matrix<_ElemT>& rMa);
      
      /// Returns sum of the elements
      _ElemT
      Sum() const; 

      /// Returns sum of the elements
      BasicVector<_ElemT>&
      AddRowSum(const Matrix<_ElemT>& rM); 

      /// Returns sum of the elements
      BasicVector<_ElemT>&
      AddColSum(const Matrix<_ElemT>& rM); 

      /// Returns log(sum(exp())) without exp overflow
      _ElemT
      LogSumExp() const;

      //########################################################################
      //########################################################################
      
      friend std::ostream & 
      operator << <> (
        std::ostream& rOut, 
        BasicVector<_ElemT>& rV);
    
      //operator _ElemT* ()
      //{
      //  return mpData;
      //}
      
    //##########################################################################
    //##########################################################################        
    //protected:
    public:
      /// data memory area
      _ElemT*   mpData;
#ifdef STK_MEMALIGN_MANUAL
      /// data to be freed (in case of manual memalignment use, see common.h)
      _ElemT*   mpFreeData;
#endif
      size_t  mLength;      ///< Number of elements
    }; // class BasicVector

}; // namespace STK


//*****************************************************************************
//*****************************************************************************
// we need to include the implementation
#include "BasicVector.tcc"
//*****************************************************************************
//*****************************************************************************


/******************************************************************************
 ******************************************************************************
 * The following section contains specialized template definitions
 * whose implementation is in BasicVector.cc
 */
 
namespace STK
{
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const BasicVector<float>& rV);

  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCVMul(const double c, const BasicVector<double>& rV);
    
      
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const float* pV);

  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCVMul(const double c, const double* pV);
    
  
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMVMul(const float c, const Matrix<float>& rV, const BasicVector<float>& rV2);

  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMVMul(const double c, const Matrix<double>& rV, const BasicVector<double>& rV2);
    

  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMtVMul(const float c, const Matrix<float>& rV, const BasicVector<float>& rV2);

  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMtVMul(const double c, const Matrix<double>& rV, const BasicVector<double>& rV2);
    

  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMtVMul(const float c, const Matrix<float>& rV, const float* pV);

  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMtVMul(const double c, const Matrix<double>& rV, const double* pV);
    

  template<>
    float
    BasicVector<float>::
    Dot(const float* pV);

  template<>
    double
    BasicVector<double>::
    Dot(const double* pV);
    
} // namespace STK


#endif // #ifndef STK_BasicVector_h
