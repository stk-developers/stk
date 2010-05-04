/** @file BasicVector.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef STK_BasicVector_tcc
#define STK_BasicVector_tcc

#include "common.h"
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <cstring>

#ifdef HAVE_ATLAS
extern "C"{
  #ifdef USE_MKL
    #include "mkl.h"
  #else
    #include <cblas.h>
  #endif
}
#endif

#include<fstream>
#include<iomanip>
 
namespace STK
{
 //******************************************************************************    
  template<typename _ElemT>
    void
    BasicVector<_ElemT>::
    Init(const size_t length)
    {
      size_t size;
      void*  data;
      void*  free_data;
      
      size = align<16>(length * sizeof(_ElemT));
      
      if (NULL != (data = stk_memalign(16, size, &free_data)))
      {
        mpData        = static_cast<_ElemT*> (data);
#ifdef STK_MEMALIGN_MANUAL
        mpFreeData    = static_cast<_ElemT*> (free_data);
#endif
        mLength = length;
        
        // set all bytes to 0
        std::memset(mpData, 0, size);
      }
      else
      {
        throw std::bad_alloc();
      }
    }
    
  
  //******************************************************************************    
  //******************************************************************************    
  template<typename _ElemT>
    BasicVector<_ElemT>::
    BasicVector(const size_t s)
    {
      Init(s);
    }
    
    
  //******************************************************************************    
  //******************************************************************************    
  template<typename _ElemT>
    BasicVector<_ElemT>::
    BasicVector(const BasicVector<_ElemT>& rV)
    {
      size_t size;
      void*  data;
      void*  free_data;
      
      size   =  rV.MSize();

#ifdef STK_MEMALIGN_MANUAL
//      if (NULL != mpFreeData) free(mpData);
#else
//      if (NULL != mpData) free(mpData);
#endif 


      if (NULL != (data = stk_memalign(16, size, &free_data)))
      {
        mpData        = static_cast<_ElemT*> (data);
#ifdef STK_MEMALIGN_MANUAL
        mpFreeData    = static_cast<_ElemT*> (free_data);
#endif
        mLength = rV.mLength;
        // copy the memory block
        memcpy(this->mpData, rV.mpData, size);
      }
      else
      {
        mpData      = NULL;
#ifdef STK_MEMALIGN_MANUAL
        mpFreeData  = NULL;
#endif        
        throw std::bad_alloc();
      }
    }


  //****************************************************************************
  //****************************************************************************
  // The destructor
  template<typename _ElemT>
    void
    BasicVector<_ElemT>::
    Destroy()
    {
      // we need to free the data block if it was defined
#ifndef STK_MEMALIGN_MANUAL
      if (NULL != mpData) free(mpData);
#else
      if (NULL != mpFreeData) free(mpFreeData);
      mpFreeData = NULL;
#endif

      mpData = NULL;
      mLength = 0;
    }

  
  
  //******************************************************************************    
  //******************************************************************************    
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    AddDiagCMMMul(const _ElemT c, const Matrix<_ElemT>& rMa, 
                  const Matrix<_ElemT>& rMb)
    {
      _ElemT tmp_val;
      
      for (size_t i = 0; i < mLength; i++)
      {
        tmp_val = 0.0;
        
        for (size_t j = 0; j < rMa.Cols(); j++)
        {
          tmp_val += rMa[i][j] * rMb[j][i];
        }
        
        mpData[i] += tmp_val * c;
      }
      
      return *this;
    }
    
  
  //******************************************************************************    
  //******************************************************************************    
  template<typename _ElemT>
    BasicVector<_ElemT>::
    BasicVector(const _ElemT* pData, const size_t s)
    {
      size_t size(align<16>(s*sizeof(_ElemT)));
      void*  data;
      void*  free_data;

      if (NULL != (data = stk_memalign(16, size, &free_data)))
      {
        mpData        = static_cast<_ElemT*> (data);
#ifdef STK_MEMALIGN_MANUAL
        mpFreeData    = static_cast<_ElemT*> (free_data);
#endif
        mLength = s;
        
        // set all bytes to 0
        std::memset(mpData, 0, size);
        // copy the memory block
        memcpy(this->mpData, pData, s);
      }
      else
      {
        mpData      = NULL;
#ifdef STK_MEMALIGN_MANUAL
        mpFreeData  = NULL;
#endif        
        throw std::bad_alloc();
      }
    }
    
  
  //****************************************************************************
  //****************************************************************************
  // The destructor
  template<typename _ElemT>
    // virtual
    BasicVector<_ElemT>::
    ~BasicVector()
    {
      // we need to free the data block if it was defined
#ifndef STK_MEMALIGN_MANUAL
      if (NULL != mpData) free(mpData);
#else
      if (NULL != mpFreeData) free(mpFreeData);
#endif
    }

  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    void 
    BasicVector<_ElemT>::
    Clear()      
    { 
      std::memset(mpData, 0, mLength*sizeof(_ElemT)); 
    }
      
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    AddCVMul(const _ElemT c, const _ElemT* pV)
    {
      Error("BasicMatrix::AddCVMul(const _ElemT,const _ElemT*,const size_t) not implemented");
      return *this;
    }
      
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    AddCVVDotMul(const _ElemT c, const _ElemT* pA, const size_t nA, 
                                 const _ElemT* pB, const size_t nB)
    {
      for (size_t i = 0; i < mLength; i++)
      {
        mpData[i] += c * pA[i] * pB[i];
      }
      
      return *this;
    }

  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    MatrixRowStack(const Matrix<_ElemT>& rMa)
    {
      assert(mLength == rMa.Cols() * rMa.Rows());
      
      _ElemT*       inc_data = mpData;
      const size_t  cols     = rMa.Cols();
      
      for (size_t i = 0; i < rMa.Rows(); i++)
      {
        // copy the data to the propper position
        memcpy(inc_data, rMa[i], cols * sizeof(_ElemT));
        
        // set new copy position
        inc_data += cols;
      }
    }
    
    
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    _ElemT
    BasicVector<_ElemT>::
    Sum() const
    {
      double sum = 0.0;

      for (size_t i = 0; i < mLength; ++i) {
        sum += mpData[i];
      } 

      return sum;
    }


  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    AddColSum(const Matrix<_ElemT>& rM)
    {
      // note the double accumulator
      double sum;

      assert(mLength == rM.Cols());

      for (size_t i = 0; i < mLength; ++i) {
        sum = 0.0;
        for (size_t j = 0; j < rM.Rows(); ++j) {
          sum += rM[j][i];
        }  
        mpData[i] += sum;
      } 

      return *this;
    }


  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    AddRowSum(const Matrix<_ElemT>& rM)
    {
      // note the double accumulator
      double sum;

      assert(mLength == rM.Rows());

      for (size_t i = 0; i < mLength; ++i) {
        sum = 0.0;
        for (size_t j = 0; j < rM.Cols(); ++j) {
          sum += rM[i][j];
        }  
        mpData[i] += sum;
      } 

      return *this;
    }


  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    _ElemT
    BasicVector<_ElemT>::
    LogSumExp() const
    {
      if (mLength > 0) {
        double sum = mpData[0];

        for (size_t i = 1; i < mLength; ++i) {
          sum = LogAdd(sum, mpData[i]);
        } 

        return sum;
      }
      else {
        return LOG_0;
      }
    }


  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    std::ostream &
    operator << (std::ostream& rOut, BasicVector<_ElemT>& rV)
    {
      for (size_t i = 0; i < rV.Length(); i++)
      {
        rOut << rV[i] << " ";
      }
      
      return rOut;
    }
    
} // namespace STK


#endif // STK_BasicVector_tcc

