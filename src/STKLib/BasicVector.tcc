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

#ifdef USE_BLAS
extern "C"{
  #include <cblas.h>
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
      size_t skip;
      size_t stride;      // leading dimension
      size_t size;
      void*  data;
      void*  free_data;
      
      // compute the size in memory
      skip   = ((16 / sizeof(_ElemT)) - length % (16 / sizeof(_ElemT))) 
             % (16 / sizeof(_ElemT));
      stride = length + skip;
      size   =  stride * sizeof(_ElemT);
      
      if (NULL != (data = stk_memalign(16, size, &free_data)))
      {
        mpData        = static_cast<_ElemT*> (data);
#ifdef STK_MEMALIGN_MANUAL
        mpFreeData    = static_cast<_ElemT*> (free_data);
#endif
        mLength = length;
        //mStride = stride;
        
        // set all bytes to 0
        memset(mpData, 0, size);
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
      if (NULL != mpFreeData) free(mpData);
#else
      if (NULL != mpData) free(mpData);
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
      memset(mpData, 0, mLength*sizeof(_ElemT)); 
    }
      
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    BasicVector<_ElemT>&
    BasicVector<_ElemT>::
    AddCVMul(const _ElemT c, const _ElemT* pV, const size_t nV)
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
    std::ostream &
    operator << (std::ostream& rOut, BasicVector<_ElemT>& rV)
    {
      for (size_t i = 0; i < rV.Length(); i++)
      {
        rOut << rV[i] << " ";
      }
      rOut << std::endl;
      
      return rOut;
    }
    
} // namespace STK


#endif // STK_BasicVector_tcc

