#include "Matrix.h"
#include "BasicVector.h"


namespace STK
{  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const float* pV, const size_t nV)
    {
      cblas_saxpy(mLength, c, pV, 1, mpData, 1);
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCVMul(const double c, const double* pV, const size_t nV)
    {
      cblas_daxpy(mLength, c, pV, 1, mpData, 1);
      return *this;
    }
        
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMVMul(const float c, const Matrix<float>& rM, const BasicVector<float>& rV)
    {
      cblas_sgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), rV.pData(), 1, 1.0F, this->pData(), 1);
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMVMul(const double c, const Matrix<double>& rM, const BasicVector<double>& rV)
    {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), rV.pData(), 1, 1.0F, this->pData(), 1);
      return *this;
    }
    
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMVMul(const float c, const Matrix<float>& rM, const float* pV)
    {
      cblas_sgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), pV, 1, 1.0F, this->pData(), 1);
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMVMul(const double c, const Matrix<double>& rM, const double* pV)
    {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), pV, 1, 1.0F, this->pData(), 1);
      return *this;
    }
        
}
