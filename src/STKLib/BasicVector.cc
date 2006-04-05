#include "Matrix.h"
#include "BasicVector.h"


namespace STK
{  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const float* pV)
    {
      cblas_saxpy(mLength, c, pV, 1, mpData, 1);
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCVMul(const double c, const double* pV)
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
                  rM.Stride(), rV.cpData(), 1, 1.0F, this->mpData, 1);
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
                  rM.Stride(), rV.cpData(), 1, 1.0F, this->mpData, 1);
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
                  rM.Stride(), pV, 1, 1.0F, this->mpData, 1);
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
                  rM.Stride(), pV, 1, 1.0F, this->mpData, 1);
      return *this;
    }
        

  //****************************************************************************
  //****************************************************************************
  template<>
    float
    BasicVector<float>::
    Dot(const float* pV)
    {
      return cblas_sdot(mLength, mpData, 1, pV, 1);
    }
  
  //****************************************************************************
  //****************************************************************************
  template<>
    double
    BasicVector<double>::
    Dot(const double* pV)
    {
      return cblas_ddot(mLength, mpData, 1, pV, 1);
    }


}

