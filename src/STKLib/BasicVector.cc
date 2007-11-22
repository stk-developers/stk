#include "Matrix.h"
#include "BasicVector.h"


namespace STK
{  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const BasicVector<float>& rV)
    {
#ifdef HAVE_ATLAS
      cblas_saxpy(mLength, c, rV.mpData, 1, mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
  
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCVMul(const double c, const BasicVector<double>& rV)
    {
#ifdef HAVE_ATLAS
      cblas_daxpy(mLength, c, rV.mpData, 1, mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
  
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const float* pV)
    {
#ifdef HAVE_ATLAS
      cblas_saxpy(mLength, c, pV, 1, mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCVMul(const double c, const double* pV)
    {
#ifdef HAVE_ATLAS
      cblas_daxpy(mLength, c, pV, 1, mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
        
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMVMul(const float c, const Matrix<float>& rM, const BasicVector<float>& rV)
    {
#ifdef HAVE_ATLAS
      cblas_sgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), rV.cpData(), 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMVMul(const double c, const Matrix<double>& rM, const BasicVector<double>& rV)
    {
#ifdef HAVE_ATLAS
      cblas_dgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), rV.cpData(), 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
    
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMtVMul(const float c, const Matrix<float>& rM, const BasicVector<float>& rV)
    {
#ifdef HAVE_ATLAS
      cblas_sgemv(CblasRowMajor, CblasTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), rV.cpData(), 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMtVMul(const double c, const Matrix<double>& rM, const BasicVector<double>& rV)
    {
#ifdef HAVE_ATLAS
      cblas_dgemv(CblasRowMajor, CblasTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), rV.cpData(), 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
    
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMtVMul(const float c, const Matrix<float>& rM, const float* pV)
    {
#ifdef HAVE_ATLAS
      cblas_sgemv(CblasRowMajor, CblasTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), pV, 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMtVMul(const double c, const Matrix<double>& rM, const double* pV)
    {
#ifdef HAVE_ATLAS
      cblas_dgemv(CblasRowMajor, CblasTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), pV, 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
    
  
  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCMVMul(const float c, const Matrix<float>& rM, const float* pV)
    {
#ifdef HAVE_ATLAS
      cblas_sgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), pV, 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }

  //****************************************************************************
  //****************************************************************************
  template<>
    BasicVector<double>&
    BasicVector<double>::
    AddCMVMul(const double c, const Matrix<double>& rM, const double* pV)
    {
#ifdef HAVE_ATLAS
      cblas_dgemv(CblasRowMajor, CblasNoTrans, rM.Rows(), rM.Cols(), c, rM.pData(), 
                  rM.Stride(), pV, 1, 1.0F, this->mpData, 1);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
        

  //****************************************************************************
  //****************************************************************************
  template<>
    float
    BasicVector<float>::
    Dot(const float* pV)
    {
#ifdef HAVE_ATLAS
      return cblas_sdot(mLength, mpData, 1, pV, 1);
#else
      Error("Method not implemented without BLAS");
      return 0.0;
#endif
    }
  
  //****************************************************************************
  //****************************************************************************
  template<>
    double
    BasicVector<double>::
    Dot(const double* pV)
    {
#ifdef HAVE_ATLAS
      return cblas_ddot(mLength, mpData, 1, pV, 1);
#else
      Error("Method not implemented without BLAS");
      return 0.0;
#endif
    }


}

