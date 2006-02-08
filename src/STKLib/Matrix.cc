/** 
 * @file Matrix.cc 
 * 
 * Implementation of specialized Matrix template methods 
 */

#include "Matrix.h"

extern "C" 
{
# include <cblas.h>
}


namespace STK
{
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSigmoid()
    {
      for(size_t row = 0; row < this->Rows(); row++)
      {
        fast_sigmoid_vec(this->Row(row), this->Row(row), this->Cols());
      }
      return *this;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    AddMatMult(Matrix<float> & a, Matrix<float> & b)
    { 
      assert(a.Cols() == b.Rows());
      assert(this->Rows() == a.Rows());
      assert(this->Cols() == b.Cols());
      assert((int)a() % 16 == 0);
      assert((int)b() % 16 == 0);
      assert((int)this->mpData % 16 == 0);
      if(b.Storage() == STORAGE_TRANSPOSED){
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Cols(), b.Rows(),
                    1.0f, a.mpData, a.mMRealCols, b.mpData, b.mMRealCols, 1.0f, this->mpData, this->mMRealCols);
      }
      else{
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.Rows(), b.Cols(), b.Rows(),
                    1.0f, a.mpData, a.mMRealCols, b.mpData, b.mMRealCols, 1.0f, this->mpData, this->mMRealCols);
      }
      return *this;
    }; // AddMatMult(const ThisType & a, const ThisType & b)

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSoftmax()
    {
      for(size_t row = 0; row < this->Rows(); row++){
        fast_softmax_vec(this->Row(row), this->Row(row), this->Cols());
      }
      return *this;
    }
    
  
}; //namespace STK
