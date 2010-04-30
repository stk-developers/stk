/** 
 * @file Matrix.cc 
 * 
 * Implementation of specialized Matrix template methods 
 */

#include "Matrix.h"

#ifdef USE_MKL
  #include "mkl_lapack.h"  // replace CLAPACK with MKL LAPACK in case of Intel MKL
#endif

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
        fast_sigmoid_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    FastRowSigmoid()
    {
      for(size_t row = 0; row < this->Rows(); row++)
      {
        fast_sigmoid_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }
  

  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    RowSigmoid()
    {
      for(size_t row = 0; row < this->Rows(); row++)
      {
        sigmoid_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    RowSigmoid()
    {
      for(size_t row = 0; row < this->Rows(); row++)
      {
        sigmoid_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCVVtMul(const float& c, const BasicVector<float>& rA, 
                               const BasicVector<float>& rB)
    {
      assert(rA.Length() == this->mMRows);
      assert(rB.Length() == this->mMCols);
      
#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "AddCVVtMul(const float& c, const BasicVector<float>& rA, const BasicVector<float>& rB)" << std::endl;
#endif 

#ifdef HAVE_ATLAS
      cblas_sger(CblasRowMajor, rA.Length(), rB.Length(), c, rA.cpData(), 1,
                 rB.pData(), 1, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }

  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    AddCVVtMul(const double& c, const BasicVector<double>& rA, 
                                const BasicVector<double>& rB)
    {
      assert(rA.Length() == this->mMRows);
      assert(rB.Length() == this->mMCols);
      
#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "AddCVVtMul(const double& c, const BasicVector<double>& rA, const BasicVector<double>& rB)" << std::endl;
#endif 

#ifdef HAVE_ATLAS
      cblas_dger(CblasRowMajor, rA.Length(), rB.Length(), c, rA.cpData(), 1,
                   rB.cpData(), 1, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }
    

          
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    AddMMMul(const Matrix<float>& a, const Matrix<float>& b)
    { 
      assert(a.Cols() == b.Rows());
      assert(this->Rows() == a.Rows());
      assert(this->Cols() == b.Cols());
      
#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "AddMMMul(const Matrix<float>& a, const Matrix<float>& b)" << std::endl;
#endif 

#ifdef HAVE_ATLAS
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }; // AddMatMult(const ThisType & a, const ThisType & b)

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    AddMMMul(const Matrix<double>& a, const Matrix<double>& b)
    { 
      assert(a.Cols() == b.Rows());
      assert(this->Rows() == a.Rows());
      assert(this->Cols() == b.Cols());
      
#ifdef HAVE_ATLAS
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }; // AddMatMult(const ThisType & a, const ThisType & b)    
    
  
  template<>
    Matrix<float> &
    Matrix<float>::
    DiagScale(const BasicVector<float>& rDiagVector)
    {
      // :TODO: 
      // optimize this
      float*  data = mpData;
      for (size_t i=0; i < Rows(); i++)
      {
        for (size_t j=0; j < Cols(); j++)
        {
          data[j] *= rDiagVector[j];
        }
        data += mStride;
      }
      
      return *this;
    }
    
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    DiagScale(const float* pDiagVector)
    {
      // :TODO: 
      // optimize this
      float*  data = mpData;
      for (size_t i=0; i < Rows(); i++)
      {
        for (size_t j=0; j < Cols(); j++)
        {
          data[j] *= pDiagVector[j];
        }
        data += mStride;
      }
      
      return *this;
    }
  
    
      
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    AddMMtMul(const Matrix<float>& a, const Matrix<float>& b)
    { 
      assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());

#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "AddMMtMul(const Matrix<float>& a, const Matrix<float>& b)" << std::endl;
#endif 

#ifdef HAVE_ATLAS
      
      /*
      if (a.mpData[0] > 0) std::cerr << a.mpData[0]     << std::endl;
      std::cerr << a.mpData[0]     << std::endl;
      if (a.mStride)       std::cerr << a.mStride       << std::endl;
      std::cerr << a.mStride       << std::endl;
      if (b.mpData[0] > 0) std::cerr << b.mpData[0]     << std::endl;
      std::cerr << b.mpData[0]     << std::endl;
      if (b.mStride)       std::cerr << b.mStride       << std::endl;
      std::cerr << b.mStride       << std::endl;
      if (this->mpData[0] > 0) std::cerr << this->mpData[0]     << std::endl;
      std::cerr << this->mpData[0] << std::endl;
      if (this->mStride)       std::cerr << this->mStride       << std::endl;
      std::cerr << this->mStride   << std::endl;
      */
      
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Rows(), b.Cols(),
                    1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
  
      return *this;
    }; 

  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    AddMMtMul(const Matrix<double>& a, const Matrix<double>& b)
    { 
      assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());
#ifdef HAVE_ATLAS
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Rows(), b.Cols(),
                    1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
  
      return *this;
    }; 

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCMMtMul(const float& c, const Matrix<float>& a, 
                               const Matrix<float>& b)
    { 
      assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());

#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "AddCMMtMul(const float& c, const Matrix<float>& a, const Matrix<float>& b)" << std::endl;
#endif 


#ifdef HAVE_ATLAS
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Rows(), b.Cols(),
                    c, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
  
      return *this;
    }; 

  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    AddCMMtMul(const double& c, const Matrix<double>& a, 
                                const Matrix<double>& b)
    { 
      assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());
#ifdef HAVE_ATLAS
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Rows(), b.Cols(),
                    c, a.mpData, a.mStride, b.mpData, b.mStride, 1.0F, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
  
      return *this;
    }; 

                  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    RepMMMul(const Matrix<float>& a, const Matrix<float>& b)
    { 
    //fprintf(stderr, "A %d %d B %d %d C %d %d \n", a.Rows(), a.Cols(), b.Rows(), b.Cols(), Rows(), Cols());
      assert(a.Rows() == this->Rows());
      assert(b.Cols() == this->Cols());
      assert(a.Cols() == b.Rows());
      
#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "RepMMMul(const Matrix<float>& a, const Matrix<float>& b)" << std::endl;
#endif 
      
      Clear();
      
#ifdef HAVE_ATLAS
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      
      //cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, C->rows, C->cols, B->rows,
      //            1.0f, A->arr, A->realCols, B->arr, B->realCols, 1.0f, C->arr, C->realCols);		  
  
      return *this;
    }; 

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    Invert()
    { 
      assert(Rows() == Cols());

#if defined(HAVE_ATLAS) || defined(USE_MKL)
  #ifdef USE_MKL
      int* pivot = new int[mMRows];
      MKL_INT info;
      MKL_INT nrows = static_cast<MKL_INT>(Rows());
      MKL_INT ncols = static_cast<MKL_INT>(Cols());
      MKL_INT stride = static_cast<MKL_INT>(mStride);
      sgetrf(&nrows, &ncols, mpData, &stride, static_cast<MKL_INT *>(pivot), &info);
      MKL_INT lwork = static_cast<MKL_INT>(nrows * 16);  // mNRows * blocksize{16:64}, see MKL manual
      float *pwork_place = new float [lwork];
      sgetri(&nrows, mpData, &stride, static_cast<MKL_INT *>(pivot), pwork_place, &lwork, &info);
      delete [] pwork_place;
      delete [] pivot;
  #else
      int* pivot = new int[mMRows];
      clapack_sgetrf(CblasColMajor, Rows(), Cols(), mpData, mStride, pivot);
      clapack_sgetri(CblasColMajor, Rows(), mpData, mStride, pivot);
      delete [] pivot;
  #endif
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }; 
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    Invert()
    { 
      assert(Rows() == Cols());

#if defined(HAVE_ATLAS) || defined(USE_MKL)
  #ifdef USE_MKL
      int* pivot = new int[mMRows];
      MKL_INT info;
      MKL_INT nrows = static_cast<MKL_INT>(Rows());
      MKL_INT ncols = static_cast<MKL_INT>(Cols());
      MKL_INT stride = static_cast<MKL_INT>(mStride);
      dgetrf(&nrows, &ncols, mpData, &stride, static_cast<MKL_INT *>(pivot), &info);
      MKL_INT lwork = static_cast<MKL_INT>(nrows * 16);  // mNRows * blocksize{16:64}, see MKL manual
      double *pwork_place = new double [lwork];
      dgetri(&nrows, mpData, &stride, static_cast<MKL_INT *>(pivot), pwork_place, &lwork, &info);
      delete [] pwork_place;
      delete [] pivot;
  #else
      int* pivot = new int[mMRows];
      clapack_dgetrf(CblasColMajor, Rows(), Cols(), mpData, mStride, pivot);
      clapack_dgetri(CblasColMajor, Rows(), mpData, mStride, pivot);
      delete [] pivot;
  #endif
#else
      Error("Method not implemented without BLAS");
#endif
      return *this;
    }; 

/*  
    STK::Matrix<FLOAT> x(5, 5);
    x[0][0] = 0.9901;
    x[0][1] = 0.6435;
    x[0][2] = 0.7446;
    x[0][3] = 0.2126;
    x[0][4] = 0.6072;
    
    x[1][0] = 0.7889;
    x[1][1] = 0.3200;
    x[1][2] = 0.2679;
    x[1][3] = 0.8392;
    x[1][4] = 0.6299;
    
    x[2][0] = 0.4387;
    x[2][1] = 0.9601;
    x[2][2] = 0.4399;
    x[2][3] = 0.6288;
    x[2][4] = 0.3705;
    
    x[3][0] = 0.4983;
    x[3][1] = 0.7266;
    x[3][2] = 0.9334;
    x[3][3] = 0.1338;
    x[3][4] = 0.5751;
    
    x[4][0] = 0.2140;
    x[4][1] = 0.4120;
    x[4][2] = 0.6833;
    x[4][3] = 0.2071;
    x[4][4] = 0.4514;
    
    for (int i = 0; i < x.Rows(); i++)
    {
      for (int j = 0; j < x.Cols(); j++)
      {
        cout << x[i][j] << " ";
      }
      cout << endl;
    }
    
    x.Invert();
    
    for (int i = 0; i < x.Rows(); i++)
    {
      for (int j = 0; j < x.Cols(); j++)
      {
        cout << x[i][j] << " ";
      }
      cout << endl;
}*/
        
        
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    RepMtMMul(const Matrix<float> & a, const Matrix<float> & b)
    { 
      /*assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());
      */
      Clear();
            
#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "RepMtMMul(const Matrix<float> & a, const Matrix<float> & b)" << std::endl;
#endif 
      
      //cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, this->Rows(), this->Cols(), b.Rows(),
      //          1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
      //printf("%d %d %d %d %d %d\n", a.Cols(), b.Cols(), b.Rows(), a.mStride, b.mStride, this->mStride);
#ifdef HAVE_ATLAS
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, a.Cols(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      
      return *this;
    }; 

  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    RepMtMMul(const Matrix<double> & a, const Matrix<double> & b)
    { 
      /*assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());
      */
      Clear();
            
      //cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, this->Rows(), this->Cols(), b.Rows(),
      //          1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
      //printf("%d %d %d %d %d %d\n", a.Cols(), b.Cols(), b.Rows(), a.mStride, b.mStride, this->mStride);
#ifdef HAVE_ATLAS
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, a.Cols(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      
      return *this;
    }; 

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    AddCMtMMul(const float& c, const Matrix<float>& a, 
                               const Matrix<float>& b)
    { 
      /*assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());
      */
      
#ifdef TRACE_MATRIX_OPERATIONS
      std::cerr << "AddCMtMMul(const float& c, const Matrix<float>& a, const Matrix<float>& b)" << std::endl;
#endif 
      
      //cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, this->Rows(), this->Cols(), b.Rows(),
      //          1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
      //printf("%d %d %d %d %d %d\n", a.Cols(), b.Cols(), b.Rows(), a.mStride, b.mStride, this->mStride);
#ifdef HAVE_ATLAS
      cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, a.Cols(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, c, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      
      return *this;
    }; 

  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    AddCMtMMul(const double& c, const Matrix<double>& a, 
                                const Matrix<double>& b)
    { 
      /*assert(a.Rows() == this->Rows());
      assert(b.Rows() == this->Cols());
      assert(a.Cols() == b.Cols());
      */
#ifdef HAVE_ATLAS
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, a.Cols(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, c, this->mpData, this->mStride);
#else
      Error("Method not implemented without BLAS");
#endif
      
      return *this;
    }; 
    
          
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    FastRowSoftmax()
    {
      for(size_t row = 0; row < this->Rows(); row++){
        fast_softmax_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }
          
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    FastRowSoftmax()
    {
      for(size_t row = 0; row < this->Rows(); row++){
        fast_softmax_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }

  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<float> &
    Matrix<float>::
    RowSoftmax()
    {
      for(size_t row = 0; row < this->Rows(); row++){
        softmax_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }
          
  
  //***************************************************************************
  //***************************************************************************
  template<>
    Matrix<double> &
    Matrix<double>::
    RowSoftmax()
    {
      for(size_t row = 0; row < this->Rows(); row++){
        softmax_vec((*this)[row], (*this)[row], this->Cols());
      }
      return *this;
    }

}; //namespace STK
