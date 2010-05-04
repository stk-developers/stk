
/** @file Matrix.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */


#ifndef STK_Matrix_tcc
#define STK_Matrix_tcc

//#pragma GCC system_header

#include "common.h"
#include "mymath.h"
#define _XOPEN_SOURCE 600
#include <cstdlib>
#include <cassert>
#include <math.h>

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
#include<typeinfo>

namespace STK
{
  
  
    
  //****************************************************************************
  // Basic constructor
  template<typename _ElemT>
    Matrix<_ElemT>::
    Matrix(const size_t r, const size_t       c)
    {
      Init(r, c);
    }


//******************************************************************************    
  template<typename _ElemT>
  void
  Matrix<_ElemT>::
  Init(const size_t rows,
       const size_t cols)
  {
    // initialize some helping vars
    size_t  skip;
    size_t  real_cols;
    size_t  size;
    void*   data;       // aligned memory block
    void*   free_data;  // memory block to be really freed

    // compute the size of skip and real cols
    skip      = ((16 / sizeof(_ElemT)) - cols % (16 / sizeof(_ElemT))) % (16 / sizeof(_ElemT));
    real_cols = cols + skip;
    size      = rows * real_cols * sizeof(_ElemT);

    // allocate the memory and set the right dimensions and parameters

    if (NULL != (data = stk_memalign(16, size, &free_data)))
    {
      this->mpData        = static_cast<_ElemT *> (data);
#ifdef STK_MEMALIGN_MANUAL
      this->mpFreeData    = static_cast<_ElemT *> (free_data);
#endif
      this->mMRows      = rows;
      this->mMCols      = cols;
      this->mStride  = real_cols;

      // set all bytes to 0
      memset(this->mpData, 0, this->MSize());
    }
    else
    {
      throw std::bad_alloc();
    }
  } //
  
  //****************************************************************************
  //****************************************************************************
  // The destructor
  template<typename _ElemT>
    void
    Matrix<_ElemT>::
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
      mMRows = mMCols = 0;
    }

  //****************************************************************************
  //****************************************************************************
  // Copy constructor
  template<typename _ElemT>
    Matrix<_ElemT>::
    Matrix(const ThisType & t)
    {
      if (NULL != t.mpData)
      {
        void* data;
        void* free_data;
  
        // first allocate the memory
        // this->mpData = new char[t.M_size];
        // allocate the memory and set the right dimensions and parameters
        if (NULL != (data = stk_memalign(16, t.MSize(), &free_data)))
        {
          this->mpData        = static_cast<_ElemT *> (data);
#ifdef STK_MEMALIGN_MANUAL
          this->mpFreeData    = static_cast<_ElemT *> (free_data);
#endif
          // copy the memory block
          memcpy(this->mpData, t.mpData, t.MSize());
  
          // set the parameters
          this->mMRows       = t.mMRows;
          this->mMCols       = t.mMCols;
          this->mStride      = t.mStride;
        }
        else
        {
          // throw bad_alloc exception if failure
          throw std::bad_alloc();
        }
      }
      else
      {
        this->mpData        = NULL;
#ifdef STK_MEMALIGN_MANUAL
        this->mpFreeData    = NULL;
#endif
        // set the parameters
        this->mMRows       = 0;
        this->mMCols       = 0;
        this->mStride      = 0;
      }
    }


  
  //****************************************************************************
  //****************************************************************************
  // The destructor
  template<typename _ElemT>
    Matrix<_ElemT>::
    ~Matrix()
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
    Matrix<_ElemT>::
    swap4b(void* a)
    {
      char* b = (char *)a;
      char c;
      c = b[0]; b[0] = b[3]; b[3] = c;
      c = b[1]; b[1] = b[2]; b[2] = c;
    }

  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    void
    Matrix<_ElemT>::
    swap2b(void* a)
    {
      char *b = (char *)a;
      char c;
      c = b[0]; b[0] = b[1]; b[1] = c;
    }

  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    bool
    Matrix<_ElemT>::
    LoadHTK(const char* pFileName)
    {
      HtkHeader htk_hdr;

      FILE *fp = fopen(pFileName, "rb");
      if(!fp)
      {
        //fprintf(stderr, "Mat<T>::loadHTK - can not open the file: file\n", file);
        //exit(1);
        return false;
      }

      read(fileno(fp), &htk_hdr, sizeof(htk_hdr));

      swap4b(&htk_hdr.nSamples);
      swap4b(&htk_hdr.sampPeriod);
      swap2b(&htk_hdr.sampSize);
      swap2b(&htk_hdr.paramKind);

      Init(htk_hdr.nSamples, htk_hdr.sampSize / sizeof(float));

      size_t i;
      size_t j;
      if (typeid(_ElemT) == typeid(float))
      {
        for (i=0; i< Rows(); ++i) {
          read(fileno(fp), (*this)[i], Cols() * sizeof(float));

          for(j = 0; j < Cols(); j++) {
            swap4b(&((*this)[i][j]));
          }
        }
      }
      else
      {
        float *pmem = new (std::nothrow) float[Cols()];
        if (!pmem)
        {
          fclose(fp);
          return false;
        }

        for(i = 0; i < Rows(); i++) {
          read(fileno(fp), pmem, Cols() * sizeof(float));

          for (j = 0; j < Cols(); ++j) {
            swap4b(&pmem[j]);
            (*this)[i][j] = static_cast<_ElemT>(pmem[j]);
          }
        }
        delete [] pmem;
      }

      fclose(fp);
      
      return true;
    }
    
      
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    DivC(const _ElemT& c)
    {
      for (size_t i = 0; i < mMRows * mStride; i++)
      {
        mpData[i] = mpData[i] / c;
      }
      
      return *this;
    }
    
      
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    AddCVVtMul(const _ElemT& c, const BasicVector<_ElemT>& rA, 
               const BasicVector<_ElemT>& rB)
    {
      for (size_t i = 0; i < rA.Length(); i++)
      {
        for (size_t j = 0; j < rB.Length(); j++)
        {
          (*this)[i][j] += rA[i] * rB[j] * c;
        }
      }
    }
          
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    AddMMMul(const ThisType& a, const ThisType& b)
    { 
      if(!(a.Cols() == b.Rows() && this->Rows() == a.Rows() &&  this->Cols() == b.Cols()))
        STK::Error("Matrix multiply: bad matrix sizes (%d %d)*(%d %d) -> (%d %d)", a.Rows(), a.Cols(), b.Rows(), b.Cols(), this->Rows(), this->Cols());
      
      // :KLUDGE: Dirty for :o)        
      for(int r = 0; r < this->Rows() ; r++){
        for(int c = 0; c < this->Cols(); c++){
          // for every out matrix cell
          for(int index = 0; index < a.Cols(); index++){
            (*this)(r, c) += a(r, index) * b(index, c);
            //std::cout << a(r, index);
          }
        }
      }
      
      STK::Warning("Could be slow, not well implemented");
      
      return *this;
    }; // AddMatMult(const ThisType & a, const ThisType & b)
  
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    RepMMMul(const Matrix<_ElemT>& a, const Matrix<_ElemT>& b)
    { 
      Clear();
      AddMMMul(a, b);
      
      return *this;
    }; // AddMatMult(const ThisType & a, const ThisType & b)
  
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    AddCMMul(const _ElemT& c, const Matrix<_ElemT>& a) 
    { 
      assert(this->Cols() == a.Cols());
      assert(this->Rows() == a.Rows());
      
      _ElemT* pA = NULL;
      _ElemT* pT = NULL;
      
      for(unsigned row = 0; row < this->Rows(); row++)
      {
        pA = a[row];
        pT = (*this)[row];
        for(unsigned col = 0; col < this->Cols(); col++)
        {
          //(*this)(row, col) += a(row, col) * c;
          (*pT) += (*pA) * c;
          pT++;
          pA++;
        }
      }
      
      return *this;
    };
  
  
  //******************************************************************************
  //******************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    AddCVVtMul(const _ElemT& c, const _ElemT* pA, const _ElemT* pB)
    {
      for (size_t i= 0; i < Rows(); i++)
      {
        for (size_t j = 0; j < Cols(); j++)
        {
          (*this)[i][j] = c * pA[i] * pB[j];
        }
      }
      
      return *this;
    }
  
  
  //******************************************************************************
  //******************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    RepMMSub(const ThisType& a, const ThisType& b)
    { 
      assert(this->Cols() == a.Cols());
      assert(this->Rows() == a.Rows());
      assert(this->Cols() == b.Cols());
      assert(this->Rows() == b.Rows());      
      
      _ElemT* pA = NULL;
      _ElemT* pB = NULL;
      _ElemT* pT = NULL;
      for(unsigned row = 0; row < this->Rows(); row++)
      {
        pA = a[row];
        pB = b[row];
        pT = (*this)[row];
        for(unsigned col = 0; col < this->Cols(); col++)
        {
          //(*this)(row, col) = a(row, col) - b(row, col);
          (*pT) = (*pA) - (*pB);
          pT++;
          pA++;
          pB++;
        }
      }
      return *this;
    }; 
  
//******************************************************************************
/*
  template<>
  inline
  Matrix<float> &
  Matrix<float>::
  AddMatMult(Matrix<float> & a, Matrix<float> & b)
  { 
    assert(a.Cols() == b.Rows());
    assert(this->Rows() == a.Rows());
    assert(this->Cols() == b.Cols());
    if(b.Storage() == STORAGE_TRANSPOSED){
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
    }
    else{
      cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.Rows(), b.Cols(), b.Rows(),
                  1.0f, a.mpData, a.mStride, b.mpData, b.mStride, 1.0f, this->mpData, this->mStride);
    }
    return *this;
  }; // AddMatMult(const ThisType & a, const ThisType & b)
*/

  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    DotMul(const ThisType& a)
    {
      size_t i;
      size_t j;

      for (i = 0; i < mMRows; ++i) {
        for (j = 0; j < mMCols; ++j) {
          (*this)[i][j] *= a[i][j];
        }
      }
      return *this;
    }
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    const _ElemT &
    Matrix<_ElemT>::
    Sum() const
    {
      double sum = 0.0;
      
      for (size_t i = 0; i < Rows(); ++i) {
        for (size_t j = 0; j < Cols(); ++i) {
          sum += (*this)[i][j];
        }
      }

      return sum; 
    }

  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    NormalizeRows()
    {
      size_t i;
      size_t j;

      for (i = 0; i < mMRows; ++i) {
        double sum = 0.0;

        for (j = 0; j < mMCols; ++j) {
          sum += (*this)[i][j];
        }

        for (j = 0; j < mMCols; ++j) {
          (*this)[i][j] /= sum;
        }
      }
      return *this;
    }
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    Log()
    {
      size_t i;
      size_t j;

      for (i = 0; i < mMRows; ++i) {
        for (j = 0; j < mMCols; ++j) {
          if ((*this)[i][j] > DBL_EPSILON) {
            (*this)[i][j] = my_log((*this)[i][j]);
          }
          else if ((*this)[i][j] < -DBL_EPSILON) {
            STK::Error("Sigmoid only implemented for float");
          }
          else {
            (*this)[i][j] = LOG_0;
          }
        }
      }
      return *this;
    }
  

  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    FastRowSigmoid()
    {
      STK::Error("Sigmoid only implemented for float");
      return *this;
    }
  

//******************************************************************************
  /*
  template<>
  inline
  Matrix<float> &
  Matrix<float>::
  FastRowSigmoid()
  {
    for(size_t row = 0; row < this->Rows(); row++)
    {
      sigmoid_vec(this->Row(row), this->Row(row), this->Cols());
    }
    return *this;
  }
  */
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    FastRowSoftmax()
    {
      STK::Error("Softmax only implemented for float");
      return *this;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    DiagScale(const _ElemT* pDiagVector)
    {
      // :TODO: 
      // optimize this
      _ElemT* data = mpData;
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
  template<typename _ElemT>
    Matrix<_ElemT>&
    Matrix<_ElemT>::
    DiagScale(const BasicVector<_ElemT>& rDiagVector)
    {
      // :TODO: 
      // optimize this
      _ElemT* data = mpData;
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
        
      
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    Matrix<_ElemT> &
    Matrix<_ElemT>::
    Clear()
    {
      memset(mpData, 0, MSize());
      return *this;
    }
  
  
  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    void 
    Matrix<_ElemT>::
    PrintOut(char* file)
    {
      FILE* f = fopen(file, "w");
      unsigned i,j;
      fprintf(f, "%dx%d\n", this->mMRows, this->mMCols);
      
      for(i=0; i<this->mMRows; i++)
      {
        _ElemT*   row = (*this)[i];
        
        for(j=0; j<this->mStride; j++){
          WriteNumber(f, row[j], 20, 17);
        }
        fprintf(f, "\n");
      }
      
      fclose(f);
    }

      template<typename _ElemT>
    void 
    Matrix<_ElemT>::
    ReadIn(char* file)
    {
      FILE* f = fopen(file, "r");
      int  i = 0;
      int j = 0;
      fscanf(f, "%dx%d\n", &i,&j);
      fprintf(stderr, "%dx%d\n", i,j);
      
      for(i=0; i<this->mMRows; i++)
      {
        _ElemT*   row = (*this)[i];
        
        for(j=0; j<this->mStride; j++){
          ReadNumber(f, &row[j]);
          fscanf(f, " ");
        }
        //fprintf(f, "\n");
      }
      
      fclose(f);
    }
  
  //****************************************************************************
  //****************************************************************************
  // Copy constructor
  template<typename _ElemT>
    WindowMatrix<_ElemT>::
    WindowMatrix(const ThisType & t) :
      Matrix<_ElemT> (t)
    {
      // we need to be sure that the storage type is defined
      if (t.mpData != NULL)
      {
        // the data should be copied by the parent constructor
        void* data;
        void* free_data;

        // first allocate the memory
        // this->mpData = new char[t.M_size];
        // allocate the memory and set the right dimensions and parameters
        if (NULL != (data = stk_memalign(16, t.MSize(), &free_data)))
        {
          Matrix<_ElemT>::mpData        = static_cast<_ElemT *> (data);
#ifdef STK_MEMALIGN_MANUAL
          Matrix<_ElemT>::mpFreeData    = static_cast<_ElemT *> (free_data);
#endif
          // copy the memory block
          memcpy(this->mpData, t.mpData, t.MSize());

          // set the parameters
          this->mMRows       = t.mMRows;
          this->mMCols       = t.mMCols;
          this->mStride   = t.mStride;
        }
        else
        {
          // throw bad_alloc exception if failure
          throw std::bad_alloc();
        }
      }
    }


//******************************************************************************
  template<typename _ElemT>
    void
    WindowMatrix<_ElemT>::
    SetSize(const size_t rowOff,
            const size_t rows)
    {
      // point to the begining of window
      Matrix<_ElemT>::mpData = mpOrigData + mOrigMRealCols * rowOff;
      Matrix<_ElemT>::mMRows = rows;
      Matrix<_ElemT>::mMCols = mOrigMCols;
      Matrix<_ElemT>::mStride = mOrigMRealCols;
    }


  //****************************************************************************
  //****************************************************************************
  template<typename _ElemT>
    void
    WindowMatrix<_ElemT>::
    SetSize(const size_t rowOff,
            const size_t rows,
            const size_t colOff,
            const size_t cols)
    {
      // point to the begining of window
      Matrix<_ElemT>::mpData = mpOrigData + mOrigMRealCols * rowOff + cols;
      Matrix<_ElemT>::mMRows = rows;
      Matrix<_ElemT>::mMCols = cols;
      Matrix<_ElemT>::mStride = this->mOrigMRealCols;
    }


  template<typename _ElemT>
    std::ostream &
    operator << (std::ostream & rOut, Matrix<_ElemT> & rM)
    {
      for (size_t i = 0; i < rM.Rows(); i++)
      {
        for (size_t j = 0; j < rM.Cols(); j++)
        {
          rOut << rM[i][j] << " ";
        }
        rOut << std::endl;
      }
      return rOut;
    }
    
    

    
  //****************************************************************************
  //****************************************************************************
  // Constructor
  template<typename _ElemT>
    MatrixRange<_ElemT>::
    MatrixRange(const Matrix<_ElemT>& rT, 
                const size_t    ro,
                const size_t    r,
                const size_t    co,
                const size_t    c)
    {
      // point to the begining of window
      Matrix<_ElemT>::mMRows = r;
      Matrix<_ElemT>::mMCols = c;
      Matrix<_ElemT>::mStride = rT.Stride();
      Matrix<_ElemT>::mpData = rT.pData() + co + ro * rT.Stride();
    }
        
    

}// namespace STK

// #define STK_Matrix_tcc
#endif 
