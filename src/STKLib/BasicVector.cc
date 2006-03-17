
#include "BasicVector.h"

namespace STK
{  
  template<>
    BasicVector<float>&
    BasicVector<float>::
    AddCVMul(const float c, const float* pV, const size_t nV)
    {
      cblas_saxpy(mLength, c, pV, 1, mpData, 1);
      return *this;
    }
}
