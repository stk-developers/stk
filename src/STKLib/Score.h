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

#include "common.h"
#include "Matrix.h"

namespace STK
{
  /**
   *  @brief Representation of gradient supervector
   */
  class GradientSV
  {
  public:
    std::string   mFileName;      ///< file name for storage
    Matrix<FLOAT> mMatrix;        ///< the matrix itselft
  };
  
  
  void
  /**
   * @brief Accumulates a gradient supervector
   * @param rMatrix matrix that is to be constructed
   * @param pFeaVector current feature vector
   * @param rHmm HMM (or better GMM) to be analyzed
   *
   * The HMM should actually be single state GMM or else assertion fails
   */
  gradient_supervector_accum( Matrix<FLOAT>& rMatrixT, 
                              FLOAT* pFeaVector,
                              const Hmm& rHmm);

}
