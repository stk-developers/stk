#include "Models.h"
#include "Decoder.h"

namespace STK
{
  //****************************************************************************
  //****************************************************************************
  void
  gradient_supervector_accum(Matrix<FLOAT>& rMatrixT, FLOAT* pFeaVector,
                             const Hmm& rHmm)
  {
    FLOAT           p;
    FLOAT           posterior_denom;
    State*          state;
    Mixture*        mixture;
    const FLOAT*    mean;
    const FLOAT*    variance;
    
    // it has to be provided that we have GMM
    assert(rHmm.NStates() == 3);
    
    state = rHmm.mpState[0];
    assert(KID_DiagC == state->mOutPdfKind);
    
    posterior_denom = DiagCGaussianMixtureDensity(state, pFeaVector, NULL);
    
    // go through the model components (mixtures)
    for (size_t j = 0; j < rMatrixT.Rows(); j++)
    {
      mixture   = state->mpMixture[j].mpEstimates;
      mean      = mixture->mpMean->mVector.cpData();
      variance  = mixture->mpVariance->mVector.cpData();
      
      // the posterior
      p = _EXP(DiagCGaussianDensity(mixture, pFeaVector, NULL)
               + state->mpMixture[j].mWeight
               - posterior_denom);
              
      // go through the feature components
      for (size_t i = 0; i < rMatrixT.Cols(); i++)
      {
        rMatrixT[j][i] += p * (pFeaVector[i] - mean[i]) * _SQRT(variance[i]);
      }
    } // go through the feature components
  } // gradient_supervector(...)
  
} // namespace STK
