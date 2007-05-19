#include "BDTree.h"

namespace STK
{
  //***************************************************************************/
  //***************************************************************************/
  const Distribution::ProbType
  Distribution::
  Entropy() const
  {
    double e = 0;
    vector<ProbType>::const_iterator i;

    for (i=mVec.begin(); i!=mVec.end(); ++i) {
      e -= *i * log(*i);
    }
    
    return e;
  }


  //***************************************************************************/
  //***************************************************************************/
  const Distribution::ProbType
  Distribution::
  SplitEntropy(const Distribution& rD1, const Distribution& rD2) const
  {
    double w1 = rD1 / (rD1.mN + rD2.mN);
    double w2 = rD2 / (rD1.mN + rD2.mN);

    return w1*rD1.Entropy() + w2*rD2.Entropy();
  }


  //***************************************************************************/
  //***************************************************************************/
  void
  Distribution::
  ComputeFromCounts(const Data& rData)
  {
    Data::const_iterator  i;
    vector<ProbType>::const_iterator  j;
    double                n = 0;
      
    // set probs to 0
    for (j=mVec.begin(); j!=mVec.end(); ++j) {
      *j = 0.0;
    }

    for (i=rData.begin(); i!=rData.end(); ++i) {
      n += *i;
    }

    for (j=mVec.begin(); j!=mVec.end(); ++j) {
      *i /= n;
    }
  }


} // namespace STK

// EOF
