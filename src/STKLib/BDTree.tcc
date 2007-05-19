#ifndef BDTree_tcc
#define BDTree_tcc

namespace STK
{
  //***************************************************************************/
  //***************************************************************************/
  template <class InputIterator>
    void
    Distribution::
    ComputeFromCounts(InputIterator first)
    {
      InputIterator               i;
      Container::iterator         j;
      double                      n = 0.0;
        
      // accumulate each dist and sum overal data mass
      for (j=mVec.begin(), i=first; j!=mVec.end(); ++i, ++j) {
        *j   = *i;
        n   += *i;
      }

      // normalize
      for (j=mVec.begin(); j!=mVec.end(); ++j) {
        *j /= n;
      }
    }


  //***************************************************************************/
  //***************************************************************************/
  template <class InputIterator>
    void
    Distribution::
    ComputeFromSamples(InputIterator first, InputIterator last)
    {
      InputIterator i;
      Container::const_iterator  j;
      double                n = 0.0;
        
      // set probs to 0
      for (j=mVec.begin(); j!=mVec.end(); ++j) {
        *j = 0.0;
      }

      // accumulate each dist and sum overal data mass
      for (i=first; i!=last; ++i) {
        mVec[i->rPtee()]  += i->Weight();
        n                 += i->Weight();
      }

      // normalize
      for (i=first; i!=last; ++i) {
        mVec[i->rPtee()] /= n;
      }
    }
} // namespace STK

#endif // BDTree_tcc
