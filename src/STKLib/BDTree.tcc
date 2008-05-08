#ifndef BDTree_tcc
#define BDTree_tcc

namespace STK
{
  //***************************************************************************/
  //***************************************************************************/
//  template <class InputIterator>
//    void
//    VecDistribution::
//    ComputeFromCounts(InputIterator first)
//    {
//      InputIterator               i;
//      Container::iterator         j;
//      double                      n = 0.0;
//        
//      // accumulate each dist and sum overal data mass
//      for (j=mVec.begin(), i=first; j!=mVec.end(); ++i, ++j) {
//        *j   = *i;
//        n   += *i;
//      }
//
//      // normalize
//      for (j=mVec.begin(); j!=mVec.end(); ++j) {
//        *j /= n;
//      }
//    }


//  //***************************************************************************/
//  //***************************************************************************/
//  template <class InputIterator>
//    void
//    VecDistribution::
//    ComputeFromNGrams(InputIterator first, InputIterator last)
//    {
//      double              mass = 0.0;
//      InputIterator       i;
//      Container::iterator j;
//        
//      // set probs to 0
//      for (j=mVec.begin(); j!=mVec.end(); ++j) {
//        *j = 0.0;
//      }
//
//      // accumulate each dist and sum overal data mass
//      for (i=first; i!=last; ++i) {
//        NGram::TokenType* p_token = &((**i)[0]);
//        NGram::ProbType   counts = (*i)->Counts();
//
//        mVec[*p_token] += counts;
//        mass           += counts;
//      }
//
//      // normalize
//      for (j=mVec.begin(); j!=mVec.end(); ++j) {
//        *j = *j / mass;
//      }
//
//      mN = mass;
//    }
//
//  //***************************************************************************/
//  //***************************************************************************/
//  template <class InputIterator>
//    void
//    VecDistribution::
//    ComputeFromNGramSubsets(InputIterator first, InputIterator last)
//    {
//      double              mass = 0.0;
//      InputIterator       i;
//      Container::iterator j;
//        
//      // set probs to 0
//      for (j=mVec.begin(); j!=mVec.end(); ++j) {
//        *j = 0.0;
//      }
//
//      // accumulate each dist and sum overal data mass
//      for (i=first; i!=last; ++i) {
//        NGram::TokenType* p_token = &((**i)[0]);
//        NGram::ProbType   counts = (*i)->Counts();
//
//        mVec[*p_token] += counts;
//        mass           += counts;
//      }
//
//      // normalize
//      for (j=mVec.begin(); j!=mVec.end(); ++j) {
//        *j = *j / mass;
//      }
//
//      mN = mass;
//    }
} // namespace STK

#endif // BDTree_tcc
