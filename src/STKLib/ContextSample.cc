#include <STKLib/ContextSample.h>

namespace STK {
  //***************************************************************************/
  //***************************************************************************/
  std::ostream&
  operator << (std::ostream& rOstr, const ConstSamplePool& rWhat) 
  {
    std::vector<ContextSampe>::iterator i;

    for (i=mData.begin(); i!=mData.end(); ++i) {
      rOstr << *i;
    }
  }
} // namespace STK

// EOF
