
#include "BQuestion.h"
#include "ContextSample.h"

namespace STK {
  // virtual
  bool
  BSetQuestion::
  Eval(const BQuestionTerm& rTerm) const
  {
    NGram::TokenType* p_token = static_cast<NGram::TokenType*>(&rTerm);

    return (mSet.end() != mSet.find(*p_token));
  }
}; // namespace STK
