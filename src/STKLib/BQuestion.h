#ifndef BQuestion_h
#define BQuestion_h

#include "BDTree_IO.h"
//#include "BDTree.h"

namespace STK {
  class BQuestionTerm {
  };

  /** 
   * @brief General binary question
   */
  class BQuestion
  {
  public:
    // constructors and destructor .............................................
    BQuestion()
    {}

    virtual
    ~BQuestion()
    {}

    // .........................................................................
    /** 
     * @brief Evaluates rTerm and returns true or false
     * @param rTerm Term (Sample) to be evaluated
     * @return true or false
     */
    virtual bool
    Eval(const BQuestionTerm& rTerm) const = 0;

    /** 
     * @brief Creates an exact copy of this
     * @return Pointer to new instance of this
     */
    virtual BQuestion*
    Clone() = 0;

    /** 
     * @brief Dumps the question to a stream
     * 
     * @param rStream std::ostream for output
     * @param rPrefix prefix which is prepended to the output
     */
    virtual void
    Dump(std::ostream& rStream, const std::string& rPrefix) const = 0;

    virtual void
    Write(std::ostream& rStream, BDTreeHeader& rHeader) = 0;

    virtual void
    Read(std::istream& rStream, BDTreeHeader& rHeader) = 0;

  }; // class Question


}; // namespace STK;

//#define BQuestion_h
#endif
