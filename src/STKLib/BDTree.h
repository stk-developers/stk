#ifndef BDTree_h
#define BDTree_h

#include <STKLib/common.h>
#include <STKLib/ContextSample.h>

#include <vector>


namespace STK
{
  typedef std::vector<FLOAT> Data;

  class Distribution
  {
  public:
    typedef FLOAT ProbType;
    typedef std::vector<FLOAT> Container;

    Distribution() 
    : mVec(), mN(0)
    {}

    Distribution(size_t n) 
    : mVec(n), mN(0)
    {}

    Distribution(const Distribution& rOrig)
    : mVec(rOrig.mVec), mN(rOrig.mN)
    {}
    
    ~Distribution()
    {}

    /** 
     * @brief Returns soft data counts
     */
    const ProbType&
    Counts() const
    { return mN; }

    /** 
     * @brief Returns the entropy of the distribution
     */
    const ProbType 
    Entropy() const;

    /** 
     * @brief Computes weighted entropy of two distributions
     */
    static const ProbType
    SplitEntropy(const Distribution& rD1, const Distribution& rD2);


    void
    BuildFromCounts(const Data& rData);

    template <class InputIterator>
      void
      ComputeFromCounts(InputIterator first);


    template <class InputIterator>
      void
      ComputeFromSamples(InputIterator first, InputIterator last);


  private:
    Container          mVec; ///< Data distribution vector
    ProbType           mN;   ///< Soft data count
    
  }; // class Distribution


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
    Eval(const ContextSample& rTerm) const = 0;

    /** 
     * @brief Creates an exact copy of this
     * @return Pointer to new instance of this
     */
    virtual BQuestion
    clone() = 0;
  }; // class Question


  /** 
   * @brief This class provides attributes to tree growing
   */
  class BDTreeAttributes 
  {
  public:


  }; //class BDTreeAttributes


  /** 
   * @brief Binary decision tree
   */
  class BDTree
  {
  public:
    
    /** 
     * @brief Plain constructor
     */
    BDTree();

    /** 
     * @brief Copy constructor
     */
    BDTree(const BDTree& rOrig);

    /** 
     * @brief Destructor
     */
    ~BDTree();

    //..........................................................................
    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param format  format 
     * @param rTraits tree building parameters
     */
    void
    Read(std::istream& rStream, int format, BDTreeAttributes& rTraits);

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param format  format 
     * @param rTraits tree building parameters
     */
    void
    Write(std::ostream& rStream, int format);

  private:
    BDTree*         mpTree0;
    BDTree*         mpTree1;
    BQuestion*      mpQuestion;
    Distribution    mDist;
  
  }; // class BDTree
} // namespace STK

#include "BDTree.tcc"

#endif
// EOF
