#ifndef BDTree_h
#define BDTree_h

#include "common.h"
#include "ContextSample.h"
#include "BQuestion.h"
//#include <STKLib/common.h>
//#include <STKLib/ContextSample.h>

#include <vector>


namespace STK
{
  typedef std::vector<FLOAT> Data;


  /** *************************************************************************
   *  *************************************************************************
   * @brief Binary "is in set" question
   */
  class BSetQuestion : public BQuestion {
  public:
    BSetQuestion() : BQuestion(), mSet(), mPred(0)
    {}

    BSetQuestion(int pred): BQuestion(), mSet(), mPred(pred)
    {};

    BSetQuestion(const BSetQuestion& rOrig): BQuestion(rOrig), mSet(rOrig.mSet), 
    mPred(rOrig.mPred)
    {};

    virtual
    ~BSetQuestion()
    {}

    virtual bool
    Eval(const BQuestionTerm& rTerm) const;

    virtual BSetQuestion*
    Clone();

    /** 
     * @brief Dumps the question to a stream
     * 
     * @param rStream std::ostream for output
     * @param rPrefix prefix which is prepended to the output
     */
    virtual void
    Dump(std::ostream& rStream, const std::string& rPrefix) const;

    void
    DumpImplicit() const;

    // Set Question specific methods ..........................................
    int
    Predictor() const
    { return mPred; }

    void
    Set(const NGram::TokenType& rToken);

    void
    Unset(const NGram::TokenType& rToken);

    bool
    EvalRawToken(const NGram::TokenType& rToken) const;

  private:
    typedef   std::set<NGram::TokenType> SetType;
    SetType   mSet;
    int       mPred;
  };


  class Distribution
  {
  public:
    typedef FLOAT ProbType;
    typedef std::map<NGram::TokenType, FLOAT> Container;

    Distribution() 
    : mVec(), mN(0)
    {}

    Distribution(size_t n) 
    //: mVec(n), mN(0)
    :mN(0)
    {}

    Distribution(const Distribution& rOrig)
    : mVec(rOrig.mVec), mN(rOrig.mN)
    {}
    
    ~Distribution()
    {}

    /** 
     * @brief Dumps the question to a stream
     * 
     * @param rStream std::ostream for output
     * @param rPrefix prefix which is prepended to the output
     */
    void
    Dump(std::ostream& rStream, const std::string& rPrefix) const;

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

    void
    ComputeFromCounts(const Data& rData);

    void
    ComputeFromNGrams(const NGramSubset& rData);

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
   * @brief This class provides attributes to tree growing
   */
  class BDTreeBuildTraits
  {
  public:
    /// Default initialization
    BDTreeBuildTraits() :
    mMinReduction(0.001),
    mMinInData(0),
    mNPred(3),
    mVerbosity(1)
    {}

    double    mMinReduction;    ///< Minimum entropy reduction
    double    mMinInData;       ///< Minimum input data mass
    int       mNPred;           ///< Number of predictors

    int       mVerbosity;       ///< Verbosity level

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
     * @brief Builds a new tree
     * 
     * @param rData self described
     * @param rTraits BDTreeBuildTraits structure specifying building parameters
     */
    BDTree(NGramSubset& rData, BDTreeBuildTraits& rTraits, 
        const std::string& rPrefix);

    /** 
     * @brief Destructor
     */
    ~BDTree();

    bool
    IsLeaf() const
    { return (NULL == mpTree0) && (NULL == mpTree1) && (NULL == mpQuestion); }
    //..........................................................................
    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param format  format 
     * @param rTraits tree building parameters
     */
    void
    Read(std::istream& rStream, int format, BDTreeBuildTraits& rTraits);

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param format  format 
     * @param rTraits tree building parameters
     */
    void
    Write(std::ostream& rStream, int format);

    /** 
     * @brief Builds the tree from the given NGram subset
     * 
     * @param rNGrams NGram counts to use
     */
    void
    BuildFromNGrams(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits,
        const std::string& rPrefix);

    BSetQuestion*
    FindSimpleQuestion(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt) ;

    /** 
     * @brief Dumps the question to a stream
     * 
     * @param rStream std::ostream for output
     * @param rPrefix prefix which is prepended to the output
     */
    void
    Dump(std::ostream& rStream, const std::string& rPrefix) const;

  private:
    BSetQuestion*
    FindSubset_Greedy(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt,
        int pred);

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
