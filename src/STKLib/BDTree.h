#ifndef BDTree_h
#define BDTree_h

#include "common.h"
#include "ContextSample.h"
#include "BQuestion.h"
#include "BDTree_IO.h"

#include <vector>


namespace STK
{
  class BSetQuestion;


  /** *************************************************************************
   *  *************************************************************************
   * @brief Binary "is in set" question
   */
  class BSetQuestion : public BQuestion {
  public:
    BSetQuestion() : BQuestion(), mSet(), mPred(0)
    {}

    BSetQuestion(int pred, int vocab_size)
    : BQuestion(), mSet(vocab_size), mPred(pred)
    {};

    BSetQuestion(const BSetQuestion& rOrig)
    : BQuestion(rOrig), mSet(rOrig.mSet), mPred(rOrig.mPred)
    {};

    BSetQuestion(std::istream& rStream, BDTreeHeader& rHeader);

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

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader header with parameters
     */
    virtual void
    Read(std::istream& rStream, BDTreeHeader& rHeader);

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader header with parameters
     */
    virtual void
    Write(std::ostream& rStream, BDTreeHeader& rHeader);


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
    //typedef   std::set<NGram::TokenType> SetType;
    typedef   std::vector<bool> SetType;
    SetType   mSet;
    int       mPred;
  };


  class Distribution
  {
  public:
    typedef FLOAT ProbType;

    Distribution() 
    :  mN(0)
    {}

    Distribution(size_t n) 
    :  mN(0)
    {}

    Distribution(const Distribution& rOrig)
    : mN(rOrig.mN)
    {}
    
    virtual
    ~Distribution()
    {}

//    /** 
//     * @brief Dumps the question to a stream
//     * 
//     * @param rStream std::ostream for output
//     * @param rPrefix prefix which is prepended to the output
//     */
//    virtual void
//    Dump(std::ostream& rStream, const std::string& rPrefix) const;

    /** 
     * @brief Returns soft data counts
     */
    const ProbType&
    Counts() const
    { return mN; }
//
//    /** 
//     * @brief Returns the entropy of the distribution
//     */
//    virtual const ProbType 
//    Entropy() const;
//
//    /** 
//     * @brief Computes weighted entropy of two distributions
//     */
//    static const ProbType
//    SplitEntropy(const Distribution& rD1, const Distribution& rD2);
//
//
//    virtual void
//    ComputeFromNGrams(const NGramSubset& rData);
//
//    template <class InputIterator>
//      void
//      ComputeFromCounts(InputIterator first);
//
//
//    template <class InputIterator>
//      void
//      ComputeFromSamples(InputIterator first, InputIterator last);
//
//
//    /** 
//     * @brief Smoothes distribution by another distribution
//     * 
//     * @param rDistr Smoothing distribution
//     * @param r smoothing kludge/fudge factor
//     *
//     * Psmoothed(s) = b*Pthis(s) + (1-b)*PrDistr(s) 
//     *
//     *   where
//     *
//     * b = #(s) / (#(s)+r)
//     */
//    virtual void
//    Smooth(const Distribution& rDistr, FLOAT r);
//
//    ProbType 
//    operator [] (const NGram::TokenType& rToken);

  protected:
    ProbType           mN;   ///< Soft data count
    
  }; // class VecDistribution


  class VecDistribution : public Distribution
  {
  public:
    typedef FLOAT ProbType;
    typedef std::vector<FLOAT> Container;

    VecDistribution() 
    : Distribution(), mVec()
    {}

    VecDistribution(size_t n) 
    : Distribution(n), mVec(n)
    {}

    VecDistribution(const VecDistribution& rOrig)
    : Distribution(rOrig), mVec(rOrig.mVec)
    {}
    
    ~VecDistribution()
    {}

    /** 
     * @brief Dumps the question to a stream
     * 
     * @param rStream std::ostream for output
     * @param rPrefix prefix which is prepended to the output
     */
    virtual void
    Dump(std::ostream& rStream, const std::string& rPrefix) const;

    //..........................................................................
    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader header with parameters
     */
    void
    Read(std::istream& rStream, BDTreeHeader& rHeader);

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader header with parameters
     */
    void
    Write(std::ostream& rStream, BDTreeHeader& rHeader);


    /** 
     * @brief Returns the entropy of the distribution
     */
    virtual const ProbType 
    Entropy() const;

    /** 
     * @brief Computes weighted entropy of two distributions
     */
    static const ProbType
    SplitEntropy(const VecDistribution& rD1, const VecDistribution& rD2);


    template <class InputIterator>
      void
      ComputeFromNGrams(InputIterator first, InputIterator last);


    /** 
     * @brief Smoothes distribution by another distribution
     * 
     * @param rDistr Smoothing distribution
     * @param r smoothing kludge/fudge factor
     *
     * Psmoothed(s) = b*Pthis(s) + (1-b)*PrDistr(s) 
     *
     *   where
     *
     * b = #(s) / (#(s)+r)
     */
    virtual void
    Smooth(const VecDistribution& rDistr, FLOAT r);


    ProbType 
    operator [] (const NGram::TokenType& rToken);

  private:
    Container          mVec; ///< Data distribution vector
    
  }; // class VecDistribution



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
    mSmoothR(0.0),
    mVerbosity(1)
    {}

    double    mMinReduction;    ///< Minimum entropy reduction
    double    mMinInData;       ///< Minimum input data mass
    int       mNPred;           ///< Number of predictors

    double    mSmoothR;         ///< Bottom-up recursive smoothing r-factor
    int       mVerbosity;       ///< Verbosity level

  }; //class BDTreeAttributes


  /** 
   * @brief Binary decision tree
   */
  class BDTree
  {
  private:
    /** 
     * @brief Plain constructor
     */
    BDTree();

  public:
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
     * @brief Loads a new tree from a stream
     * 
     * @param rStream stream to read
     * @param rHeader header containing file info
     */
     BDTree(std::istream& rStream, BDTreeHeader& rHeader);


    /** 
     * @brief Destructor
     */
    ~BDTree();


    BDTree*
    Clone()
    { return new BDTree(*this); }


    bool
    IsLeaf() const
    { return (NULL == mpTree0) && (NULL == mpTree1) && (NULL == mpQuestion); }


    //..........................................................................
    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader tree building parameters
     */
    void
    Read(std::istream& rStream, BDTreeHeader& rHeader);

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader tree building parameters
     */
    void
    Write(std::ostream& rStream, BDTreeHeader& rHeader);


    /** 
     * @brief Builds the tree from the given NGram subset
     * 
     * @param rNGrams NGram counts to use
     */
    void
    BuildFromNGrams(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits,
        const std::string& rPrefix);


    /** 
     * @brief Finds "optimal" split question for this node
     * 
     * @param rNGrams data set
     * @param rTraits tree building parameters
     * @param pSplitEnt fill with split entropy
     * 
     * @return pointer to new question if succeede, otherwise NULL
     */
    BSetQuestion*
    FindSimpleQuestion(NGramSubset& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt) ;


    /** 
     * @brief Scores a single n-gram
     * @param rNGram data to score
     * @return Score based on the decision tree structure and precomputed smoothed probabilities
     *
     * NOTE: Smoothing is not done on-line. It has to be precomputed using ComputeBackoffDists
     */
    FLOAT 
    ScoreNGram(const NGram& rNGram);

    /** 
     * @brief Scores data
     * 
     * @param rNGrams utterance data
     * @param r smoothing factor
     * 
     * @return utterance score
     */
    FLOAT
    ScoreNGramSubset(const NGramSubset& rNGrams);

    /** 
     * @brief Computes smoothed probabilities using smoothing factor r
     * 
     * @param pParent parent node, NULL for root
     * @param r smoothing factor
     *
     * See Distribution::Smooth(...) for details
     */
    void
    ComputeBackoffDists(BDTree* pParent, FLOAT r);

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
    BDTree*           mpTree0;
    BDTree*           mpTree1;
    BQuestion*        mpQuestion;
    VecDistribution   mDist;
    VecDistribution*  mpBackoffDist;
  }; // class BDTree
} // namespace STK

#include "BDTree.tcc"

#endif
// EOF
