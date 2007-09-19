#ifndef BDTree_h
#define BDTree_h

#include "common.h"
#include "ContextSample.h"
#include "BQuestion.h"
#include "BDTree_IO.h"

#include <vector>
#include <limits>


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


  /** 
   * @brief Generalized distribution representation
   */
  class Distribution
  {
  public:
    typedef double ProbType;

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
//    Smooth(const Distribution& rDistr, double r);
//
//    ProbType 
//    operator [] (const NGram::TokenType& rToken);

  protected:
    ProbType           mN;   ///< Soft data count
    
  }; // class VecDistribution


  class VecDistribution : public Distribution
  {
  public:
    typedef double ProbType;
    typedef std::vector<double> Container;

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

    /** 
     * @brief Computes distribution from collection of grams
     * @param first iterator to the begining of the NGram* container
     * @param last iterator to the endo of th NGram* container
     */
    template <class InputIterator>
      void
      ComputeFromNGrams(InputIterator first, InputIterator last);

    /** 
     * @brief Computes distribution from collection of grams
     * @param first iterator to the begining of the NGram* container
     * @param last iterator to the endo of th NGram* container
     */
    template <class InputIterator>
      void
      ComputeFromNGramSubsets(InputIterator first, InputIterator last);

    void
    ComputeFromNGramSubsets(const NGramSubsets& rSubsets);

    /** 
     * @brief Normalizes distribution to probabilities
     */
    void
    Fix();

    /** 
     * @brief Sets all records to 0
     */
    void 
    Reset();

    /** 
     * @brief Merges another distribution
     * @param rDistr distribution to merge with
     */
    void
    Merge(const VecDistribution& rDistr);

      
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
    Smooth(const VecDistribution& rDistr, double r);


    ProbType
    operator [] (const NGram::TokenType& rToken);

    friend class BDTree;

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
    mOrder(3),
    mMaxDepth(std::numeric_limits<int>::max()),
    mCurDepth(0),
    mSmoothR(0.0),
    mAdaptR(0.0),
    mVerbosity(1),
    mMMItoggle(false),
    mMultiScript(false)
    {}

    double    mMinReduction;    ///< Minimum entropy reduction
    double    mMinInData;       ///< Minimum input data mass
    int       mOrder;           ///< Number of predictors
    int       mMaxDepth;        ///< Maximum tree depth (root is depth 0)
    int       mCurDepth;        ///< Maximum tree depth (root is depth 0)

    double    mSmoothR;         ///< Bottom-up recursive smoothing r-factor
    double    mAdaptR;          ///< Adapt smoothing r-factor
    int       mVerbosity;       ///< Verbosity level

    bool      mMMItoggle;
    double    mMinMMIIncrease;

    bool      mMultiScript;
  }; //class BDTreeAttributes


  /** 
   * @brief Holds some statistics about the tree
   */
  struct BDTreeInfo
  {
    BDTreeInfo()
    : mAverageEntropy(0.0), mTotalData(0.0), mTotalNodes(0), mTotalLeaves(0)
    {}

    ~BDTreeInfo()
    {}

    double  mAverageEntropy;
    double  mTotalData;
    size_t  mTotalNodes;
    size_t  mTotalLeaves;
  };


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
    BDTree(NGramSubsets& rData, BDTreeBuildTraits& rTraits, BDTree* pParent, 
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
    void
    LoadFile(const std::string& rName);

    void
    SaveFile(const std::string& rName);


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
    BuildFromNGrams(NGramSubsets& rNGrams,  BDTreeBuildTraits& rTraits, 
        BDTree* pParent, const std::string& rPrefix);


    void
    GetInfo(BDTreeInfo& rInfo, const BDTree* pParent) const;

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
    FindSimpleQuestion(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitEnt) ;


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
    FindSimpleQuestion_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitMMI) ;


    /** 
     * @brief Scores a single n-gram
     * @param rNGram data to score
     * @return Score based on the decision tree structure and precomputed smoothed probabilities
     *
     * NOTE: Smoothing is not done on-line. It has to be precomputed using ComputeBackoffDists
     */
    double 
    ScoreNGram(const NGram& rNGram);

    /** 
     * @brief Scores data
     * 
     * @param rNGrams utterance data
     * @param r smoothing factor
     * 
     * @return utterance score
     */
    double
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
    ComputeBackoffDists(BDTree* pParent, double r);

    /** 
     * @brief Recomputes distribution in non-leaf nodes
     */
    void
    RecomputeDists();

    /** 
     * @brief Adapts leaf distributions using new data and smoothing factor
     * 
     * @param rNGrams adaptation data
     * @param r 
     *
     * Pnew(s) = b*Pthis(s) + (1-b)*Porig(s) 
     *   where
     * Pthis = #(s) / #(all_symbols_in_leaf), 
     *   and
     * b = #(s) / (#(s)+r)
     */
    void
    AdaptLeafDists(const NGramSubset& rNGrams, double r);

    /** 
     * @brief Flatens distributions in all leaves
     */
    void
    ResetLeaves();

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
    FindSubset_Greedy(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitEnt,
        int pred);

    BSetQuestion*
    FindSubset_Greedy_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, double* pSplitMMI,
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
