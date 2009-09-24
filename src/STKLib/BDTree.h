#ifndef BDTree_h
#define BDTree_h

#include "common.h"
#include "ContextSample.h"
#include "BQuestion.h"
#include "BDTree_IO.h"
#include "BasicVector.h"

#include <vector>
#include <limits>
#include <list>


namespace STK
{
  class SparceBigramMatrix;
  class BSetQuestion;
  class BDTree;


  /** *************************************************************************
   *  *************************************************************************
   * @brief Cell in sparse bigram matrix
   */
  class SparseBigramCell {
  public:
    typedef FLOAT CountType;

    SparseBigramCell() : mColumn(0), mCellCount(0)
    {};

    SparseBigramCell(int ColumnNumber, CountType CurrentCount)
      : mColumn(ColumnNumber), mCellCount(CurrentCount)
    {};

    ~SparseBigramCell()
    {};

    /** 
     * @brief Gets column index in sparse bigram matrix
     */
    int
    GetColumn() const
    { return mColumn; }

    /** 
     * @brief Gets count value in sparse bigram matrix
     */
    CountType
    GetCount() const
    { return mCellCount; }

    /** 
     * @brief Increments count of a cell in sparse bigram matrix
     * 
     * @param AddCount value count should be incremented by
     */
    CountType
    IncrementCount(CountType AddCount)
    { mCellCount += AddCount; return mCellCount; }

    /** 
     * @brief Move an element from one class to another (update both sparse matrix classes)
     * 
     * @param ColumnNumber
     * @param CurrentCount
     */
    void
    SetCell(int ColumnNumber, CountType CurrentCount)
    { mColumn = ColumnNumber; mCellCount = CurrentCount; }

  private:
    int          mColumn;
    CountType    mCellCount;
  };

  /** *************************************************************************
   *  *************************************************************************
   * @brief Sparce bigram matrix representation
   */
  class SparseBigramMatrix {
  public:
    //typedef double CountType;
    typedef FLOAT ProbType;
    typedef std::vector<NGram::ProbType> ContainerProb;
    typedef std::vector<int> ContainerInt;
    typedef SparseBigramCell* SparseBigramCellPntr;

    SparseBigramMatrix() 
    : mPresentVocabSize(0), mPredictorVocabSize(0), mSizeVector(), 
      mMarginalCountsPresent(), mMarginalCountsPredictor(), mPointerVector(), 
      mTotalSum(0)
    {}

    SparseBigramMatrix(int vocab_size)
    : mPresentVocabSize(vocab_size), mPredictorVocabSize(vocab_size), 
      mSizeVector(vocab_size), mMarginalCountsPresent(vocab_size), 
      mMarginalCountsPredictor(vocab_size), mPointerVector(vocab_size), 
      mTotalSum(0)
    {
      int i;
      for(i = 0; i < vocab_size; i++)
      {
	mSizeVector[i] = 0;
	mMarginalCountsPresent[i] = 0;
	mMarginalCountsPredictor[i] = 0;
	//mPointerVector[i] = 0;
      }
    }

    SparseBigramMatrix(int present_size, int predictor_size)
    : mPresentVocabSize(present_size), mPredictorVocabSize(predictor_size), 
      mSizeVector(predictor_size), mMarginalCountsPresent(present_size), 
      mMarginalCountsPredictor(predictor_size), mPointerVector(predictor_size), 
      mTotalSum(0)
    {
      int i;
      for(i = 0; i < predictor_size; i++)
      {
	mSizeVector[i] = 0;
	mMarginalCountsPredictor[i] = 0;
      }
      for(i = 0; i < present_size; i++)
	mMarginalCountsPresent[i] = 0;
    }

    ~SparseBigramMatrix();

    /** 
     * @brief Creates vector that tells how many different N-grams we have for 
     *        each fixed predictor position
     * 
     * @param rNGrams N-gram set we have at the node
     * @param rQuestion Current question at given node
     * @param pred Predictor position
     * @param YesAnswer 1 if we work with  left (answer=yes) subset, 0 otherwise
      */
    void
    CreateSizeVector(const NGramSubset& rNGrams, const BSetQuestion& rQuestion, 
        int pred, bool YesAnswer);

    /** 
     * @brief Allocate memory for sparce bigram matrix
     * 
     * @param rNGrams N-gram set we have at the node
     * @param rQuestion Current question at given node
     * @param YesAnswer 1 if we work with  left (answer=yes) subset, 0 otherwise
     */
    void
    AllocateMem();

    /** 
     * @brief Fills sparse bigram matrix from N-gram counts
     * 
     * @param rNGrams N-gram set we have at the node
     * @param rQuestion Current question at given node
     * @param pred Predictor position
     * @param YesAnswer 1 if we work with  left (answer=yes) subset, 0 otherwise
     */
    void
    Fill(const NGramSubset& rNGrams, const BSetQuestion& rQuestion, const int pred, const bool YesAnswer);

    /** 
     * @brief Move an element from one class to another (update both sparse matrix classes)
     * 
     * @param MoveFromThisSet the set from which we move basic element to current object
     * @param Predictor Current word in predictor we ask question about
     */
    ProbType
    InsertBasicElement(SparseBigramMatrix& MoveFromThisSet, const int Predictor);

    /** 
     * @brief Count Log-likelihood for the split
     * 
     * @param MoveFromThisSet Second set with which we count log-likelihood
     */
    ProbType
    CountLogLikelihood(const SparseBigramMatrix& MoveFromThisSet) const;

    /** 
     * @brief Count Log-likelihood for the split
     * 
     * @param MoveFromThisSet Second set with which we count log-likelihood
     */
    int
    GetSizeVectorCell(int position) const
    { return mSizeVector[position]; }

    int
    GetPresVocabSize() const
    { return mPresentVocabSize; }

    int
    GetPredVocabSize() const
    { return mPredictorVocabSize; }

    NGram::ProbType
    GetTotalCountSum() const
    { return mTotalSum; }

  private:
    int             mPresentVocabSize;
    int             mPredictorVocabSize;
    ContainerInt    mSizeVector;
    ContainerProb   mMarginalCountsPresent;
    ContainerProb   mMarginalCountsPredictor;
    std::vector<SparseBigramCellPntr> mPointerVector;
    NGram::ProbType mTotalSum;
  };

  /** *************************************************************************
   *  *************************************************************************
   * @brief Info on different sections of vocabulary containing different morphologicalunits
   */
  class MorphVocabSection
  {
  public:
    MorphVocabSection(): mFirst(0), mLast(0)
      {};

    MorphVocabSection(int frst, int lst): mFirst(frst), mLast(lst)
      {};

    MorphVocabSection(NGramSubsets& rNGrams, int pred);

    ~MorphVocabSection()
      {};

    int Size()
    { return mLast - mFirst + 1; }

    int FirstElem()
    { return mFirst; }

    int LastElem()
    { return mLast; }

  private:
    int mFirst;
    int mLast;
  };


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

    void
    RandomShuffle(const NGram::TokenType& rToken);

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
    typedef FLOAT ProbType;

    Distribution() 
    :  mN(0)
    {}

    Distribution(size_t n) 
    : mN(0)
    {}

    Distribution(const Distribution& rOrig)
    : mN(rOrig.mN)
    {}
    
    virtual
    ~Distribution()
    {}


    virtual Distribution*
    Clone() = 0;


    /** 
     * @brief Dumps the question to a stream
     * 
     * @param rStream std::ostream for output
     * @param rPrefix prefix which is prepended to the output
     */
    virtual void
    Dump(std::ostream& rStream, const std::string& rPrefix) const = 0;


    /** 
     * @brief Dumps the distribution to STDOUT
     */
    virtual void
    DumpImplicit() const;


    //..........................................................................
    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader header with parameters
     */
    virtual void
    Read(std::istream& rStream, BDTreeHeader& rHeader) = 0;

    /** 
     * @brief Reads the tree definition from a stream
     * 
     * @param rStream std::stream to read from
     * @param rHeader header with parameters
     */
    virtual void
    Write(std::ostream& rStream, BDTreeHeader& rHeader) = 0;


    /** 
     * @brief Returns the size of the distribution
     */
    virtual size_t
    Size() const = 0;
    
    /** 
     * @brief Returns soft data counts
     */
    const ProbType&
    Counts() const
    { return mN; }

    /** 
     * @brief Returns the entropy of the distribution
     */
    virtual const ProbType 
    Entropy() const = 0;

    /** 
     * @brief Computes weighted entropy of two distributions
     */
    static const ProbType
    SplitEntropy(const Distribution& rD1, const Distribution& rD2);

    /** 
     * @brief Normalizes distribution to probabilities
     */
    virtual void
    Fix() = 0;

    /** 
     * @brief Sets all records to 0
     */
    virtual void 
    Reset() = 0;

    /** 
     * @brief Merges another distribution
     * @param rDistr distribution to merge with
     */
    virtual void
    Merge(const Distribution& rDistr) = 0;



    virtual ProbType
    operator [] (const NGram::TokenType& rToken) const = 0;


    
    virtual ProbType& 
    operator [] (const NGram::TokenType& rToken) = 0;

    
    /** 
     * @brief Computes distribution from collection of grams
     */
    virtual void
    ComputeFromNGramSubsets(const NGramSubsets& rSubsets) = 0;


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
    Smooth(const Distribution& rDistr, FLOAT r) = 0;


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
     * b = N / (N+r)
     */
    virtual void
    MapAdapt(const Distribution& rDistr, FLOAT r) = 0;

    /** 
     * @brief Smoothes distribution by another distribution
     * 
     * @param rDistr Smoothing distribution
     * @param r smoothing kludge/fudge factor
     *
     * Psmoothed(s) = r*Pthis(s) + (1-r)*PrDistr(s) 
     */
    virtual void
    Interpolate(const Distribution& rDistr, FLOAT r) = 0;


    friend class BDTree;
  protected:
    ProbType           mN;   ///< Soft data count
    
  }; // class VecDistribution

  

  class MapDistribution : public Distribution
  {
  public:
    typedef FLOAT ProbType;
    typedef std::map<int, ProbType> Container;

    MapDistribution() 
    : Distribution(), mMap(), mVocabSize(0)
    {}

    MapDistribution(size_t n) 
    : Distribution(n), mMap(), mVocabSize(n)
    {}

    MapDistribution(const MapDistribution& rOrig)
    : Distribution(rOrig), mMap(rOrig.mMap), mVocabSize(rOrig.mVocabSize)
    {}

    /** 
     * @brief Loads a new distr from a stream
     * 
     * @param rStream stream to read
     * @param rHeader header containing file info
     */
    MapDistribution(std::istream& rStream, BDTreeHeader& rHeader)
    { }

    
    virtual
    ~MapDistribution()
    {}

    virtual MapDistribution*
    Clone() 
    {
      return new MapDistribution(*this);
    }

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


    virtual size_t
    Size() const 
    { return mMap.size(); }

    /** 
     * @brief Returns the entropy of the distribution
     */
    virtual const ProbType 
    Entropy() const;

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
    virtual void
    Merge(const Distribution& rDistr);

      
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
    Smooth(const Distribution& rDistr, FLOAT r);


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
     * b = N / (N+r)
     */
    virtual void
    MapAdapt(const Distribution& rDistr, FLOAT r);


    /** 
     * @brief Smoothes distribution by another distribution
     * 
     * @param rDistr Smoothing distribution
     * @param r smoothing kludge/fudge factor
     *
     * Psmoothed(s) = r*Pthis(s) + (1-r)*PrDistr(s) 
     */
    virtual void
    Interpolate(const Distribution& rDistr, FLOAT r);


    virtual ProbType
    operator [] (const NGram::TokenType& rToken) const;
    
    virtual ProbType&
    operator [] (const NGram::TokenType& rToken);


    friend class BDTree;
  private:
    Container          mMap; ///< Data distribution
    size_t             mVocabSize;
  };


  class VecDistribution : public Distribution
  {
  public:
    typedef FLOAT ProbType;
    typedef std::vector<ProbType> Container;

    VecDistribution() 
    : Distribution(), mVec()
    {}

    VecDistribution(size_t n) 
    : Distribution(n), mVec(n)
    {}

    VecDistribution(const VecDistribution& rOrig)
    : Distribution(rOrig), mVec(rOrig.mVec)
    {}

    /** 
     * @brief Loads a new distr from a stream
     * 
     * @param rStream stream to read
     * @param rHeader header containing file info
     */
    VecDistribution(std::istream& rStream, BDTreeHeader& rHeader);
    
    virtual
    ~VecDistribution()
    {}

    virtual VecDistribution*
    Clone();

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


    virtual size_t
    Size() const 
    { return mVec.size(); }

    /** 
     * @brief Returns the entropy of the distribution
     */
    virtual const ProbType 
    Entropy() const;

    // /** 
    //  * @brief Computes weighted entropy of two distributions
    //  */
    // static const ProbType
    // SplitEntropy(const VecDistribution& rD1, const VecDistribution& rD2);

    // /** 
    //  * @brief Computes distribution from collection of grams
    //  * @param first iterator to the begining of the NGram* container
    //  * @param last iterator to the endo of th NGram* container
    //  */
    // template <class InputIterator>
    //   void
    //   ComputeFromNGrams(InputIterator first, InputIterator last);

    // /** 
    //  * @brief Computes distribution from collection of grams
    //  * @param first iterator to the begining of the NGram* container
    //  * @param last iterator to the endo of th NGram* container
    //  */
    // template <class InputIterator>
    //   void
    //   ComputeFromNGramSubsets(InputIterator first, InputIterator last);

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
    virtual void
    Merge(const Distribution& rDistr);

      
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
    Smooth(const Distribution& rDistr, FLOAT r);


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
     * b = N / (N+r)
     */
    virtual void
    MapAdapt(const Distribution& rDistr, FLOAT r);


    /** 
     * @brief Smoothes distribution by another distribution
     * 
     * @param rDistr Smoothing distribution
     * @param r smoothing kludge/fudge factor
     *
     * Psmoothed(s) = r*Pthis(s) + (1-r)*PrDistr(s) 
     */
    virtual void
    Interpolate(const Distribution& rDistr, FLOAT r);



    virtual ProbType
    operator [] (const NGram::TokenType& rToken) const;

    virtual ProbType& 
    operator [] (const NGram::TokenType& rToken)
    { return mVec[rToken]; }

    
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
    mMMIAlpha(1.0), 
    mMMIEta(1.0), 
    mMinMMIIncrease(0.0),
    mMapAdapt(false),
    mLVCSR(false),
    mRandomizeTree(false),
    mRandomizePredictors(false),
    mRandomizeQuestions(false),
    mMorphologicalPredictors(0),
    mUseMapDistribution(false)
    {}

    FLOAT     mMinReduction;    ///< Minimum entropy reduction
    FLOAT     mMinInData;       ///< Minimum input data mass
    int       mOrder;           ///< Number of predictors
    int       mMaxDepth;        ///< Maximum tree depth (root is depth 0)
    int       mCurDepth;        ///< Maximum tree depth (root is depth 0)

    FLOAT     mSmoothR;         ///< Bottom-up recursive smoothing r-factor
    FLOAT     mAdaptR;          ///< Adapt smoothing r-factor
    int       mVerbosity;       ///< Verbosity level

    bool      mMMItoggle;
    FLOAT     mMMIAlpha;        ///< MMI alpha constant
    FLOAT     mMMIEta;          ///< MMI eta constant
    FLOAT     mMinMMIIncrease;  ///< 

    bool      mMapAdapt;
    
    bool      mLVCSR;
    bool      mRandomizeTree;
    bool      mRandomizePredictors;
    bool      mRandomizeQuestions;
    int       mMorphologicalPredictors;
    
    bool      mUseMapDistribution;
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
     * @brief Builds a new tree using heldout data as stopping criterium
     * 
     * @param rData        self described
     * @param rHeldoutData NGram heldout data
     * @param rTraits      BDTreeBuildTraits structure specifying building parameters
     */
    BDTree(NGramSubsets& rData, NGramSubsets* rHeldoutData, BDTreeBuildTraits& rTraits, BDTree* pParent, 
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
    BuildFromNGrams(NGramSubsets& rNGrams, NGramSubsets* rHeldoutNGrams,  BDTreeBuildTraits& rTraits, 
        BDTree* pParent, const std::string& rPrefix);

    /** 
     * @brief Builds the tree from the given NGram subset using heldout data as stopping criterium
     * 
     * @param rNGrams NGram counts to use
     * @param rHeldoutData  heldout N-grams
     */
    void 
    BuildFromNGrams(NGramSubsets& rNGrams,  NGramPool& rHeldoutData, BDTreeBuildTraits& rTraits, 
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
    FindSimpleQuestion(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt) ;


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
    FindSimpleQuestion_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitMMI) ;


    /** 
     * @brief Scores a single n-gram
     * @param rNGram data to score
     * @return Score based on the decision tree structure and precomputed 
     * smoothed probabilities
     *
     * NOTE: Smoothing is not done on-line. It has to be precomputed using 
     * ComputeBackoffDists
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
     * @brief For given set of N-grams scores the probabilities and outputs in ARPA-like format (though unnormalized in conventional way since only needed subset of all trigrams is scored)
     * 
     * @param rNGrams N-gram data
     * @param lm_stream output stream
     * @param rOrder N-gram order
     * 
     */
    void 
    OutputNGramSubsetARPA(const NGramSubset& rNGrams, std::ostream& lm_stream, int rOrder);

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
     * @brief Recomputes distribution in non-leaf nodes
     */
    void
    RecomputeDists();

    /** 
     * @brief Fills leaf distributions with clustered data counts
     * 
     * @param rNGrams the data
     *
     * NOTE: Use FillLeafSupervector function to stack the leaves to form one
     * super vector
     */
    void
    CollectCounts(const NGramSubset& rNGrams);

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
    AdaptLeafDists(const NGramSubset& rNGrams, BDTreeBuildTraits& rTraits);

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

    void
    DumpImplicit() const;


    /** 
     * @brief Returns this tree's leaves in a vector container
     * 
     * @param rCollection 
     * 
     * @return number of leaves
     */
    int
    GetLeaves(std::vector<BDTree*>& rCollection);

    void
    FillLeafSupervector(BasicVector<FLOAT>& rVector, bool backoff, 
        bool includeCounts);

    void
    PushLeafSupervector(const std::vector<FLOAT>& rVector, bool backoff, 
        bool includeCounts);

    const Distribution*
    cpDistribution()
    { return mpDist; }

  private:
    BSetQuestion*
    FindSubset_Greedy(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitEnt,
        int pred);

    BSetQuestion*
    FindSubset_Greedy_MMI(NGramSubsets& rNGrams, BDTreeBuildTraits& rTraits, FLOAT* pSplitMMI,
        int pred);


  private:
    BDTree*           mpTree0;
    BDTree*           mpTree1;
    BQuestion*        mpQuestion;
    Distribution*     mpDist;
    Distribution*     mpBackoffDist;
  }; // class BDTree
} // namespace STK

#include "BDTree.tcc"

#endif
// EOF
