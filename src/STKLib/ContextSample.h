#ifndef ContextSample_h
#define ContextSample_h

#include "common.h"
#include "BQuestion.h"

#include <iostream>
#include <map>
#include <string>
#include <set>


namespace STK
{
  const int LINE_BUF_DEFAULT_MIN_LENGTH = 2048;
  const int CONTEXT_SAMPLE_DEFAULT_MIN_LENGTH = 256;

  class NGram;
  class NGramPool;
  class NGramSubset;
  class VocabularyTable;

  /** 
   * @brief 
   */
  class NGram : public BQuestionTerm
  {
  public:
    typedef FLOAT ProbType;
    typedef int   TokenType;

    //..........................................................................
    NGram()
    {}

    ~NGram()
    {}

    //..........................................................................
    /// Sample weight, soft counts, whatever
    const ProbType&
    Counts() const
    { return mCounts; }

    const TokenType&
    operator[](const size_t index) const
    { return mpTokens[index]; }

    TokenType&
    operator[](const size_t index)
    { return mpTokens[index]; }

    friend std::ostream&
    operator<<(std::ostream& rOstr, const NGram& rCS);
    
    friend class NGramPool;

  private:
    ProbType      mCounts;    ///< Sample's weight
    TokenType*    mpTokens;   ///< Pointer to token array 
  }; // class NGram


  /** *************************************************************************
   *  *************************************************************************
   * @brief Maps real token names to their integer value
   */
  class VocabularyTable 
  : public std::map<std::string, int>
  {
  public:
    typedef std::map<std::string,int> StringMapType;
    typedef std::vector<std::string>  IntMapType;

    VocabularyTable()
    : mStrMap(), mIntMap()
    {}

    /// int to string conversion
    const std::string&
    IToA(const int) const;

    /// string to int conversion
    int
    AToI(const std::string& rString);
    
    /// loads the vocab from a file
    VocabularyTable&
    LoadFromFile(const std::string& rFName);

    /// loads the vocab from a stream
    VocabularyTable&
    LoadFromStream(std::istream& rIStream);

  private:
    StringMapType   mStrMap;
    IntMapType      mIntMap;
    bool            mStrictMode; ///< if true, throw errors
  }; // class VocabularyTable

  
  /** *************************************************************************
   *  *************************************************************************
   * @brief NGram storage.
   * 
   * This structure just holds pointers to the NGram pool. No memory allocation 
   * is done here
   */
  class NGramSubset
  {
  public:
    typedef std::vector<NGram*> NGramContainer;

    NGramSubset()
    : mData(), mOrder(0), mpPool(NULL)
    {}

    NGramSubset(size_t order)
    : mData(), mOrder(order), mpPool(NULL)
    {}

    NGramSubset(const NGramSubset& rOrig)
    : mData(rOrig.mData), mOrder(rOrig.mOrder), mpPool(rOrig.mpPool)
    {}

    virtual
    ~NGramSubset()
    {}

    /// Returns context size
    const size_t
    Order() const
    { return mOrder; }

    /// gives access to the parrent pool
    const NGramPool&
    Parent() const
    { return *mpPool;}

    /// computes data entropy
    FLOAT
    Entropy() const;

    /// computes split data entropy, provided question rQuestion
    FLOAT
    SplitEntropy(const BQuestion& rQuestion) const;

    /// returns total data mass
    FLOAT
    Mass() const;

    /** 
     * @brief Splits data according to a given question
     * 
     * @param rQuestion question that decides about the split
     * @param rData0 negative answer data
     * @param rData1 positive answer data
     */
    void
    Split(const BQuestion& rQuestion, NGramSubset& rData0, NGramSubset& rData1);

    friend class Distribution;
    friend class BDTree;

  protected:
    NGramContainer       mData;
    size_t               mOrder;
    NGramPool*           mpPool;
  };


  /** *************************************************************************
   *  *************************************************************************
   * @brief Parent storage for NGrams
   */
  class NGramPool : public NGramSubset
  {
  private:
    NGram*
    pNewNGram();

  public:
    typedef std::vector<NGram*> NGramContainer;

    NGramPool()
    : NGramSubset(), 
      mpPredictorTable(NULL), mpTargetTable(NULL)
    { mpPool = this; }

    NGramPool(size_t n)
    : NGramSubset(n),
      mpPredictorTable(NULL), mpTargetTable(NULL)
    { mpPool= this; }

    virtual
    ~NGramPool();

    void 
    AddFromFile(const std::string& rFileName, const double& weight);

    void 
    AddFromStream(std::istream& rStream, const double& weight);

    void
    setPredictorVocab(VocabularyTable* pVocab)
    { mpPredictorTable = pVocab;} 

    void
    setTargetVocab(VocabularyTable* pVocab)
    { mpTargetTable = pVocab;} 


    /** 
     * @brief Search for a given NGram
     * @param pNGram pointer to the NGram to look for
     * @return Pointer to the found NGram or NULL if not found
     */
    NGram*
    FindNGram(const NGram* pNGram);

    /** 
     * @brief Compares two ngram structures
     * @return true if different, false if equal
     */
    const bool
    CompareNGram(const NGram* pFirst, const NGram* pSecond) const;
    
    /// Returns number of samples in the pool
    const size_t
    Size() const
    { return mData.size(); }

    /// Dumps the whole pool to the stream
    friend std::ostream&
    operator <<(std::ostream& rOstr, const NGramPool& rWhat);

  private:
    VocabularyTable*     mpPredictorTable;
    VocabularyTable*     mpTargetTable;
  }; // class NGramPool
  
  
} // namespace STK

//#include <STKLib/ContextSample.tcc>

#endif // ContextSample_h
// EOF
