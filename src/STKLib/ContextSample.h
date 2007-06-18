#ifndef ContextSample_h
#define ContextSample_h

#include "common.h"

#include <iostream>
#include <map>
#include <string>


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
  class NGram
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


  /** 
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

  
  /** 
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
    {}

    NGramSubset(const NGramSubset& rOrig)
    : mData(rOrig.mData), mpPool(rOrig.mpPool)
    {}

    ~NGramSubset()
    {}

    /// gives access to the parrent pool
    const NGramPool&
    Parrent() const
    { return *mpPool;}

    /// returns total data mass
    FLOAT
    Mass() const;

  protected:
    NGramContainer       mData;
    NGramPool*           mpPool;
  };


  /** 
   * @brief Parent storage for NGrams
   */
  class NGramPool
  {
  private:
    NGram*
    pNewNGram();

  public:
    typedef std::vector<NGram*> NGramContainer;

    NGramPool()
    : mData(), mOrder(0), mpPredictorTable(NULL), mpTargetTable(NULL)
    {}

    NGramPool(size_t n)
    : mData(), mOrder(n), mpPredictorTable(NULL), mpTargetTable(NULL)
    {}

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


    /// Returns context size
    const size_t
    Order() const
    { return mOrder; }

    /// Returns number of samples in the pool
    const size_t
    Size() const
    { return mData.size(); }

    /// Dumps the whole pool to the stream
    friend std::ostream&
    operator <<(std::ostream& rOstr, const NGramPool& rWhat);

  private:
    NGramContainer       mData;
    size_t               mOrder;
    VocabularyTable*     mpPredictorTable;
    VocabularyTable*     mpTargetTable;
  }; // class NGramPool
  
  
} // namespace STK

//#include <STKLib/ContextSample.tcc>

#endif // ContextSample_h
// EOF
