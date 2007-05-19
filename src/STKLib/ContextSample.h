#ifndef ContextSample_h
#define ContextSample_h

#include <iostream>

namespace STK
{
  class ContextSample;
  class ContextSamplePool;

  /** 
   * @brief 
   */
  class ContextSample
  {
  public:
    typedef FLOAT ProbType;
    typedef int   Predictee;
    typedef int   Predictor;

    //..........................................................................
    ContextSample();

    ContextSample(const ContextSample& rOrig);

    ~ContextSample()
    {}

    //..........................................................................
    /// Sample weight, soft counts, whatever
    const ProbType&
    Weight() const
    { return mWeight; }

    /// Returns number of predictors
    const size_t
    NPtors() const;

    /// Returns reference to the i'th predictor
    Predictor&
    rPtor(size_t i)
    { return mpPredictors[i]; }

    /// Returns reference to the predictee
    Predictee&
    rPtee()
    { return mValue; }

    friend std::ostream&
    operator<<(std::ostream& rOstr, const ContextSample& rCS);

    
    friend class ContextSamplePool;

  private:
    Predictee     mValue;         ///< Value of the sample
    Predictor*    mpPredictors;   ///< Pointer to array. ContexSamplePool (de)allocates this one
    ProbType      mWeight;        ///< Sample's weight
    
  }; // class ContextSample


  /** 
   * @brief Holds samples
   */
  class ContextSamplePool
  {
  public:
    ContextSamplePool();

    ContextSamplePool(const ContextSamplePool& rOrig);

    ContextSamplePool(size_t n);

    ~ContextSamplePool();

    void 
    LoadFrom(std::istream& rStream);

    /// Returns context size
    const size_t
    NPtors() const;

    /// Returns number of samples in the pool
    const size_t
    Size() const
    { return mData.size(); }


    /// Dumps the whole pool to the stream
    friend std::ostream&
    operator <<(std::ostream& rOstr, const ContextSamplePool& rWhat);

    /// Reads data from a stream
    friend std::istream&
    operator >>(std::istream& rIstr, ContextSamplePool& rWhat);

  private:
    std::vector<ContextSample>  mData;
  }; // class ContextSamplePool
} // namespace STK

//#include <STKLib/ContextSample.tcc>

#endif // ContextSample_h
// EOF
