#include <list>
#include <string>

namespace STK {
  /** 
   * @brief General string tokenizer
   */
  class Tokenizer 
  : public std::list<std::string>
  {
  public:
    // Constructors and Destructors ............................................
    Tokenizer(const char* pSeparator, bool skipEmpty = false)
    : std::list<std::string>(), mSeparator(pSeparator), mSkipEmpty(skipEmpty)
    {}

    Tokenizer(const char* pString, const char* pSeparator, bool skipEmpty = false)
    : std::list<std::string>(), mSeparator(pSeparator), mSkipEmpty(skipEmpty)
    { AddString(pString); }

    ~Tokenizer()
    {}

    /** 
     * @brief Parses a string and appends the tokens to the list
     * @param pString string to parse
     */
    void
    AddString(const char* pString);

    /** 
     * @brief Constant accessor to the separators string
     * @return Const refference
     */
    const std::string&
    Separator() const
    {return mSeparator;}

  private:
    std::string mSeparator;   ///< holds the list of separators
    bool        mSkipEmpty;   ///< if true, multiple separators will be regarded as one
  }; // class Tokenizer
}; // namespace STK


