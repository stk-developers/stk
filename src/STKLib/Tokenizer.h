#include <list>
#include <string>

namespace STK {
  class Tokenizer 
  : public std::list<std::string>
  {
  public:
    typedef std::list<std::string> ListBaseType;

    // Constructors and Destructors ............................................
    Tokenizer(const char* pSeparator)
    : std::list<std::string>(), mSeparator(pSeparator)
    {}

    Tokenizer(const char* pString, const char* pSeparator)
    : std::list<std::string>(), mSeparator(pSeparator)
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
  }; // class Tokenizer
}; // namespace STK


