/** @file MlfStream.h
 *  This is an STK C++ Library header.
 *
 *  The naming convention in this file coppies the std::* naming as well as STK
 */


#ifndef STK_MlfStream_h
#define STK_MlfStream_h

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <set>


namespace STK
{
  class LabelRecord;
  class LabelContainer;


  /// this container stores the lables in linear order as they came
  /// i.e. they cannot be hashed
  typedef  std::list< std::pair<std::string,LabelRecord> *> LabelListType;

  /// type of the container used to store the labels
  typedef  std::map<std::string, LabelRecord>               LabelHashType;



  /**
   *  @brief Describes type of MLF definition
   *
   *  See HTK book for MLF structure. Terms used in STK are
   *  compatible with those in HTK book.
   */
  enum MlfDefType
  {
    MLF_DEF_UNKNOWN = 0,              ///< unknown definition
    MLF_DEF_IMMEDIATE_TRANSCRIPTION,  ///< immediate transcription
    MLF_DEF_SUB_DIR_DEF               ///< subdirectory definition
  };



  /** **************************************************************************
   *  @brief Holds association between label and stream
   */
  class LabelRecord
  {

  public:
    LabelRecord() : miLabelListLimit(NULL)
    { }

    ~LabelRecord()
    { }

    /// definition type
    MlfDefType                mDefType;

    /// position of the label in the stream
    std::streampos            mStreamPos;

    /**
     *  @brief points to the current end of the LabelList
     *
     *  The reason for storing this value is to know when we inserted
     *  a label into the hash. It is possible, that the hash label came
     *  after list label, in which case the list label is prefered
     */
    LabelListType::iterator   miLabelListLimit;

  };




  /**
   *  @brief Provides an interface to label hierarchy and searching
   *
   *  This class stores label files in a map structure. When a wildcard
   *  convence is used, the class stores the labels in separate maps according
   *  to level of wildcard abstraction. By level we mean the directory structure
   *  depth.
   */
  class LabelContainer
  {
  public:
    /// The constructor
    LabelContainer() : mUseHashedSearch(true) {}

    /// The destructor
    ~LabelContainer();

    /**
     *  @brief Inserts new label to the hash structure
     */
    void
    Insert(
      const std::string &      rLabel,
      std::streampos           Pos);


    /**
     *  @brief Looks for a record in the hash
     */
    bool
    FindInHash(
      const std::string&        rLabel,
      LabelRecord&              rLS);

    /**
     *  @brief Looks for a record in the list
     *  @param rLabel Label to look for
     *  @param rLS    Structure to fill with found data
     *  @param limitSearch If true @p rLS's @c mLabelListLimit gives the limiting position in the list
     */
    bool
    FindInList(
      const std::string&        rLabel,
      LabelRecord&              rLS,
      bool                      limitSearch = false);

    /**
     *  @brief Looks for a record
     */
    bool
    Find(
      const std::string &       rLabel,
      LabelRecord &             rLS);

    /**
     *  @brief Returns the matched pattern
     */
    const std::string &
    MatchedPattern() const
    {
      return mMatchedPattern;
    }

    /**
     *  @brief Returns the matched pattern mask (%%%)
     */
    const std::string &
    MatchedPatternMask() const
    {
      return mMatchedPatternMask;
    }

    /** 
     * @brief Writes contents to stream (text)
     * @param rOStream stream to write to
     */
    void
    Write(std::ostream& rOStream);

  private:
    /// type used for directory depth notation
    typedef  size_t                 DepthType;


    /// this set stores depths of * labels observed at insertion
    std::set<DepthType>             mDepths;

    /// stores the labels
    LabelHashType                   mLabelMap;
    LabelListType                   mLabelList;

    /// true if labels are to be sought by hashing function (fast) or by
    /// sequential search (slow)
    bool                            mUseHashedSearch;

    /// if Find matches the label, this var stores the pattern that matched the
    /// query
    std::string                     mMatchedPattern;

    /// if Find matches the label, this var stores the the masked characters.
    /// The mask is given by '%' symbols
    std::string                     mMatchedPatternMask;

    /**
     *  @brief Returns the directory depth of path
     */
    size_t
    DirDepth(const std::string & path);


  };


  /** 
   * @brief MLF output buffer definition
   */
  template<
    typename _CharT, 
    typename _Traits = std::char_traits<_CharT>,
    typename _CharTA = std::allocator<_CharT>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT> 
  > 
    class BasicOMlfStreamBuf 
    : public std::basic_streambuf<_CharT, _Traits> 
    {
    public:
      // necessary typedefs ....................................................
      typedef BasicOMlfStreamBuf<_CharT,_Traits,_CharTA,ByteT,ByteAT>
                            this_type; 
      typedef std::basic_ostream<_CharT, _Traits>& 
                            OStreamReference;
      typedef std::basic_streambuf<_CharT, _Traits>
                            StreamBufType;
      typedef _CharTA       char_allocator_type;
      typedef _CharT        char_type;
      typedef typename _Traits::int_type int_type;
      typedef typename _Traits::pos_type pos_type;
      typedef ByteT         byte_type;
      typedef ByteAT        byte_allocator_type; 
      typedef byte_type*    byte_buffer_type;
      typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
      typedef std::vector<char_type, char_allocator_type > char_vector_type;


      BasicOMlfStreamBuf(OStreamReference rOStream, size_t bufferSize);

      ~BasicOMlfStreamBuf();

      // virtual functions inherited from basic_streambuf.......................
      int 
      sync();

      /** 
       * @brief Write character in the case of overflow
       * @param c Character to be written.
       * @return A value different than EOF (or traits::eof() for other traits) 
       *         signals success.  If the function fails, either EOF 
       *         (or traits::eof() for other traits) is returned or an 
       *         exception is thrown.
       */
      int_type
      overflow(int_type c = _Traits::eof());


      // MLF specific functions ................................................
      /** 
       * @brief Creates a new MLF block
       * @param rFileName filename to be opened
       */
      this_type*
      Open(const std::string& rFileName);

      /** 
       * @brief Closes MLF block
       */
      void
      Close();

      /** 
       * @brief Returns true if the MLF is now in open state
       */
      bool
      IsOpen() const
      { return mIsOpen; }

      LabelContainer&
      rLabels()
      { return mLabels; }

    private:
      bool             mIsOpen;
      char_type        mLastChar;
      OStreamReference mOStream;
      LabelContainer   mLabels;
    }; // class BasicOMlfStreamBuf



  /** 
   * @brief MLF input buffer definition
   */
  template<
    typename _CharT, 
    typename _Traits = std::char_traits<_CharT>,
    typename _CharTA = std::allocator<_CharT>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT> 
  > 
    class BasicIMlfStreamBuf 
    : public std::basic_streambuf<_CharT, _Traits> 
    {
    private:
      // internal automaton states
      static const int IN_HEADER_STATE   = 0;
      static const int OUT_OF_BODY_STATE = 1;
      static const int IN_TITLE_STATE    = 2;
      static const int IN_BODY_STATE     = 3;


    public: // necessary typedefs ..............................................
      typedef BasicIMlfStreamBuf<_CharT,_Traits,_CharTA,ByteT,ByteAT>
                            this_type; 
      typedef std::basic_istream<_CharT, _Traits>& IStreamReference;
      typedef std::basic_streambuf<_CharT, _Traits>
                            StreamBufType;
      typedef _CharTA       char_allocator_type;
      typedef _CharT        char_type;
      typedef typename _Traits::int_type int_type;
      typedef typename _Traits::pos_type pos_type;
      typedef ByteT         byte_type;
      typedef ByteAT        byte_allocator_type; 
      typedef byte_type*    byte_buffer_type;
      typedef std::vector<byte_type, byte_allocator_type > byte_vector_type;
      typedef std::vector<char_type, char_allocator_type > char_vector_type;


    public:
      // constructors and destructors ..........................................
      BasicIMlfStreamBuf(IStreamReference rIStream, size_t bufferSize = 1024);

      ~BasicIMlfStreamBuf();

      // virtual functions inherited from basic_streambuf.......................
      /** 
       * @brief Get character in the case of underflow
       * 
       * @return The new character available at the get pointer position, if 
       *         any. Otherwise, traits::eof() is returned.  
       */
      int_type
      underflow();


      // MLF specific functions ................................................
      /** 
       * @brief Creates a new MLF block
       * @param rFileName filename to be opened
       */
      this_type*
      Open(const std::string& rFileName);

      /** 
       * @brief Closes MLF block
       */
      this_type*
      Close();

      /** 
       * @brief Returns true if the MLF is now in open state
       */
      bool
      IsOpen() const
      { return mIsOpen; }

      /** 
       * @brief Parses the stream (if possible) and stores positions to the 
       *        label titles
       */
      void
      Index();

      const bool
      IsHashed() const
      { return mIsHashed; }

      /** 
       * @brief Jumps to next label definition
       * @param rName std::string to be filled with the label name
       * @return true on success
       *
       * The procedure automatically tries to hash the labels.
       */
      bool
      JumpToNextDefinition(std::string& rName);


    private: // auxillary functions ............................................
      /** 
       * @brief Fills the line buffer with next line and updates the internal
       * state of the finite automaton
       */
      void
      FillLineBuffer();


    private: // atributes ......................................................
      // some flags
      bool              mIsOpen;
      bool              mIsHashed;
      bool              mIsEof;

      /// internal state of the finite automaton
      int               mState;

      IStreamReference  mIStream;
      LabelContainer    mLabels;

      std::vector<char_type>  mLineBuffer;
    }; // class BasicIMlfStreamBuf




  /** 
   * @brief Base class with type-independent members for the Mlf Output 
   *        Stram class
   *
   * This is a derivative of the basic_ios class. We derive it as we need 
   * to override some member functions
   */
  template<
    typename Elem, 
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
  >	
    class BasicOMlfStreamBase
    : virtual public std::basic_ios<Elem,Tr>
    {
    public:
      typedef std::basic_ostream<Elem, Tr>& OStreamReference;
      typedef BasicOMlfStreamBuf <
        Elem,Tr,ElemA,ByteT,ByteAT> OMlfStreamBufType;

      /** 
       * @brief constructor
       * 
       * @param rOStream user defined output stream 
       */
      BasicOMlfStreamBase(OStreamReference rOStream, 
          size_t bufferSize)
      : mBuf(rOStream, bufferSize)
      { init(&mBuf); };
      
      /** 
       * @brief Returns a pointer to the buffer object for this stream
       */
      OMlfStreamBufType* 
      rdbuf() 
      { return &mBuf; };

    private:
      OMlfStreamBufType mBuf;
    };  


  template<
    typename Elem, 
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
  >	
    class BasicIMlfStreamBase
    : virtual public std::basic_ios<Elem,Tr>
    {
    public:
      typedef std::basic_istream<Elem, Tr>& IStreamReference;
      typedef BasicIMlfStreamBuf <
        Elem,Tr,ElemA,ByteT,ByteAT> IMlfStreamBufType;

      BasicIMlfStreamBase( IStreamReference rIStream,
          size_t bufferSize)
      : mBuf(rIStream, bufferSize)
      { init(&mBuf ); };
      
      IMlfStreamBufType* 
      rdbuf() 
      { return &mBuf; };

    private:
      IMlfStreamBufType mBuf;
    };


  template<
    typename Elem, 
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
  >
    class BasicOMlfStream 
    : public BasicOMlfStreamBase<Elem,Tr,ElemA,ByteT,ByteAT>, 
      public std::basic_ostream<Elem,Tr>
    {
    public:
      typedef BasicOMlfStreamBase< Elem,Tr,ElemA,ByteT,ByteAT> 
                                          BasicOMlfStreamBaseType;
      typedef std::basic_ostream<Elem,Tr> OStreamType;
      typedef OStreamType&                OStreamReference;

      BasicOMlfStream(OStreamReference rOStream, size_t bufferSize = 32)
      : BasicOMlfStreamBaseType(rOStream, bufferSize), 
        OStreamType(BasicOMlfStreamBaseType::rdbuf())
      { }

      /** 
       * @brief Destructor closes the stream
       */
      ~BasicOMlfStream()
      { }


      /** 
       * @brief Creates a new MLF block
       * @param rFileName filename to be opened
       */
      void
      Open(const std::string& rFileName)
      { BasicOMlfStreamBaseType::rdbuf()->Open(rFileName); }

      /** 
       * @brief Closes MLF block
       */
      void
      Close()
      { BasicOMlfStreamBaseType::rdbuf()->Close(); }

      /** 
       * @brief Returns true if the MLF is now in open state
       */
      bool
      IsOpen() const
      { return BasicOMlfStreamBaseType::rdbuf()->IsOpen(); }

      /** 
       * @brief Accessor to the label container
       * @return Reference to the label container
       */
      LabelContainer&
      rLabels()
      { return BasicOMlfStreamBaseType::rdbuf()->rLabels(); }
    };



  template<
    typename Elem, 
    typename Tr = std::char_traits<Elem>,
    typename ElemA = std::allocator<Elem>,
    typename ByteT = char,
    typename ByteAT = std::allocator<ByteT>
  >	
    class BasicIMlfStream 
    : public BasicIMlfStreamBase<Elem,Tr,ElemA,ByteT,ByteAT>, 
      public std::basic_istream<Elem,Tr>
    {
    public:
      typedef BasicIMlfStreamBase <Elem,Tr,ElemA,ByteT,ByteAT> 
                                          BasicIMlfStreamBaseType;
      typedef std::basic_istream<Elem,Tr> IStreamType;
      typedef IStreamType&                IStreamReference;
      typedef unsigned char               byte_type;

      BasicIMlfStream(IStreamReference rIStream, size_t bufferSize = 32)
      : BasicIMlfStreamBaseType(rIStream, bufferSize), 
        IStreamType(BasicIMlfStreamBaseType::rdbuf())
      {};


      /** 
       * @brief Creates a new MLF block
       * @param rFileName filename to be opened
       */
      void
      Open(const std::string& rFileName)
      { 
        std::basic_streambuf<Elem, Tr>* p_buf;

        p_buf = BasicIMlfStreamBaseType::rdbuf()->Open(rFileName);

        if (NULL == p_buf) {
          IStreamType::clear(IStreamType::rdstate() | std::ios::failbit);
        }
        else {
          IStreamType::clear();
        }
      }

      /** 
       * @brief Closes MLF block.
       * In fact, nothing is done
       */
      void 
      Close()
      { 
        if (NULL == BasicIMlfStreamBaseType::rdbuf()->Close()) {
          IStreamType::clear(IStreamType::rdstate() | std::ios::failbit);
        }
      }

      void
      Index()
      { BasicIMlfStreamBaseType::rdbuf()->Index(); }

      const bool
      IsHashed() const
      { return BasicIMlfStreamBaseType::rdbuf()->IsHashed(); }

    };



  // MAIN TYPEDEFS..............................................................
  typedef BasicOMlfStream<char>     OMlfStream;
  typedef BasicOMlfStream<wchar_t>  WOMlfStream;
  typedef BasicIMlfStream<char>     IMlfStream;
  typedef BasicIMlfStream<wchar_t>  WIMlfStream;


#ifdef PATH_MAX
  const size_t MAX_LABEL_DEPTH = PATH_MAX;
#else
  const size_t MAX_LABEL_DEPTH = 1024;
#endif


}; // namespace STK

#include "MlfStream.tcc"

#endif
