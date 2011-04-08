

/** @file stkstream.h
 *  This is an STK C++ Library header.
 */


#ifndef STK_StkStream_h
#define STK_StkStream_h

#include <fstream>
#include <string>


#pragma GCC system_header


//extern const char * gpFilterWldcrd;

namespace STK
{

  /**
   *  @brief Expands a filter command into a runnable form
   *
   *  This function replaces all occurances of *filter_wldcard in *command by
   *  *filename
   */
  //char * ExpandFilterCommand(const char *command, const char *filename);

  /**
   *  @brief Provides a layer of compatibility for C/POSIX.
   *
   *  This GNU extension provides extensions for working with standard C
   *  FILE*'s and POSIX file descriptors.  It must be instantiated by the
   *  user with the type of character used in the file stream, e.g.,
   *  basic_stkbuf<char>.
   */
  template<
    typename _CharT, 
    typename _Traits = std::char_traits<_CharT>
  > 
    class basic_stkbuf : public std::basic_filebuf<_CharT, _Traits>
    {
    public:

      typedef basic_stkbuf<_CharT, _Traits>     this_type;

      // Types:
      typedef _CharT                        char_type;
      typedef _Traits                       traits_type;

      typedef typename traits_type::int_type        int_type;
      typedef typename traits_type::pos_type        pos_type;
      typedef typename traits_type::off_type        off_type;
      typedef std::size_t                           size_t;

    public:

      /// @{
      /// Type of streambuffer
      static const unsigned int t_undef  = 0; ///< undefined
      static const unsigned int t_file   = 1; ///< file stream
      static const unsigned int t_pipe   = 2; ///< pipe
      static const unsigned int t_filter = 4; ///< filter
      static const unsigned int t_stdio  = 8; ///< standard input/output
      /// @}

    public:

      /**
      * deferred initialization
      */
      basic_stkbuf() : std::basic_filebuf<_CharT, _Traits>(),
        mFilename(""), mpFilePtr(0), mStreamType(t_undef){}

      /**
       *  @brief  Opens a stream.
       *  @param  fName  The name of the file.
       *  @param  m      The open mode flags.
       *  @param  pFilter The pFilter command to use
       *  @return  @c this on success, NULL on failure
       *
       *  If a file is already open, this function immediately fails.
       *  Otherwise it tries to open the file named @a s using the flags
       *  given in @a mode.
       *
       *  [Table 92 gives the relation between openmode combinations and the
       *  equivalent fopen() flags, but the table has not been copied yet.]
       */
      basic_stkbuf(const char* fName, std::ios_base::openmode m, const char* pFilter="");

      
      /**
      *  @return  The underlying FILE*.
      *
      *  This function can be used to access the underlying "C" file pointer.
      *  Note that there is no way for the library to track what you do
      *  with the file, so be careful.
      */
      std::__c_file*
      file() { return this->_M_file.file(); }


      /**
      *  @return  The underlying FILE*.
      *
      *  This function can be used to access the underlying "C" file pointer.
      *  Note that there is no way for the library to track what you do
      *  with the file, so be careful.
      */
      std::__c_file*
      fp() { return this->_M_file.file(); }
      
      
      /**
       *  @brief  Opens an external file.
       *  @param  fName  The name of the file.
       *  @param  m      The open mode flags.
       *  @param  pFilter The pFilter command to use
       *  @return  @c this on success, NULL on failure
       *
       *  If a file is already open, this function immediately fails.
       *  Otherwise it tries to open the file named @a s using the flags
       *  given in @a mode.
       *
       *  [Table 92 gives the relation between openmode combinations and the
       *  equivalent fopen() flags, but the table has not been copied yet.]
       */
      this_type*
      open(const char* pFName, std::ios_base::openmode m, const char* pFilter="");
      
      /**
       *  @brief  Closes the currently associated file.
       *  @return  @c this on success, NULL on failure
       *
       *  If no file is currently open, this function immediately fails.
       *
       *  If a "put buffer area" exists, @c overflow(eof) is called to flush
       *  all the characters.  The file is then closed.
       *
       *  If any operations fail, this function also fails.
       */
      this_type*
      close();

      /**
      *  Closes the external data stream if the file descriptor constructor
      *  was used.
      */
      virtual
      ~basic_stkbuf() 
      {close();};

      /// Returns the file name
      const std::string 
      name() const 
      {return mFilename;}


    private:
      /// converts the ios::xxx mode to stdio style
      static void open_mode(std::ios_base::openmode __mode, int&, int&,  char* __c_mode);

      /**
       *  @param  __f  An open @c FILE*.
       *  @param  __mode  Same meaning as in a standard filebuf.
       *  @param  __size  Optimal or preferred size of internal buffer, in chars.
       *                Defaults to system's @c BUFSIZ.
       *
       *  This method associates a file stream buffer with an open
       *  C @c FILE*.  The @c FILE* will not be automatically closed when the
       *  basic_stkbuf is closed/destroyed. It is equivalent to one of the constructors
       *  of the stdio_filebuf class defined in GNU ISO C++ ext/stdio_filebuf.h
      */
      void superopen(std::__c_file* __f, std::ios_base::openmode __mode,
            size_t __size = static_cast<size_t>(BUFSIZ));


    private:
      /// Holds the full file name
      std::string           mFilename;

      std::ios_base::openmode  mMode;

      /// Holds a pointer to the main FILE structure
      FILE *                mpFilePtr;

      /// tells what kind of stream we use (stdio, file, pipe)
      unsigned int          mStreamType;

    };



  /**
   *  @brief This extension wraps stkbuf stream buffer into the standard ios class.
   *
   *  This class is inherited by (i/o)stkstream classes which make explicit use of
   *  the custom stream buffer
   */
  template<
    typename _CharT, 
    typename _Traits = std::char_traits<_CharT>
  >	
    class BasicStkIos 
    : virtual public std::basic_ios<_CharT, _Traits>
    {
    public:
      typedef basic_stkbuf <_CharT,_Traits>        StkBufType;

      BasicStkIos() 
      : mBuf() 
      { init(&mBuf) ;};

      BasicStkIos(const char* fName, std::ios::openmode m, const char* pFilter) 
      : mBuf(fName, m, pFilter) 
      { init(&mBuf) ; }

      StkBufType*
      rdbuf() 
      { return &mBuf; }

    protected:
      StkBufType  mBuf;
    };


  /**
   *  @brief  Controlling input for files.
   *
   *  This class supports reading from named files, using the inherited
   *  functions from std::istream.  To control the associated
   *  sequence, an instance of std::stkbuf is used.
  */
  template<
    typename _CharT, 
    typename _Traits = std::char_traits<_CharT>
  >	
    class BasicIStkStream 
    : public BasicStkIos<_CharT, _Traits>, 
      public std::basic_istream<_CharT, _Traits>
    {
    public:
      typedef BasicStkIos<_CharT, _Traits> BasicStkIosType;
      typedef std::basic_istream<_CharT,_Traits>  IStreamType;


      // Constructors:
      /**
       *  @brief  Default constructor.
       *
       *  Initializes @c mBuf using its default constructor, and passes
       *  @c &sb to the base class initializer.  Does not open any files
       *  (you haven't given it a filename to open).
       */
      BasicIStkStream() 
      : BasicStkIosType(),
        IStreamType(BasicStkIosType::rdbuf())
      {};

      /**
       *  @brief  Create an input file stream.
       *  @param  fName  String specifying the filename.
       *  @param  m      Open file in specified mode (see std::ios_base).
       *  @param  pFilter String specifying pFilter command to use on fName
       *
       *  @c ios_base::in is automatically included in
       *  @a m.
       *
       *  Tip:  When using std::string to hold the filename, you must use
       *  .c_str() before passing it to this constructor.
      */
      BasicIStkStream(const char* pFName, std::ios::openmode m=std::ios::out, const char* pFilter="")
      : BasicStkIosType(),
        IStreamType(BasicStkIosType::rdbuf())
      {this->open(pFName, std::ios::in, pFilter);}

      ~BasicIStkStream() 
      {
        this->close();
      }
        
      /**
      *  @brief  Opens an external file.
      *  @param  s  The name of the file.
      *  @param  mode  The open mode flags.
      *  @param  pFilter The pFilter command to use
         *
      *  Calls @c std::basic_filebuf::open(s,mode|in).  If that function
      *  fails, @c failbit is set in the stream's error state.
      *
      *  Tip:  When using std::string to hold the filename, you must use
      *  .c_str() before passing it to this constructor.
      */
      void open(const char* pFName, std::ios::openmode m=std::ios::in, const char* pFilter = "")
      {
        if (!BasicStkIosType::mBuf.open(pFName, m | std::ios_base::in, pFilter)) {
          this->setstate(std::ios_base::failbit);
        }
        else {
        // Closing an fstream should clear error state
          BasicStkIosType::clear();
        }
      }

      /**
      *  @brief  Returns true if the external file is open.
      */
      bool is_open() const {return BasicStkIosType::mBuf.is_open();}


      /**
      *  @brief  Closes the stream
      */
      void close() {BasicStkIosType::mBuf.close();}

      /**
      *  @brief  Returns the filename
      */
      const std::string name() const {return BasicStkIosType::mBuf.name();}

      /// Returns a pointer to the main FILE structure
      std::__c_file*
      file() {return BasicStkIosType::mBuf.file();}

      /// Returns a pointer to the main FILE structure
      std::__c_file*
      fp() {return BasicStkIosType::mBuf.fp();}

      // /**
      //  *  @brief  Reads a single line
      //  *
      //  *  This is a specialized function as std::getline does not provide a way to
      //  *  read multiple end-of-line symbols (we need both '\n' and EOF to delimit
      //  *  the line)
      //  */
      // void
      // GetLine(string& rLine);

    }; // class BasicIStkStream


  /**
   *  @brief  Controlling output for files.
   *
   *  This class supports reading from named files, using the inherited
   *  functions from std::ostream.  To control the associated
   *  sequence, an instance of STK::stkbuf is used.
  */
  template<
    typename _CharT, 
    typename _Traits = std::char_traits<_CharT>
  >	
    class BasicOStkStream 
    : public BasicStkIos<_CharT, _Traits>,
      public std::basic_ostream<_CharT, _Traits>
    {
    public:
      typedef BasicStkIos<_CharT, _Traits> BasicStkIosType;
      typedef std::basic_ostream<_CharT,_Traits>  OStreamType;

      // Constructors:
      /**
       *  @brief  Default constructor.
       *
       *  Initializes @c sb using its default constructor, and passes
       *  @c &sb to the base class initializer.  Does not open any files
       *  (you haven't given it a filename to open).
       */
      BasicOStkStream() 
      : BasicStkIosType(),
        OStreamType(BasicStkIosType::rdbuf())
      {};

      /**
       *  @brief  Create an output file stream.
       *  @param  fName  String specifying the filename.
       *  @param  m      Open file in specified mode (see std::ios_base).
       *  @param  pFilter String specifying pFilter command to use on fName
       *
       *  @c ios_base::out|ios_base::trunc is automatically included in
       *  @a mode.
       *
       *  Tip:  When using std::string to hold the filename, you must use
       *  .c_str() before passing it to this constructor.
       */
      BasicOStkStream(const char* pFName, std::ios::openmode m=std::ios::out, const char* pFilter="")
      : BasicStkIosType(),
        OStreamType(BasicStkIosType::rdbuf())
      {this->open(pFName, std::ios::out, pFilter);}

      /**
      *  @brief  Opens an external file.
      *  @param  fName  The name of the file.
      *  @param  m  The open mode flags.
      *  @param  pFilter String specifying pFilter command to use on fName
      *
      *  Calls @c std::basic_filebuf::open(s,mode|out).  If that function
      *  fails, @c failbit is set in the stream's error state.
      *
      *  Tip:  When using std::string to hold the filename, you must use
      *  .c_str() before passing it to this constructor.
      */
      void open(const char* pFName, std::ios::openmode m=std::ios::out, const char* pFilter="")
      {
        if (!BasicStkIosType::mBuf.open(pFName, m | std::ios_base::out, pFilter))
          this->setstate(std::ios_base::failbit);
        else
        // Closing an fstream should clear error state
          this->clear();
      }

      /**
      *  @brief  Returns true if the external file is open.
      */
      bool is_open() const 
      { return BasicStkIosType::mBuf.is_open();}

      /**
      *  @brief  Closes the stream
      */
      void close() 
      { BasicStkIosType::mBuf.close();}

      /**
      *  @brief  Returns the filename
      */
      const std::string name() const 
      { return BasicStkIosType::mBuf.name();}

      /// Returns a pointer to the main FILE structure
      std::__c_file*
      file() 
      { return BasicStkIosType::mBuf.file();}

      /// Returns a pointer to the main FILE structure
      std::__c_file*
      fp() 
      { return BasicStkIosType::mBuf.fp();}

    }; // class BasicOStkStream


  /**
   * We define some implicit stkbuf class
   */
  ///@{
#ifndef _GLIBCPP_USE_WCHAR_T
  typedef BasicOStkStream<char>      OStkStream;
  typedef BasicOStkStream<wchar_t>  WOStkStream;
  typedef BasicIStkStream<char>      IStkStream;
  typedef BasicIStkStream<wchar_t>  WIStkStream;
#else 
  typedef BasicOStkStream<char>     WOStkStream;
  typedef BasicOStkStream<wchar_t>  WOStkStream;
  typedef BasicIStkStream<char>     WIStkStream;
  typedef BasicIStkStream<wchar_t>  WIStkStream;
#endif
  /// @}


}; // namespace STK

# include "stkstream.tcc"

// STK_StkStream_h
#endif 
