/***************************************************************************
 *   copyright            : (C) 2004 by Lukas Burget,UPGM,FIT,VUT,Brno     *
 *   email                : burget@fit.vutbr.cz                            *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
  
#ifndef COMMON_H
#define COMMON_H

// include config.h if desired
#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Error.h"
#include "mymath.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <malloc.h>

#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#include <limits.h>
#include <string>
#include <sstream>



#if defined __CYGWIN32__ && !defined __CYGWIN__
   /* For backwards compatibility with Cygwin b19 and
      earlier, we define __CYGWIN__ here, so that
      we can rely on checking just for that macro. */
#  define __CYGWIN__  __CYGWIN32__
#endif

#if defined _WIN32 && !defined __CYGWIN__
   /* Use Windows separators on all _WIN32 defining
      environments, except Cygwin. */
#  define DIR_SEPARATOR_CHAR		'\\'
#  define DIR_SEPARATOR_STR		"\\"
#  define PATH_SEPARATOR_CHAR		';'
#  define PATH_SEPARATOR_STR		";"
#endif
#ifndef DIR_SEPARATOR_CHAR
   /* Assume that not having this is an indicator that all
      are missing. */
#  define DIR_SEPARATOR_CHAR		'/'
#  define DIR_SEPARATOR_STR		"/"
#  define PATH_SEPARATOR_CHAR		':'
#  define PATH_SEPARATOR_STR		":"
#endif /* !DIR_SEPARATOR_CHAR */


#if !HAVE_REENTRANT_SEARCH || !HAVE_SEARCH_H
#  include <gnusearch.h>
#else
#  include <search.h>
#endif

#include <vector>



/* Alignment of critical dynamic data structure
 *
 * Not all platforms support memalign so we provide a stk_memalign wrapper
 * void *stk_memalign( size_t align, size_t size, void **pp_orig )
 * *pp_orig is the pointer that has to be freed afterwards.
 */
#ifdef HAVE_POSIX_MEMALIGN
#  define stk_memalign(align,size,pp_orig) \
     ( !posix_memalign( pp_orig, align, size ) ? *(pp_orig) : NULL )
#  ifdef STK_MEMALIGN_MANUAL
#    undef STK_MEMALIGN_MANUAL
#  endif
#elif defined(HAVE_MEMALIGN)
   /* Some systems have memalign() but no declaration for it */
   //void * memalign( size_t align, size_t size );
#  define stk_memalign(align,size,pp_orig) \
     ( *(pp_orig) = memalign( align, size ) )
#  ifdef STK_MEMALIGN_MANUAL
#    undef STK_MEMALIGN_MANUAL
#  endif
#else /* We don't have any choice but to align manually */
#  define stk_memalign(align,size,pp_orig) \
     (( *(pp_orig) = malloc( size + align - 1 )) ? \
     (void *)( (((uintptr_t)*(pp_orig)) + 15) & ~0xFUL ) : NULL )
#  define STK_MEMALIGN_MANUAL
#endif

/// Define this, if new matrix will be used instead of plain pointer array
#define USE_NEW_MATRIX                                                  

#ifdef USE_NEW_MATRIX
#  define FEATURES_16_ALIGNED
#endif


// some extra code to handle strchr and strrchr 
#if !STDC_HEADERS
#  if !HAVE_STRCHR
#    define strchr index
#    define strrchr rindex
#  endif
#endif


// The following declaration assumes that SSE instructions are enabled
// and that we are using GNU C/C++ compiler, which defines the __attribute__ (...)
// notation.
//
// ENABLE_SSE is defined in <config.h>. Its value depends on options given
// in the configure phase of builidng the library
#if ENABLE_SSE && defined(__GNUC__ )
// vector of four single floats
typedef float  v4sf __attribute__((vector_size(16))); 
// vector of two single doubles
typedef double v2sd __attribute__((vector_size(16))); 

typedef union 
{
  v4sf    v;
  float   f[4];
} f4vector; 

typedef union 
{
  v2sd    v;
  double  f[2];
} d2vector; 
#endif // ENABLE_SSE && defined(__GNUC__ )





#ifndef M_PI
#  define M_PI 3.1415926535897932384626433832795
#endif

#define M_LOG_2PI 1.8378770664093454835606594728112

#define swap8(a) { \
  int t=((char*)&a)[0]; ((char*)&a)[0]=((char*)&a)[7]; ((char*)&a)[7]=t;\
      t=((char*)&a)[1]; ((char*)&a)[1]=((char*)&a)[6]; ((char*)&a)[6]=t;\
      t=((char*)&a)[2]; ((char*)&a)[2]=((char*)&a)[5]; ((char*)&a)[5]=t;\
      t=((char*)&a)[3]; ((char*)&a)[3]=((char*)&a)[4]; ((char*)&a)[4]=t;}
#define swap4(a) { \
  int t=((char*)&a)[0]; ((char*)&a)[0]=((char*)&a)[3]; ((char*)&a)[3]=t;\
      t=((char*)&a)[1]; ((char*)&a)[1]=((char*)&a)[2]; ((char*)&a)[2]=t;}
#define swap2(a) { \
  int t=((char*)&a)[0]; ((char*)&a)[0]=((char*)&a)[1]; ((char*)&a)[1]=t;}


#ifndef DOUBLEPRECISION
#  define DOUBLEPRECISION 1
#endif

#if DOUBLEPRECISION
#  define FLOAT double
#  define EPSILON DBL_EPSILON
#  define FLOAT_FMT "%lg"
#  define swapFLOAT swap8
#  define _ABS  my_fabs
#  define _EXP  my_exp
#  define _LOG  my_log
#  define _SQRT my_sqrt
#else
#  define FLOAT float
#  define EPSILON FLT_EPSILON
#  define FLOAT_FMT "%g"
#  define swapFLOAT swap4
#  define _ABS  my_fabsf
#  define _EXP  my_expf
#  define _LOG  my_logf
#  define _SQRT my_sqrtf
#endif


#define LOG_0     (-1.0e10)
#define LOG_MIN   (0.5 * LOG_0)
#define MIN_WEGIHT (1e-5)
#define MIN_LOG_WEGIHT (-11.5)
#define HIGHER_OF(a, b) ((a) > (b) ? (a) : (b))
#define LOWER_OF(a, b)  ((a) < (b) ? (a) : (b))
#define UNDEF_TIME -10000

#ifdef WIN32
#  define access _access
#endif


#define UINT_16  unsigned short
#define UINT_32  unsigned
#define INT_16   short
#define INT_32   int
#define FLOAT_32 float
#define DOUBLE_64 double

#define LOG_INC(a, b) ((a) = LogAdd((a),(b)))

#define PARAMKIND_WAVEFORM  0
#define PARAMKIND_LPC       1
#define PARAMKIND_LPREFC    2
#define PARAMKIND_LPCEPSTRA 3
#define PARAMKIND_LPDELCEP  4
#define PARAMKIND_IREFC     5
#define PARAMKIND_MFCC      6
#define PARAMKIND_FBANK     7
#define PARAMKIND_MELSPEC   8
#define PARAMKIND_USER      9
#define PARAMKIND_DISCRETE 10
#define PARAMKIND_PLP      11
#define PARAMKIND_ANON     12

#define PARAMKIND_E   0000100 /// has energy
#define PARAMKIND_N   0000200 /// absolute energy suppressed
#define PARAMKIND_D   0000400 /// has delta coefficients
#define PARAMKIND_A   0001000 /// has acceleration coefficients
#define PARAMKIND_C   0002000 /// is compressed
#define PARAMKIND_Z   0004000 /// has zero mean static coef.
#define PARAMKIND_K   0010000 /// has CRC checksum
#define PARAMKIND_0   0020000 /// has 0'th cepstral coef.
#define PARAMKIND_V   0040000 /// has VQ codebook index
#define PARAMKIND_T   0100000 /// has triple delta coefficients


/// Defines the white chars for string trimming
#if !defined(WHITE_CHARS)
#  define WHITE_CHARS " \t"
#endif


#include <sys/stat.h>

#if !defined(S_IFDIR)
#  if defined(_IFDIR)
#    define S_IFDIR _IFDIR
#  elif defined(_S_IFDIR)
#    define S_IFDIR _S_IFDIR
#  else
#    error "_IFDIR flag not found. Please check your sys/stat.h header file for compatible flag"
#  endif
#endif


//using namespace STK;

    
//namespace STK
//{
  /// This type will be used for flag passing
  typedef unsigned int  FlagType;

  
  /// Type for representing node's start and stop times
  typedef long long  TimingType;

  
  /** **************************************************************************
   ** **************************************************************************
   * @brief Aligns a number to a specified base
   * @param n Number of type @c _T to align
   * @return Aligned value of type @c _T
   */
  template<const size_t _align, typename _T>
    inline _T 
    align(const _T n)
    {
      const _T x(_align - 1); 
      return (n + x) & ~(x);
    }
  
  
  /** **************************************************************************
   ** **************************************************************************
   */
  enum PropagDirectionType
  {
    FORWARD,
    BACKWARD
  };
  
  
  /** **************************************************************************
   ** **************************************************************************
   */
  struct MyHSearchData 
  {
    ENTRY **            mpEntry;
    size_t              mNEntries;
  
    struct hsearch_data mTab;
    size_t              mTabSize;
  };
  
  
  /** **************************************************************************
   ** **************************************************************************
   */
  struct ReadlineData 
  {
    char*               buffer;
    int                 size;
  };
  
  
namespace STK 
{

  /** **************************************************************************
   ** **************************************************************************
   */
  class FileListElem
  {
  private:
    std::string         mLogical;     ///< Logical file name representation
    std::string         mPhysical;    ///< Pysical file name representation
    FLOAT               mWeight;
    
  public:
    FileListElem(const std::string & rFileName);
    ~FileListElem() {}
    
    const std::string &
    Logical() const { return mLogical; }

    const std::string &
    Physical() const { return mPhysical; }

    const FLOAT&
    Weight() const { return mWeight; }
        
    FileListElem *      mpNext;
    char *              mpPhysical;
    char                logical[1];
  };
  
  
  /** **************************************************************************
   ** **************************************************************************
   */
  typedef struct 
  {
    double logvalue; 
    unsigned negative:1;
  } FloatInLog;
  
  
  int     ReadParmKind(const char *str, bool checkBrackets);
  int     ParmKind2Str(unsigned parmKind, char *outstr);
  void    MakeFileName(char *outFileName, const char* inFileName,
                    const char *out_dir, const char *out_ext);
  
  char*   strtoupper(char *str);
  int     qsstrcmp(const void *a, const void *b);
  int     qsptrcmp(const void *a, const void *b);
  
  void    InitLogMath(void);
  double  LogAdd(double x, double y);
  double  LogSub(double x, double y);
  
  
  
  FloatInLog FIL_Add(FloatInLog a, FloatInLog b);
  FloatInLog FIL_Sub(FloatInLog a, FloatInLog b);
  FloatInLog FIL_Mul(FloatInLog a, FloatInLog b);
  FloatInLog FIL_Div(FloatInLog a, FloatInLog b);
  
  bool isBigEndian();
  
  void exp_vec(FLOAT *in, FLOAT *out, int size);
  void log_vec(FLOAT *in, FLOAT *out, int size);
  void sqrt_vec(FLOAT *in, FLOAT *out, int size);

  void softmax_vec(float *in, float *out, int size);
  void sigmoid_vec(float *in, float *out, int size);

  void softmax_vec(double *in, double *out, int size);
  void sigmoid_vec(double *in, double *out, int size);
  
  void fast_softmax_vec(float *in, float *out, int size);
  void fast_sigmoid_vec(float *in, float *out, int size);
  
  void fast_softmax_vec(double *in, double *out, int size);
  void fast_sigmoid_vec(double *in, double *out, int size);


  /** 
   * @brief 
   * 
   * @param f1 
   * @param f2 
   * @param nRounds 
   * 
   * @return 
   */
  bool 
  close_enough(const FLOAT f1, const FLOAT f2, int nRounds);
  
  FLOAT 
  float_safe_substract(const FLOAT& f1, const FLOAT &f2, int nRounds);

  int my_hcreate_r(size_t nel,
                  MyHSearchData *tab);
  int my_hsearch_r(ENTRY item, ACTION action, ENTRY **ret,
                  MyHSearchData *tab);
  void my_hdestroy_r (MyHSearchData *tab, int freeKeys);
  
  FILE *my_fopen(
    const char *file_name,
    const char *type,
    const char *filter);

  int IsSpace(int c);
  
  int my_fclose(FILE *fp);
  
  int fprintHTKstr(FILE *fp, const char *str);
  int getHTKstr(char *str, char **endPtrOrErrMsg);
  int skipHTKstr(const char *str, const char **endPtrOrErrMsg);
  char *expandFilterCommand(const char *command, const char *filename);
  char *readline(FILE *fp, struct ReadlineData *data);
  
  void InsertConfigParam(
    MyHSearchData *config_hash,
    const char *param_name,
    const char *value,
    int optionChar);
  
  void ReadConfig(const char *file_name, MyHSearchData *config_hash);
  void PrintConfig(MyHSearchData *config_hash);
  void CheckCommandLineParamUse(MyHSearchData *config_hash);
  
  const char *GetParamStr(
    MyHSearchData *config_hash,
    const char *param_name,
    const char *default_value);
  
  long GetParamInt(
    MyHSearchData *config_hash,
    const char *param_name,
    long default_value);
  
  FLOAT GetParamFlt(
    MyHSearchData *config_hash,
    const char *param_name,
    FLOAT default_value);
  
  bool GetParamBool(
    MyHSearchData *config_hash,
    const char *param_name,
    bool default_value);
  
  int GetParamEnum(
    MyHSearchData *config_hash,
    const char *param_name,
    int default_value, ...);
  
  int GetDerivParams(
    MyHSearchData *config_hash,
    int *derivOrder,
    int **derivWinLens,
    int *startFrmExt,
    int *endFrmExt,
    char **CMNPath,
    char **CMNFile,
    const char **CMNMask,
    char **CVNPath,
    char **CVNFile,
    const char **CVNMask,
    const char **CVGFile,
    const char *toolName,
    int pseudoModeule);
  
  int ParseOptions(
    int             argc,
    char*           argv[],
    const char*     optionMapping,
    const char*     toolName,
    MyHSearchData*  cfgHash);
  
  FileListElem **AddFileElem(FileListElem **last, char *fileElem);

  /**
   * @brief Converts initial portion of a string to a number using 'C' locale
   *
   * @param pStr String to convert
   * @param pValue Pointer where to store converted number
   * @param pRetStr A pointer to the character after the last character used in
   *                the conversion or @c pStr if no conversion is performed
   *
   * @return true if the conversion was successful, false otherwise
   *
   * Use this function instead of standard strtod(), strtof(), etc. if you want
   * it to behave independently of current locale settings
   */
  template <typename _T>
  bool Str2Number(const char *pStr, _T *pValue, char **pRetStr = 0)
  {
    std::streampos offset = 0;
    bool result = true;
  
    try
    {
      _T value;
      std::istringstream stream(pStr);
      stream.exceptions(std::ios::failbit | std::ios::badbit);
      stream.imbue(std::locale("C"));
      stream >> value;
      *pValue = value;
      offset = stream.tellg();
    }
    catch (...)
    {
      result = false;
    }
  
    if (pRetStr)
      *pRetStr = const_cast<char *>(pStr) + offset;
  
    return result;
  }

  // FIXME: specific function for float, this is necessary because using
  //        istringstream on Win32 in multithread process causes SIGSEGV
  bool Str2Number(const char *pStr, float *pValue, char **pRetStr = 0);

  // FIXME: specific function for double, this is necessary because using
  //        istringstream on Win32 in multithread process causes SIGSEGV
  bool Str2Number(const char *pStr, double *pValue, char **pRetStr = 0);

  /**
   * @brief Converts number to string using 'C' locale
   *
   * @param value Number to convert
   * @param pStr String where to store converted number
   * @param maxLen Maximum number of characters that can be stored to @c pStr
   * @param fieldWidth Minimum number of characters used in output string,
   *                   if converted number has less characters it is padded
   *                   with spaces from left
   * @param precision Number of digits after the decimal point
   *
   * @return String containing converted number, i.e. @c pStr
   *
   * Use this function instead of standard snprintf(), e.g. for float:
   * snprintf(pStr, maxLen, "%W.Pf", value); (W ... fieldWidth, P ... precision)
   * if you want it to behave independently of current locale settings
   */
  template <typename _T>
  char *Number2Str(_T value, char *pStr, size_t maxLen, int fieldWidth, int precision)
  {
    try
    {
      std::ostringstream stream;
      stream.width(fieldWidth);
      stream.precision(precision);
      stream.setf(std::ios::fixed, std::ios::floatfield);
      stream.exceptions(std::ios::failbit | std::ios::badbit);
      stream.imbue(std::locale("C"));
      stream << value;
      pStr[maxLen - 1] = '\0';
      strncpy(pStr, stream.str().c_str(), maxLen - 1);
    }
    catch (...)
    {
      pStr[0] = '\0';
    }
    return pStr;
  }

  /**
   * @brief Converts number to string using 'C' locale
   *
   * @param value Number to convert
   * @param pStr String where to store converted number
   * @param maxLen Maximum number of characters that can be stored to @c pStr
   *
   * @return String containing converted number, i.e. @c pStr
   *
   * Use this function instead of standard snprintf(), e.g. for float:
   * snprintf(pStr, maxLen, "%g", value);
   * if you want it to behave independently of current locale settings
   */
  template <typename _T>
  char *Number2Str(_T value, char *pStr, size_t maxLen)
  {
    try
    {
      std::ostringstream stream;
      stream.exceptions(std::ios::failbit | std::ios::badbit);
      stream.imbue(std::locale("C"));
      stream << value;
      pStr[maxLen - 1] = '\0';
      strncpy(pStr, stream.str().c_str(), maxLen - 1);
    }
    catch (...)
    {
      pStr[0] = '\0';
    }
    return pStr;
  }

  /**
   * @brief Reads number from the stream using 'C' locale
   *
   * @param pFile Pointer to the open stream
   * @param pValue Pointer where to store the number
   *
   * @return true if the operation was successful, false otherwise
   *
   * Use this function instead of standard fscanf() if you want
   * it to behave independently of current locale settings
   */
  template <typename _T>
  bool ReadNumber(FILE *pFile, _T *pValue)
  {
    char pnumber[256];
    pnumber[sizeof(pnumber) - 1] = '\0';
    if (fscanf(pFile, " %255[-+.0-9eE]", pnumber) != 1)
      return false;
    char *pend_ptr;
    bool result = Str2Number(pnumber, pValue, &pend_ptr);
    if (*pend_ptr != '\0')
      fseek(pFile, pend_ptr - pnumber - strlen(pnumber), SEEK_CUR);
    return result;
  }

  /**
   * @brief Writes number to the stream using 'C' locale
   *
   * @param pFile Pointer to the open stream
   * @param value Number to be written
   * @param fieldWidth Minimum number of characters used in output,
   *                   if converted number has less characters it is padded
   *                   with spaces from left
   * @param precision Number of digits after the decimal point
   *
   * @return Number of characters written. If an output error is encountered,
   *         a negative value is returned.
   *
   * Use this function instead of standard fprintf(), e.g. for float:
   * fprintf(pFile, "%W.Pf", value); (W ... fieldWidth, P ... precision)
   * if you want it to behave independently of current locale settings
   */
  template <typename _T>
  int WriteNumber(FILE *pFile, _T value, int fieldWidth, int precision)
  {
    char pnumber[256];
    Number2Str(value, pnumber, sizeof(pnumber), fieldWidth, precision);
    return fprintf(pFile, "%s", pnumber);
  }

  /**
   * @brief Writes number to the stream using 'C' locale
   *
   * @param pFile Pointer to the open stream
   * @param value Number to be written
   *
   * @return Number of characters written. If an output error is encountered,
   *         a negative value is returned.
   *
   * Use this function instead of standard fprintf(), e.g. for float:
   * fprintf(pFile, "%g", value);
   * if you want it to behave independently of current locale settings
   */
  template <typename _T>
  int WriteNumber(FILE *pFile, _T value)
  {
    char pnumber[256];
    Number2Str(value, pnumber, sizeof(pnumber));
    return fprintf(pFile, "%s", pnumber);
  }
}
  
  int npercents(const char *str);
  int process_mask(const char *normstr, const char *wildcard, char *substr);
  
  void fprintf_ll(FILE* fp, long long n);
  
namespace STK 
{

  //////////////////////////////////////////////////////////////////////////////
  // THE C++ COMMON ROUTINES
  //////////////////////////////////////////////////////////////////////////////
  /**
  *  @brief Returns true if rString matches rWildcard and fills substr with
  *         corresponding %%% matched pattern
  *  @param rString    String to be parsed
  *  @param rWildcard  String containing wildcard pattern
  *  @param Substr     The mathced %%% pattern is stored here
  *
  *  This is a C++ extension to the original process_mask function.
  *
  */
  bool
  ProcessMask(const std::string & rString,
              const std::string & rWildcard,
                    std::string & rSubstr);

  /**
  *  @brief Returns true if rString matches rWildcard and fills substr with
  *         corresponding %%% matched pattern
  *  @param rString    String to be parsed
  *  @param rWildcard  String containing wildcard pattern
  *  @param Substr     The mathced %%% pattern is stored here
  *
  *  This is a C++ extension to the original process_mask function.
  *
  */
  bool
  ProcessMask(const std::string & rString,
              const std::string & rWildcard,
              char*               pSubstr);
  
  /** 
   * @brief Wrap to the getHTKstr
   * @param rIn instring
   * @param rOut outstring
   */
  void
  ParseHTKString(const std::string & rIn, std::string & rOut);
  
  /**
   * @brief Splits the string into tokens based on the specified separator
   *
   * @param rString The original string
   * @param rTokens STL vector (not necesarilly empty) of strings
   * @param rSeparators String containing separators
   *
   * @return Returns a refference to the original vector
   */
  std::vector<std::string>&
  TokenizeString(
      const std::string&              rString, 
            std::vector<std::string>& rTokens, 
      const std::string&              rSeparators = ",");


  /** 
   * @brief Removes the leading and trailing white chars
   * 
   * @param rStr Refference to the string to be processed
   *
   * @return Refference to the original string
   *
   * The white characters are determined by the @c WHITE_CHARS macro defined 
   * above.
   */
  std::string&
  Trim(std::string& rStr);


  
  /// Sets application HTK compatibility. If true, all functions work to
  /// be HTK compatible.
  extern bool           gHtkCompatible;
  
  extern const char*    gpFilterWldcrd;
  extern const char*    gpScriptFilter;
  extern const char*    gpParmFilter;
  extern const char*    gpMmfFilter;
  extern const char*    gpHListOFilter;
  extern const char*    gpMmfOFilter;
  extern const char*    gpParmOFilter;


  struct ParameterRecord
  {
    std::string    mName;
    std::string    mType;
    std::string    mHelp;
  };

  extern std::vector<ParameterRecord>    gRegisteredParameters;

  void
  print_registered_parameters();
}
  
  int breakpoint(); // Empty function, just for debuging purposes

  
//}; //namespace STK

#endif // COMMON_H
