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

#include "Error.h"
#include "Math.h"

#include <stdio.h>
#include <float.h>
#include <math.h>
//#define __USE_GNU
#include <limits.h>
#include <search.h>
#include <string>




#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
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
#define DOUBLEPRECISION 1
#endif

#if DOUBLEPRECISION
#define FLOAT double
#define EPSILON DBL_EPSILON
#define FLOAT_FMT "%lg"
#define swapFLOAT swap8
#define _EXP  exp
#define _LOG  log
#define _SQRT sqrt
#else
#define FLOAT float
#define EPSILON FLT_EPSILON
#define FLOAT_FMT "%g"
#define swapFLOAT swap4
#define _EXP  expf
#define _LOG  logf
#define _SQRT sqrtf
#endif

//#define BOOL  int
//#define TRUE  1
//#define FALSE 0

#define LOG_0     (-1.0e10)
#define LOG_MIN   (0.5 * LOG_0)
#define MIN_WEGIHT (1e-5)
#define MIN_LOG_WEGIHT (-11.5)
#define HIGHER_OF(a, b) ((a) > (b) ? (a) : (b))
#define LOWER_OF(a, b)  ((a) < (b) ? (a) : (b))
#define UNDEF_TIME -10000

#ifdef WIN32
#define access _access
#endif


#define UINT_16 unsigned short
#define UINT_32 unsigned long
#define INT_16 short
#define INT_32 long
#define FLOAT_32 float

using namespace STK;

//namespace STK
//{
  /// This type will be used for flag passing
  typedef unsigned int    FlagType;
  
  
  typedef enum 
  {
    FALSE = 0,
    TRUE  = 1
  } BOOL;
  
  
  enum PropagDirectionType
  {
    FORWARD,
    BACKWARD
  };
  
  
  struct MyHSearchData 
  {
    ENTRY **            mpEntry;
    size_t              mNEntries;
  
  //private
    struct hsearch_data mTab;
    size_t              mTabSize;
  };
  
  
  struct ReadlineData 
  {
    char *              buffer;
    int                 size;
  };
  
  
  class FileListElem
  {
    std::string         mLogical;     ///< Logical file name representation
    std::string         mPhysical;    ///< Pysical file name representation
    
  public:
    FileListElem(const std::string & rFileName);
    ~FileListElem() {}
    
    const std::string &
    Logical() const { return mLogical; };

    const std::string &
    Physical() const { return mPhysical; };
        
    FileListElem *      mpNext;
    char *              mpPhysical;
    char                logical[1];
    
    
  };
  
  
  typedef struct 
  {
    double logvalue; 
    unsigned negative:1;
  } FloatInLog;
  
  //#define Error(...) _Error_(__func__, __FILE__, __LINE__, __VA_ARGS__)
  //#define Warning(...) _Warning_(__func__, __FILE__, __LINE__, __VA_ARGS__)
  
  //void _Error_(const char *func, const char *file, int line, char *msg, ...);
  //void _Warning_(const char *func, const char *file, int line, char *msg, ...);
  //void TraceLog(char *msg, ...);
  int     ReadParmKind(const char *str, bool checkBrackets);
  int     ParmKind2Str(int parmKind, char *outstr);
  void    MakeFileName(char *outFileName, const char* inFileName,
                    const char *out_dir, const char *out_ext);
  
  char *  strtoupper(char *str);
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
  
  void sigmoid_vec(FLOAT *in, FLOAT *out, int size);
  void exp_vec(FLOAT *in, FLOAT *out, int size);
  void log_vec(FLOAT *in, FLOAT *out, int size);
  void sqrt_vec(FLOAT *in, FLOAT *out, int size);
  void softmax_vec(FLOAT *in, FLOAT *out, int size);
  
  void fast_softmax_vec(float *in, float *out, int size);
  void fast_sigmoid_vec(float *in, float *out, int size);
  
  int my_hcreate_r(size_t nel,
                  MyHSearchData *tab);
  int my_hsearch_r(ENTRY item, ACTION action, ENTRY **ret,
                  MyHSearchData *tab);
  void my_hdestroy_r (MyHSearchData *tab, int freeKeys);
  
  FILE *my_fopen(
    const char *file_name,
    const char *type,
    const char *filter);
  
  int my_fclose(FILE *fp);
  
  int fprintHTKstr(FILE *fp, const char *str);
  int getHTKstr(char *str, char **endPtrOrErrMsg);
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
    int argc,
    char *argv[],
    const char *optionMapping,
    const char *toolName,
    MyHSearchData *cfgHash);
  
  FileListElem **AddFileElem(FileListElem **last, char *fileElem);
  
  int npercents(const char *str);
  int process_mask(const char *normstr, const char *wildcard, char *substr);
  
  
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
  #define PARAMKIND_0   0020000 /// has 0'th cepstral coef.  cccccc
  #define PARAMKIND_V   0040000 /// has VQ codebook index
  #define PARAMKIND_T   0100000 /// has triple delta coefficients
  
  extern const char *gpFilterWldcrd;
  extern const char *script_filter;
  extern const char *parm_filter;
  extern const char *MMF_filter;
  extern const char *parm_ofilter;
  extern const char *hlist_ofilter;
  extern const char *MMF_ofilter;
  
  /// Sets application HTK compatibility. If true, all functions work to
  /// be HTK compatible.
  extern bool       gHtkCompatible;
//}; //namespace STK

#endif // COMMON_H
