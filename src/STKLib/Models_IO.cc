#include "Models.h"
#include "mymath.h"
#include "SigP.h"
#include <stdarg.h>
 
namespace STK
{
//  int             gCurrentMmfLine       = 1;
//  const char *    gpCurrentMmfName;
//  bool            gHmmReadBinary;
//  bool            gStringUnget          = false;
//  const char *    gpKwds[KID_MaxKwdID]  = {0};
  
  
  
  //###########################################################################
  //###########################################################################
  // GENERAL FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  InitKwdTable()
  {
    mpKwds[KID_BeginHMM    ] = "BeginHMM";      mpKwds[KID_Use        ] = "Use";
    mpKwds[KID_EndHMM      ] = "EndHMM";        mpKwds[KID_NumMixes   ] = "NumMixes";
    mpKwds[KID_NumStates   ] = "NumStates";     mpKwds[KID_StreamInfo ] = "StreamInfo";
    mpKwds[KID_VecSize     ] = "VecSize";       mpKwds[KID_NullD      ] = "NullD";
    mpKwds[KID_PoissonD    ] = "PoissonD";      mpKwds[KID_GammaD     ] = "GammaD";
    mpKwds[KID_RelD        ] = "RelD";          mpKwds[KID_GenD       ] = "GenD";
    mpKwds[KID_DiagC       ] = "DiagC";         mpKwds[KID_FullC      ] = "FullC";
    mpKwds[KID_XformC      ] = "XformC";        mpKwds[KID_State      ] = "State";
    mpKwds[KID_TMix        ] = "TMix";          mpKwds[KID_Mixture    ] = "Mixture";
    mpKwds[KID_Stream      ] = "Stream";        mpKwds[KID_SWeights   ] = "SWeights";
    mpKwds[KID_Mean        ] = "Mean";          mpKwds[KID_Variance   ] = "Variance";
    mpKwds[KID_InvCovar    ] = "InvCovar";      mpKwds[KID_Xform      ] = "Xform";
    mpKwds[KID_GConst      ] = "GConst";        mpKwds[KID_Duration   ] = "Duration";
    mpKwds[KID_InvDiagC    ] = "InvDiagC";      mpKwds[KID_TransP     ] = "TransP";
    mpKwds[KID_DProb       ] = "DProb";         mpKwds[KID_LLTC       ] = "LLTC";
    mpKwds[KID_LLTCovar    ] = "LLTCovar";      mpKwds[KID_XformKind  ] = "XformKind";
    mpKwds[KID_ParentXform ] = "ParentXform";   mpKwds[KID_NumXforms  ] = "NumXforms";
    mpKwds[KID_XformSet    ] = "XformSet";      mpKwds[KID_LinXform   ] = "LinXform";
    mpKwds[KID_Offset      ] = "Offset";        mpKwds[KID_Bias       ] = "Bias";
    mpKwds[KID_BlockInfo   ] = "BlockInfo";     mpKwds[KID_Block      ] = "Block";
    mpKwds[KID_BaseClass   ] = "BaseClass";     mpKwds[KID_Class      ] = "Class";
    mpKwds[KID_XformWgtSet ] = "XformWgtSet";   mpKwds[KID_ClassXform ] = "ClassXform";
    mpKwds[KID_MMFIDMask   ] = "MMFIDMask";     mpKwds[KID_Parameters ] = "Parameters";
    mpKwds[KID_NumClasses  ] = "NumClasses";    mpKwds[KID_AdaptKind  ] = "AdaptKind";
    mpKwds[KID_Prequal     ] = "Prequal";       mpKwds[KID_InputXform ] = "InputXform";
    mpKwds[KID_RClass      ] = "RClass";        mpKwds[KID_RegTree    ] = "RegTree";
    mpKwds[KID_Node        ] = "Node";          mpKwds[KID_TNode      ] = "TNode";
    mpKwds[KID_HMMSetID    ] = "HMMSetID";      mpKwds[KID_ParmKind   ] = "ParmKind";
  
    /* Non-HTK keywords */
    mpKwds[KID_FrmExt      ] = "FrmExt";        mpKwds[KID_PDFObsVec  ] = "PDFObsVec";
    mpKwds[KID_ObsCoef     ] = "ObsCoef";       mpKwds[KID_Input      ] = "Input";
    mpKwds[KID_NumLayers   ] = "NumLayers";     mpKwds[KID_NumBlocks  ] = "NumBlocks";
    mpKwds[KID_Layer       ] = "Layer";         mpKwds[KID_Copy       ] = "Copy";
    mpKwds[KID_Stacking    ] = "Stacking";      mpKwds[KID_Transpose  ] = "Transpose";    
    mpKwds[KID_Constant    ] = "Constant";
    mpKwds[KID_XformPredef ] = "XformPredef";   mpKwds[KID_Window     ] = "Window";  
    mpKwds[KID_WindowPredef] = "WindowPredef";  mpKwds[KID_BlockCopy  ] = "BlockCopy";       
    mpKwds[KID_BiasPredef  ] = "BiasPredef";
        
    mpKwds[KID_ExtendedXform  ] = "ExtendedXform";
    mpKwds[KID_RegionDependent] = "RegionDependentXform";
    mpKwds[KID_GmmPosteriors]   = "GmmPosteriors";
    
    /* Numeric functions - FuncXform*/
    mpKwds[KID_Sigmoid     ] = "Sigmoid";       mpKwds[KID_Log        ] = "Log";
    mpKwds[KID_Exp         ] = "Exp";           mpKwds[KID_Sqrt       ] = "Sqrt";
    mpKwds[KID_SoftMax     ] = "SoftMax";
    
    mpKwds[KID_Weights     ] = "Weights";
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  ModelSet::
  CheckKwd(const char *str, KeywordID kwdID)
  {
    const char *chptr;
    
    if (str[0] == ':') {
      mHmmReadBinary = true;
      return static_cast<unsigned char>(str[1]) == kwdID;
    }
  
    if (str[0] != '<') return 0;
    for (chptr = mpKwds[kwdID], str++; *chptr; chptr++, str++) {
      if (toupper(*chptr) != toupper(*str)) return 0;
    }
  
    if (str[0] != '>') return 0;
  
    assert(str[1] == '\0');
    mHmmReadBinary = false;
    return 1;
  }
  
  //***************************************************************************
  //***************************************************************************
  char *
  ModelSet::
  GetString(FILE *fp, int eofNotExpected)
  {
    char* chptr = mpBuffer;
    int ch;
    int lines = 0;
  
  //  fputs("GetString: ", stdout);
    if (mStringUnget) {
      mStringUnget = false;
  
  //    puts(buffer);
      return mpBuffer;
    }
  
    RemoveSpaces(fp);
  
    ch = getc(fp);
    if (ch == '\"' || ch == '\'' ) {
      int termChar = ch;
  
      while (((ch = getc(fp)) != EOF) && 
            (ch != termChar) && 
            ((chptr-mpBuffer) < static_cast<int>(sizeof(mpBuffer)-1))) 
      {
        if (ch == '\n') {
          ++lines;
        }
        *chptr++ = ch;
      }
  
      if (ch == EOF && ferror(fp)) {
        Error("Cannot read input file %s", mpCurrentMmfName);
      }
  
      if (ch != termChar) {
        Error("Unterminated string constant (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      }
      mCurrentMmfLine += lines;
    } else if (ch == '<') {
      *chptr++ = '<';
      while (((ch = getc(fp)) != EOF) &&
            !isspace(ch) &&
            (ch != '>')  &&
	    (chptr-mpBuffer) < static_cast<int>(sizeof(mpBuffer)-1)) 
      {
        *chptr++ = ch;
      }
  
      if (ch == EOF && ferror(fp)) {
        Error("Cannot read input file %s", mpCurrentMmfName);
      }
  
      if (ch != '>') {
        Error("Unterminated keyword %s (%s:%d)", mpBuffer, mpCurrentMmfName, mCurrentMmfLine);
      }
  
      *chptr++ = '>';
    } else if (ch == ':') {
      *chptr++ = ':';
      *chptr++ = ch = getc(fp);
  
      if (ch == EOF){
      if (ferror(fp)) Error("Cannot read input file %s", mpCurrentMmfName);
      else           Error("Unexpected end of file %s", mpCurrentMmfName);
      }
    } else {
      while ((ch != EOF) && (ch != '\0') &&
            !isspace(ch) &&
            (chptr-mpBuffer) < static_cast<int>(sizeof(mpBuffer)-1)) 
      {
        *chptr++ = ch;
        ch = getc(fp);
      }
  
      if (ch != EOF) {
        if(ch != ' ' || !mHmmReadBinary) {	
          ungetc(ch, fp);
        }
      } else if (ferror(fp)) {
        Error("Cannot read input file %s", mpCurrentMmfName);
      }
  
      if (chptr == mpBuffer) {
        if (eofNotExpected) {
          Error("Unexpected end of file %s", mpCurrentMmfName);
        }
        return NULL;
      }
    }
  
    *chptr = '\0';
  //  puts(buffer);
    return mpBuffer;
  
  }

  //***************************************************************************
  //***************************************************************************  
  void
  ModelSet::
  PutString(FILE *fp, bool binary, const char *msg, ...)
  {
    va_list ap;
    va_start(ap, msg);
    vfprintf(fp, msg, ap);
    if(binary)
    {
      fputc(' ', fp);
    }    
    va_end(ap);    
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  UngetString(void)
  {
    mStringUnget = true;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  ModelSet::
  GetInt(FILE *fp)
  {
    int   cc;
    int   ret;
  //puts("GetInt");
    
    if (mHmmReadBinary) {
      INT_16 i;
      cc = fread(&i, sizeof(INT_16), 1, fp);
      if (!isBigEndian()) swap2(i);
      ret = i;      
    } else {
      RemoveSpaces(fp);
      cc = fscanf(fp, "%d", &ret);
    }
    
    if (cc != 1) {
      if (ferror(fp)) {
        Error("Cannot read input file %s", mpCurrentMmfName);
      }
      Error("Integral number expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  ModelSet::
  GetFloat(FILE *fp)
  {
    int cc;
    FLOAT ret;
  //puts("GetFloat");
    
    if (mHmmReadBinary) {
      FLOAT_32 f;
      cc = fread(&f, sizeof(FLOAT_32), 1, fp);
      ret = f;
      if (!isBigEndian()) swap4(ret);  
    } else {
      RemoveSpaces(fp);
      cc = ReadNumber(fp, &ret);
    }
    
    if (cc != 1) {
      if (ferror(fp)) {
        Error("Cannot read input file %s", mpCurrentMmfName);
      }
      Error("Float number expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  RemoveSpaces(FILE *fp)
  {
    int ch;
  
  //  puts("RemoveSpaces");
    while (isspace(ch = getc(fp))) {
      if (ch == '\n') {
        ++mCurrentMmfLine;
      }
    }
    if (ch != EOF) {
      ungetc(ch, fp);
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  unsigned int 
  faddfloat(FLOAT *vec, size_t size, float mul_const, FILE *fp) 
  {
    size_t    i;
    FLOAT     f;
  
    for (i = 0; i < size; i++) 
    {
      if (fread(&f, sizeof(FLOAT), 1, fp) != 1) break;
      vec[i] += f * mul_const;
    }
    
    return i;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  KeywordID
  ModelSet::
  ReadOutPDFKind(char *str)
  {
    if (     CheckKwd(str, KID_DiagC))    return KID_DiagC;
    else if (CheckKwd(str, KID_InvDiagC)) return KID_InvDiagC;
    else if (CheckKwd(str, KID_FullC))    return KID_FullC;
    else if (CheckKwd(str, KID_XformC))   return KID_XformC;
    else if (CheckKwd(str, KID_LLTC))     return KID_LLTC;
    else if (CheckKwd(str, KID_PDFObsVec))return KID_PDFObsVec;
    else return KID_UNSET;
  }
  
  //***************************************************************************
  //***************************************************************************
  KeywordID
  ModelSet::
  ReadDurKind(char *str)
  {
    if (     CheckKwd(str, KID_NullD))    return KID_NullD;
    else if (CheckKwd(str, KID_PoissonD)) return KID_PoissonD;
    else if (CheckKwd(str, KID_GammaD))   return KID_GammaD;
    else if (CheckKwd(str, KID_GenD))     return KID_GenD;
    else return KID_UNSET;
  }

  
  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  PutKwd(FILE *fp, bool binary, KeywordID kwdID)
  {
    if (binary) {
      putc(':', fp);
      putc(kwdID, fp);
    } else {
      putc('<', fp);
      fputs(mpKwds[kwdID], fp);
      fputs("> ", fp);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  PutInt(FILE *fp, bool binary, int i)
  {
    if (binary) {
      INT_16 b = i;      
      if (!isBigEndian()) swap2(b);
      fwrite(&b, sizeof(INT_16), 1, fp);
    } else {
      fprintf(fp, "%d ", i);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  PutFlt(FILE *fp, bool binary, FLOAT f)
  {
    if (binary) {
      FLOAT_32 b = f;
      if (!isBigEndian()) swap4(b);
      fwrite(&b, sizeof(FLOAT_32), 1, fp);
    } else {
      WriteNumber(fp, f);
      fputc(' ', fp);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  PutNLn(FILE *fp, bool binary)
  {
    if (!binary) putc('\n', fp);
  }
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  PutSpace(FILE *fp, bool binary)
  {
    if (!binary) putc(' ', fp);
  }
  
  //###########################################################################################################
  
  //###########################################################################
  //###########################################################################
  // HMM LIST INPUT
  //###########################################################################
  //###########################################################################
  
  ModelSet::
  ModelSet() : 
    mCurrentMmfLine(1),
    mpCurrentMmfName(0),
    mHmmReadBinary(false),
    mStringUnget(false),
    mHmmsIgnoreMacroRedefinition(true),
    mpHListFilter(0)
  {
  }

  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  ReadHMMList(const char * pFileName, 
              const char * pInMmfDir, 
              const char * pInMmfExt)
  {
    struct ReadlineData       rld = {0};
    char *                    lhmm;
    char *                    fhmm;
    char *                    chptr;
    Macro *                   macro;
    Macro *                   macro2;
    int                       line_no = 0;
    char                      mmfile[1024];
    FILE *                    fp;
  
    if ((fp = my_fopen(pFileName, "rt", mpHListFilter)) == NULL) 
    {
      Error("Cannot open file: '%s'", pFileName);
    }
    
    while ((lhmm = fhmm = readline(fp, &rld)) != NULL) 
    {
      line_no++;
      
      if (getHTKstr(lhmm, &chptr)) 
        Error("%s (%s:%d)", chptr, pFileName, line_no);
      
      if (*chptr && getHTKstr(fhmm = chptr, &chptr)) 
        Error("%s (%s:%d)", chptr, pFileName, line_no);
      
      if ((macro = FindMacro(&mHmmHash, fhmm)) == NULL) 
      {
        mmfile[0] = '\0';
        
        if (pInMmfDir) 
          strcat(strcat(mmfile, pInMmfDir), "/");
        
        strcat(mmfile, fhmm);
        
        if (pInMmfExt) 
          strcat(strcat(mmfile, "."), pInMmfExt);
        
        ParseMmf(mmfile, fhmm, false);
        
        if ((macro = FindMacro(&mHmmHash, fhmm)) == NULL)
          Error("Definition of model '%s' not found in file '%s'", fhmm, mmfile);
      }
      
      if (lhmm != fhmm) 
      {
        mpCurrentMmfName = NULL; // Global variable; macro will not be written to any output MMF
        macro2 = this->pAddMacro('h', lhmm);
        assert(macro2 != NULL);
        
        if (macro2->mpData != NULL) 
        {
          if (mHmmsIgnoreMacroRedefinition == 0) 
          {
            Error("Redefinition of HMM %s (%s:%d)", lhmm, pFileName, line_no);
          } 
          else 
          {
            Warning("Redefinition of HMM %s (%s:%d) is ignored",
                    lhmm, pFileName, line_no);
          }
        } 
        else 
        {
          macro2->mpData = macro->mpData;
        }
      }
    }
    if (ferror(fp) || my_fclose(fp))
      Error("Cannot read HMM list file %s", pFileName);
  }

  //###########################################################################################################
  
  //###########################################################################
  //###########################################################################
  // MMF INPUT
  //###########################################################################
  //###########################################################################
  
  //***************************************************************************
  //***************************************************************************
  void
  ModelSet::
  ParseMmf(const char * pFileName, char * expectHMM, bool readOnly)
  {
    IStkStream    input_stream;
    FILE*         fp;
    char*         keyword;
    Macro*        macro;
    MacroData*    data = NULL;
    std::string   filter;
    
    mCurrentMmfLine = 1;
    mpCurrentMmfName = pFileName;
    
    try
    {
      filter = (NULL != gpMmfFilter) ? gpMmfFilter : "";
      
      // try to open the stream
      input_stream.open(pFileName, std::ios::in|std::ios::binary, gpMmfFilter);
      if (!input_stream.good())
      {
        Error("Cannot open input MMF %s", pFileName);
      }
      fp = input_stream.file();    
      
      for (;;) 
      {
        if ((keyword = GetString(fp, 0)) == NULL) 
        {
          if (ferror(fp)) 
          {
            Error("Cannot read input MMF", pFileName);
          }
          input_stream.close();
          return;
        }
    
        
        if (keyword[0] == '~' && keyword[2] == '\0' ) 
        {
          int type = keyword[1];
    
          if (type == 'o') {
            if (!ReadGlobalOptions(fp, readOnly)) {
              Error("No global option defined (%s:%d)", pFileName, mCurrentMmfLine);
            }
          } else {
            keyword = GetString(fp, 1);
            if ((macro = pAddMacro(type, keyword)) == NULL) {
              Error("Unrecognized macro type ~%c (%s:%d)", type, pFileName, mCurrentMmfLine);
            }
            
            // this is perhaps the place, where we say whether we want 
            // to write the macro in future
            macro->mReadOnly = readOnly;
    
            if (macro->mpData != NULL) {
              if (mHmmsIgnoreMacroRedefinition == 0) {
                Error("Redefinition of macro ~%c %s (%s:%d)", type, keyword, pFileName, mCurrentMmfLine);
              } else {
                Warning("Redefinition of macro ~%c %s (%s:%d)", type, keyword, pFileName, mCurrentMmfLine);
    //            Warning("Redefinition of macro ~%c %s (%s:%d) is ignored", type, keyword, mmFileName, mCurrentMmfLine);
              }
            }
            switch (type)
            {
              case 'h': data =  ReadHMM          (fp, macro); break;
              case 's': data =  ReadState        (fp, macro); break;
              case 'm': data =  ReadMixture      (fp, macro); break;
              case 'u': data =  ReadMean         (fp, macro); break;
              case 'v': data =  ReadVariance     (fp, macro); break;
              case 't': data =  ReadTransition   (fp, macro); break;
              case 'j': data =  ReadXformInstance(fp, macro); break;
              case 'x': data =  ReadXform        (fp, macro); break;
              default : data =  NULL; break;
            }  
    
            assert(data != NULL);
            
            if (macro->mpData == NULL) {
              macro->mpData = data;
            } else {
    
              // Macro is redefined. New item must be checked for compatibility with
              // the old one (vector sizes, delays, memory sizes) All references to
              // old item must be replaced and old item must be released
              // !!! How about AllocateAccumulatorsForXformStats() and ResetAccumsForHMMSet()
              // !!! Replacing HMM will not work correctly after attaching network nodes to models    
    
              unsigned int i;
              MyHSearchData *hash = NULL;
              ReplaceItemUserData ud;
              ud.mpOldData = macro->mpData;
              ud.mpNewData = data;
              ud.mType     = type;
              
              switch (type) {
              case 'h':
                ud.mpOldData->Scan(MTM_REVERSE_PASS | MTM_ALL, NULL,ReleaseItem,NULL);
                hash = &mHmmHash;
                break;
              case 's':
                this->Scan(MTM_HMM, NULL,ReplaceItem, &ud);
                ud.mpOldData->Scan(MTM_REVERSE_PASS | MTM_ALL, NULL,ReleaseItem,NULL);
                hash = &mStateHash;
                break;
              case 'm':
                this->Scan(MTM_STATE, NULL,ReplaceItem, &ud);
                ud.mpOldData->Scan(MTM_REVERSE_PASS | MTM_ALL, NULL,ReleaseItem,NULL);
                hash = &mMixtureHash;
                break;
              case 'u':
                this->Scan(MTM_MIXTURE, NULL,ReplaceItem, &ud);
                delete ud.mpOldData;
                hash = &mMeanHash;
                break;
              case 'v':
                this->Scan(MTM_MIXTURE, NULL,ReplaceItem, &ud);
                delete ud.mpOldData;
                hash = &mVarianceHash;
                break;
              case 't':
                this->Scan(MTM_HMM, NULL,ReplaceItem, &ud);
                delete ud.mpOldData;
                hash = &mTransitionHash;
                break;
              case 'j':
                this->Scan(MTM_XFORM_INSTANCE | MTM_MIXTURE, NULL, ReplaceItem, &ud);
                ud.mpOldData->Scan(MTM_REVERSE_PASS|MTM_ALL, NULL, ReleaseItem, NULL);
                hash = &mXformInstanceHash;
                break;
              case 'x':
                //:TODO:
                // Fix this
                if (mUpdateMask & UM_CWEIGHTS
                &&  XT_BIAS == static_cast<Xform*>(ud.mpNewData)->mXformType)
                {
                  for (int i = 0; i < mNClusterWeightVectors; i++)
                  {
                    if (mpClusterWeightVectors[i] == static_cast<BiasXform*>(ud.mpOldData))
                    {
                      // Recompute the weights vector
                      ComputeClusterWeightsVector(i);
                      WriteClusterWeightsVector(i);
                      
                      // clear the accumulators                          
                      mpGw[i].Clear();
                      mpKw[i].Clear();
                      mpClusterWeightVectors[i] = static_cast<BiasXform*>(ud.mpNewData);
                      
                      free(ud.mpNewData->mpMacro->mpFileName);
                      if ((ud.mpNewData->mpMacro->mpFileName = strdup(mpCurrentMmfName)) == NULL)
                      {
                        Error("Insufficient memory");
                      }
                    }
                  }
                }              
                
                this->Scan(MTM_XFORM | MTM_XFORM_INSTANCE | MTM_MEAN | MTM_MIXTURE ,NULL,ReplaceItem, &ud);
                ud.mpOldData->Scan(MTM_REVERSE_PASS | MTM_ALL, NULL,ReleaseItem,NULL);
                hash = &mXformHash;
                break;
              }
              
              for (i = 0; i < hash->mNEntries; i++) {
                if (reinterpret_cast<Macro *>(hash->mpEntry[i]->data)->mpData == ud.mpOldData) {
                  reinterpret_cast<Macro *>(hash->mpEntry[i]->data)->mpData = ud.mpNewData;
                }
              }
  /*            for (macro = mpFirstMacro; macro != NULL; macro->nextAll)
              {
                macro
              }*/
            }
          }
        } else if (CheckKwd(keyword, KID_BeginHMM)) {
          UngetString();
          if (expectHMM == NULL) {
            Error("Macro definition expected (%s:%d)",pFileName,mCurrentMmfLine);
          }
          //macro = AddMacroToHMMSet('h', expectHMM, hmm_set);
          macro = pAddMacro('h', expectHMM);
          
          // this is perhaps the place, where we say whether we want 
          // to write the macro in future
          macro->mReadOnly = readOnly;
    
          macro->mpData = ReadHMM(fp, macro, readOnly);
        } else {
          Error("Unexpected keyword %s (%s:%d)",keyword,pFileName,mCurrentMmfLine);
        }
      }
      
      if (input_stream.is_open())
      {
        input_stream.close();
      }
      
    } // try
    
    catch (...)
    {
      if (input_stream.is_open())
      {
        input_stream.close();
      }

      // if cluster weights update
      if ((mUpdateMask & UM_CWEIGHTS) && mClusterWeightsStream.is_open())
      {
        mClusterWeightsStream.close();
      }
      throw;
    }
  }
  

  //***************************************************************************
  //***************************************************************************
  Hmm *
  ModelSet::
  ReadHMM(FILE * fp, Macro * macro, bool readOnly)
  {
    Hmm *   ret;
    char *  keyword;
    size_t  nstates;
    size_t  i;
    int     state_id;
  
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~h")) 
    {    
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mHmmHash, keyword)) == NULL) 
      {
        Error("Undefined reference to macro ~h %s (%s:%d)",
        keyword, mpCurrentMmfName, mCurrentMmfLine);
      }
      return (Hmm *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_BeginHMM))
      Error("Keyword <BeginHMM> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    ReadGlobalOptions(fp, readOnly);
  
    if (mInputVectorSize == -1)
      Error("<VecSize> is not defined yet (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    if (mDurKind == -1) 
      mDurKind = KID_NullD;
  
    keyword = GetString(fp, 1);
    
    
    if (!CheckKwd(keyword, KID_NumStates))
      Error("Keyword <NumStates> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    nstates = GetInt(fp);
    
    if (nstates < 3)
      Error("HMM must have at least 3 states (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      
    ret = new Hmm(nstates);  
    
  
    for (i=0; i<nstates-2; i++) 
      ret->mpState[i] = NULL;
  
    for (i=0; i<nstates-2; i++) 
    {
      keyword = GetString(fp, 1);
      
      if (!CheckKwd(keyword, KID_State)) 
        Error("Keyword <State> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
      state_id = GetInt(fp);
  
  //    printf("%d\n", state_id);
      if (state_id < 2 || state_id >= static_cast<int>(nstates)) 
        Error("State number out of the range (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
      if (ret->mpState[state_id-2] != NULL) 
        Error("Redefinition of state (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
      ret->mpState[state_id-2] = ReadState(fp, NULL);
  //    printf("\n%d: %x\n", state_id-2, ret->mpState[state_id-2]);
    }
  
    ret->mpTransition    = ReadTransition(fp, NULL);
    ret->mpMacro         = macro;
  
    if (ret->mpTransition->mNStates != nstates) 
      Error("Invalid transition matrix size (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    keyword = GetString(fp, 1);
    
    if (!CheckKwd(keyword, KID_EndHMM))
      Error("Keyword <EndHMM> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    
    return ret;
  }

  
  //***************************************************************************
  //***************************************************************************
  State*
  ModelSet::
  ReadState(FILE* fp, Macro* macro)
  {
    State *     ret;
    char *      keyword;
    int         mixture_id;
    int         i;
    int         num_mixes = 1;
    FLOAT       mixture_weight;
  
  //  puts("ReadState");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~s")) 
    {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mStateHash, keyword)) == NULL) 
      {
        Error("Undefined reference to macro ~s %s (%s:%d)", keyword, mpCurrentMmfName, mCurrentMmfLine);
      }
      return (State *) macro->mpData;
    }
  
    if (mOutPdfKind == -1) 
    {
      mOutPdfKind = CheckKwd(keyword, KID_ObsCoef)
                            ? KID_PDFObsVec : KID_DiagC;
    }
  
    if (mOutPdfKind == KID_PDFObsVec) 
    {
      num_mixes = 0;
    } 
    else if (CheckKwd(keyword, KID_NumMixes)) 
    {
      num_mixes = GetInt(fp);
  
      keyword = GetString(fp, 1);
    }
  
    // ***
    // old malloc
    // ret = (State*) malloc(sizeof(State) + (num_mixes-1)*sizeof(ret->mpMixture[0]));
    //
    // if (ret == NULL) 
    //  Error("Insufficient memory");
    ret = new State(num_mixes);
  
    if (mOutPdfKind == KID_PDFObsVec) 
    {
      int range;
  
      if (!CheckKwd(keyword, KID_ObsCoef)) 
      {
        Error("Keyword <ObsCoef> expected (%s:%d)",
              mpCurrentMmfName, mCurrentMmfLine);
      }
      ret->PDF_obs_coef = GetInt(fp) - 1;
      range = mpInputXform ? mpInputXform->OutSize()
                           : mInputVectorSize;
      if (ret->PDF_obs_coef < 0 || ret->PDF_obs_coef >= range) 
      {
        Error("Parameter <ObsCoef> is out of the range 1:%d (%s:%d)",
              range, mpCurrentMmfName, mCurrentMmfLine);
      }
    } 
    else 
    {
      ret->mNMixtures = num_mixes;
    //  printf("ptr: %x num_mixes: %d\n", ret, num_mixes);
  
      for (i=0; i<num_mixes; i++) 
        ret->mpMixture[i].mpEstimates = NULL;
  
      if (CheckKwd(keyword, KID_Stream)) 
      {
        if (GetInt(fp) != 1) 
        {
          Error("Stream number out of the range (%s:%d)",
                mpCurrentMmfName, mCurrentMmfLine);
        }
      } 
      else 
      {
        UngetString();
      }
  
      for (i=0; i<num_mixes; i++) 
      {
        keyword = GetString(fp, 1);
        if (!CheckKwd(keyword, KID_Mixture)) 
        {
          if (num_mixes > 1) 
          {
            Error("Keyword <Mixture> expected (%s:%d)",
                  mpCurrentMmfName, mCurrentMmfLine);
          }
          UngetString();
          mixture_id = 1;
          mixture_weight = 1.0;
        } 
        else 
        {
          mixture_id = GetInt(fp);
          mixture_weight = GetFloat(fp);
        }
  
        if (mixture_id < 1 || mixture_id > num_mixes) 
        {
          Error("Mixture number out of the range (%s:%d)",
                mpCurrentMmfName, mCurrentMmfLine);
        }
  
        if (ret->mpMixture[mixture_id-1].mpEstimates != NULL) 
        {
          Error("Redefinition of mixture %d (%s:%d)",
                mixture_id, mpCurrentMmfName, mCurrentMmfLine);
        }
  
        ret->mpMixture[mixture_id-1].mpEstimates  = ReadMixture(fp, NULL);
        ret->mpMixture[mixture_id-1].mWeight      = my_log(mixture_weight);
        ret->mpMixture[mixture_id-1].mWeightAccum = 0.0;
      }
    }
    
    ret->mOutPdfKind = mOutPdfKind;
    ret->mID = mNStates;
    mNStates++;
    ret->mpMacro = macro;
  
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Mixture*
  ModelSet::
  ReadMixture(FILE *fp, Macro *macro)
  {
    Mixture* ret = NULL;
    char*    keyword;
  
  //  puts("ReadMixture");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~m")) 
    {
      keyword = GetString(fp, 1);
      
      if ((macro = FindMacro(&mMixtureHash, keyword)) == NULL)
        Error("Undefined reference to macro ~m %s (%s:%d)",
              keyword, mpCurrentMmfName, mCurrentMmfLine);
      
      return (Mixture *) macro->mpData;
    }
  
    // create new mixture object
    ret = new Mixture;  
  
    if (CheckKwd(keyword, KID_InputXform)) 
    {
      ret->mpInputXform = ReadXformInstance(fp, NULL);
    } 
    else 
    {
      ret->mpInputXform = this->mpInputXform;
      UngetString();
    }
  
    //
    ret->mpMean = ReadMean(fp, NULL);
    
    // allocate accumulators for Cluster Adaptive Training cluster parameters
    if (NULL != ret->mpMean->mpClusterWeightVectors)
    {
      // Cluster Parameter Update section
      ret->mAccumG.Init(ret->mpMean->mClusterMatrixT.Rows(),
                        ret->mpMean->mClusterMatrixT.Rows());
      ret->mAccumK.Init(ret->mpMean->mClusterMatrixT.Rows(),
                        ret->mpMean->mClusterMatrixT.Cols());
      ret->mAccumL.Init(ret->mpMean->mClusterMatrixT.Cols());
      
      // per-speaker partial accumulators
      ret->mPartialAccumG = 0.0;
      ret->mPartialAccumK.Init(ret->mpMean->mClusterMatrixT.Cols());
      
      // for discriminative training
      if (UT_EBW == mUpdateType)
      {
        ret->mPartialAccumGd = 0.0;
        ret->mAccumGd.Init(ret->mpMean->mClusterMatrixT.Rows(),
                           ret->mpMean->mClusterMatrixT.Rows());        
      }
    }
    
    //
    ret->mpVariance = ReadVariance(fp, NULL);
  
    if (!ret->mpInputXform && (mInputVectorSize == -1))
    {
      Error("<VecSize> is not defined yet (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    } 
    
    size_t size = ret->mpInputXform ? ret->mpInputXform->OutSize() : mInputVectorSize;
    
    if (ret->mpMean->VectorSize()     != size || 
        ret->mpVariance->VectorSize() != size) 
    {
      Error("Invalid mean or variance vector size (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    if ((keyword = GetString(fp, 0)) != NULL && CheckKwd(keyword, KID_GConst)) 
    {
      ret->rGConst() = GetFloat(fp);
    } 
    else 
    {
      ret->ComputeGConst();
      if (keyword != NULL) UngetString();
    }
  
    ret->mID = mNMixtures;
    mNMixtures++;
    ret->mpMacro = macro;
  
    return ret;
  }


  //***************************************************************************
  //***************************************************************************
  Mean*
  ModelSet::
  ReadMean(FILE* fp, Macro* macro)
  {
    Mean*     ret = NULL;
    char*     keyword;
    int       vec_size;
    int       i;
//    int       accum_size = 0;
    
  //  puts("ReadMean");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~u")) 
    {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mMeanHash, keyword)) == NULL) 
      {
        Error("Undefined reference to macro ~u %s (%s:%d)", 
          keyword, 
          mpCurrentMmfName, 
          mCurrentMmfLine);
      }
      return  (Mean*) macro->mpData;
    }
  
    // Mean definition can be defined by a set of weighted means
    if (CheckKwd(keyword, KID_Weights))
    {
      Macro*        tmp_macro;
      FLOAT         tmp_val;  // temporary read value
      int           n_xforms; // number of xform refferences
      BiasXform**   xforms = NULL; // temporary storage of xform refferences
      int           total_means = 0; 
      bool          init_Gwkw = false;
      ret = NULL;
      
      // now we expect number of Xform references
      n_xforms = GetInt(fp);
      // allocate space for refferences to these xforms
      xforms = new BiasXform*[n_xforms];
      
      // we need to know globally which bias xforms are the weights
      if ((mUpdateMask & UM_CWEIGHTS)
      && (NULL == mpClusterWeightVectors))
      {
        mNClusterWeightVectors = n_xforms;      
        mpClusterWeightVectors = new BiasXform*[n_xforms];
        mpGw             = new Matrix<FLOAT>[n_xforms];
        mpKw             = new BasicVector<FLOAT>[n_xforms];
        init_Gwkw        = true;
      }
      
      
      for (int xform_i = 0; xform_i < n_xforms; xform_i++)
      {
        // expect Xform macro reference
        keyword = GetString(fp, 1);
        if (!strcmp(keyword, "~x"))
        {
          keyword = GetString(fp, 1);
          if (NULL == (tmp_macro = (FindMacro(&mXformHash, keyword))))
          {
            Error("Undefined reference to macro ~x %s (%s:%d)", 
              keyword, 
              mpCurrentMmfName, 
              mCurrentMmfLine);
          }
          
          xforms[xform_i] = static_cast<BiasXform *>(tmp_macro->mpData);
          total_means    += xforms[xform_i]->mInSize;
          
          if (mUpdateMask & UM_CWEIGHTS)
          {
            mpClusterWeightVectors[xform_i] = xforms[xform_i];
          }          
        }
        else
        {
          Error("Reference to macro ~x expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        }
      } //for (int xform_i; xform_i < n_xforms; xform_i++)
        
      // at this moment we know how many means we are going to use
      if (init_Gwkw)
      {
        for (int xform_i = 0; xform_i < n_xforms; xform_i++)
        {
          mpGw[xform_i].Init(total_means, total_means);
          mpKw[xform_i].Init(total_means);
        }
      }
      
      // read total_means means
      for (int i = 0; i < total_means; i++)
      {
        keyword = GetString(fp, 1);
        if (!CheckKwd(keyword, KID_Mean))
          Error("Keyword <Mean> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      
        vec_size = GetInt(fp);
        // maybe we haven't created any Mean object yet
        if (NULL == ret)
        {
          // create new Mean object
          ret = new Mean(vec_size, mAllocAccums);
          // initialize matrix for original mean vectors definition
          ret->mClusterMatrixT.Init(total_means, vec_size);
          // set the weights vectors
          ret->mpClusterWeightVectors = xforms;
          ret->mNClusterWeightVectors = n_xforms;
          
          if (mUpdateMask & UM_CWEIGHTS)
          {
            ret->mCwvAccum.Init(n_xforms, vec_size);
            ret->mpOccProbAccums = new FLOAT[n_xforms];
            memset(ret->mpOccProbAccums, 0, n_xforms * sizeof(FLOAT));
          }
        }
        
        // read the values and fill appropriate structures
        for (int j = 0; j < vec_size; j++)
        {
          tmp_val = GetFloat(fp);
          // fill the matrix 
          ret->mClusterMatrixT[i][j] = tmp_val;
        }
        
      } // for (i = 0; i < bx->mInSize; i++)
      
                  

      // recalculate the real vector using cluster mean vectors
      ret->RecalculateCAT();
    } // if (CheckKwd(keyword, KID_Weights))
    
    else if (CheckKwd(keyword, KID_Mean))
    {
      vec_size = GetInt(fp);
          
      ret = new Mean(vec_size, mAllocAccums);
  
      for (i=0; i<vec_size; i++) 
      {
        // old vector
        //ret->mpVectorO[i] = GetFloat(fp);
        ret->mVector[i] = GetFloat(fp);
      }      
    } // else if (CheckKwd(keyword, KID_Mean))
    
    else
    {
      Error("Keyword <Mean> or <Weights> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
    
    ret->mpMacro = macro;
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Variance *
  ModelSet::
  ReadVariance(FILE *fp, Macro *macro)
  {
    Variance *ret;
    char *keyword;
    int vec_size, i;
    //int accum_size = 0;
  
  //  puts("ReadVariance");
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~v")) {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mVarianceHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~v %s (%s:%d)", keyword, mpCurrentMmfName, mCurrentMmfLine);
  
      }
      return (Variance *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_Variance)) {
      Error("Keyword <Variance> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    vec_size = GetInt(fp);
  
    ret = new Variance(vec_size, mAllocAccums);
  
    for (i=0; i < vec_size; i++) 
    {
      ret->mVector[i] = 1.0 / GetFloat(fp);
    }
  
    ret->mpMacro = macro;
    return ret;
  }
    
  
  //***************************************************************************
  //***************************************************************************
  XformInstance *
  ModelSet::
  ReadXformInstance(FILE *fp, Macro *macro) 
  {
    XformInstance * ret;
    XformInstance * input = NULL;
    char *          keyword;
    int             out_vec_size;
    int             i;
  
  //  puts("ReadXformInstance");
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~j")) {
      keyword = GetString(fp, 1);
      
      if ((macro = FindMacro(&mXformInstanceHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~j %s (%s:%d)", keyword, 
            mpCurrentMmfName, mCurrentMmfLine);
      }
      
      return (XformInstance *) macro->mpData;
    }
  
    if (CheckKwd(keyword, KID_Input)) 
    {
      input = ReadXformInstance(fp, NULL);
      keyword = GetString(fp, 1);
    }
  
    if (CheckKwd(keyword, KID_MMFIDMask)) 
    {
      keyword = GetString(fp, 1);
      if (strcmp(keyword, "*")) {
        Error("<MMFIdMask> different than '*' is not supported (%s:%d)", 
            mpCurrentMmfName, mCurrentMmfLine);
      }
  
      keyword = GetString(fp, 1);
    }
  
    if ((i = ReadParmKind(keyword, true)) != -1) {
      if (mParamKind != -1 && mParamKind != i) {
        Error("ParamKind mismatch (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      }
  
      keyword = GetString(fp, 1);
    }
  
    if (CheckKwd(keyword, KID_LinXform)) {
      keyword = GetString(fp, 1);
  //    Error("Keyword <LinXform> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    if (!CheckKwd(keyword, KID_VecSize)) {
      Error("Keyword <VecSize> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    out_vec_size = GetInt(fp);
  
    // create new instance
    ret = new XformInstance(out_vec_size);

    // read instance from file    
    ret->mpXform = ReadXform(fp, NULL);
    
    
    if (input == NULL && mInputVectorSize == -1) {
      Error("<VecSize> has not been defined yet (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    if (static_cast<size_t>(out_vec_size) != ret->mpXform->mOutSize) {
      Error("XformInstance <VecSize> must equal to Xform "
            "output size (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    if (input == NULL) {
      if (ret->mpXform->mInSize != static_cast<size_t>(mInputVectorSize)) {
        Error("Xform input size must equal to ~o <VecSize> (%s:%d)",
              mpCurrentMmfName, mCurrentMmfLine);
      }
    } else {
      if (ret->mpXform->mInSize != input->OutSize()) {
        Error("Xform input size must equal to <Input> <VecSize> (%s:%d)",
              mpCurrentMmfName, mCurrentMmfLine);
      }
    }
  
    ret->mpMemory = NULL;
    
    if (ret->mpXform->mMemorySize > 0)
      ret->mpMemory = new char[ret->mpXform->mMemorySize];
    
    memset(ret->mpMemory, 0, ret->mpXform->mMemorySize);
  
    ret->mpNext = mpXformInstances;
    mpXformInstances = ret;
    ret->mpInput = input;
    ret->mpMacro = macro;
  //  puts("ReadXformInstance exit");
  
    ret->mNumberOfXformStatCaches = 0;
    ret->mpXformStatCache   = NULL;
    ret->mTotalDelay = ret->mpXform->mDelay + (input ? input->mTotalDelay : 0);
  
    mTotalDelay = HIGHER_OF(mTotalDelay, ret->mTotalDelay);
    return ret;
  }; //ReadXformInstance(FILE *fp, Macro *macro) 
  
  

  //****************************************************************************
  //****************************************************************************
  Xform *
  ModelSet::
  ReadXform(FILE *fp, Macro *macro) 
  {
    char *keyword;
    unsigned int i;
  
  //  puts("ReadXform");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~x")) {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mXformHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~x %s (%s:%d)", keyword, 
            mpCurrentMmfName, mCurrentMmfLine);
      }
      return (Xform *) macro->mpData;
    }
    if (CheckKwd(keyword, KID_Xform)) {
      return (Xform *) ReadLinearXform(fp, macro, false);
    }

    if (CheckKwd(keyword, KID_XformPredef)) {
      return (Xform *) ReadLinearXform(fp, macro, true);
    }    
      
    if (CheckKwd(keyword, KID_Bias)) {
      return (Xform *) ReadBiasXform(fp, macro, false);
    }

    if (CheckKwd(keyword, KID_BiasPredef)) {
      return (Xform *) ReadBiasXform(fp, macro, true);
    }    
      
    if (CheckKwd(keyword, KID_Copy)) {
      return (Xform *) ReadCopyXform(fp, macro);
    }

    if (CheckKwd(keyword, KID_BlockCopy)) {
      return (Xform *) ReadBlockCopyXform(fp, macro);
    }    
    
    if (CheckKwd(keyword, KID_Transpose)) {
      return (Xform *) ReadTransposeXform(fp, macro);
    }

    if (CheckKwd(keyword, KID_Window)) {
      return (Xform *) ReadWindowXform(fp, macro, false);
    }
    
    if (CheckKwd(keyword, KID_WindowPredef)) {
      return (Xform *) ReadWindowXform(fp, macro, true);
    }
    
    if (CheckKwd(keyword, KID_Stacking)) {
      return (Xform *) ReadStackingXform(fp, macro);
    }
    if (CheckKwd(keyword, KID_Constant)) {
      return (Xform *) ReadConstantXform(fp, macro);
    }
    if (CheckKwd(keyword, KID_GmmPosteriors)) {
      return (Xform *) ReadGmmPosteriorsXform(fp, macro);
    }
    if (CheckKwd(keyword, KID_RegionDependent)) {
      return (Xform *) ReadRegionDependentXform(fp, macro);
    }
    if (CheckKwd(keyword, KID_NumLayers) ||
      CheckKwd(keyword, KID_NumBlocks) ||
      CheckKwd(keyword, KID_BlockInfo))
    {
      UngetString();
      return (Xform *) ReadCompositeXform(fp, macro);
    }
  
    for (i=0; i < gFuncTableSize; i++) 
    {
      if (CheckKwd(keyword, gFuncTable[i].KID)) {
        return (Xform *) ReadFuncXform(fp, macro, i);
      }
    }
  
    // check for XForm extension
    if (CheckKwd(keyword, KID_ExtendedXform))
    {
      // what extension???
      keyword = GetString(fp, 1);
      
      if (!strcmp(keyword, "FeatureMapping"))
      {
        return (Xform *) ReadFeatureMappingXform(fp, macro);
      }
      
      if (!strcmp(keyword, "FrantaProduct"))
      {
        return static_cast<Xform*>(ReadFrantaProductXform(fp, macro));
      }
      
      if (!strcmp(keyword, "Matlab"))
      {
        return static_cast<Xform*>(ReadMatlabXform(fp, macro));
      }
      
      Error("Invalid Extended Xform definition; unexpected keyword %s (%s:%d)", 
        keyword, mpCurrentMmfName, mCurrentMmfLine);
    }
    
    Error("Invalid Xform definition; unexpected keyword %s (%s:%d)", 
      keyword, mpCurrentMmfName, mCurrentMmfLine);
      
    return NULL;
  }; //ReadXform(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************
  //***************************************************************************
  CompositeXform *
  ModelSet::
  ReadCompositeXform(FILE *fp, Macro *macro) 
  {
    CompositeXform *    ret;
    Xform **            block;
    char *              keyword;
    size_t              i;
    size_t              j;
    int                 layer_delay;
    int                 layer_id;
    size_t              nlayers;
    int                 block_id;
    size_t              nblocks;
    size_t              prev_out_size = 0;
  
    keyword = GetString(fp, 1);
    if (CheckKwd(keyword, KID_NumLayers)) {
      nlayers = GetInt(fp);
    } else {
      nlayers = 1;
      UngetString();
    }
  
    //if ((ret = (CompositeXform *) malloc(sizeof(CompositeXform) +
    //                            (nlayers-1) * sizeof(ret->mpLayer[0]))) == NULL) {
    //  Error("Insufficient memory");
    //}
    //
    //ret->mMemorySize = 0;
    //ret->mDelay      = 0;
    //ret->mNLayers    = nlayers;
    //
    //for (i=0; i<nlayers; i++) 
    //  ret->mpLayer[i].mpBlock = NULL;
  
    // create new Composite Xform
    ret = new CompositeXform(nlayers);    
    
    // load each layer
    for (i=0; i<nlayers; i++) 
    {
      keyword = GetString(fp, 1);
      
      if (!CheckKwd(keyword, KID_Layer)) 
      {
        if (nlayers > 1) 
          Error("Keyword <Layer> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        
        layer_id = 1;
      } 
      else 
      {
        layer_id = GetInt(fp);
        keyword = GetString(fp, 1);
      }
  
      if (layer_id < 1 || static_cast<size_t>(layer_id) > nlayers)
        Error("Layer number out of the range (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      
      if (ret->mpLayer[layer_id-1].mpBlock != NULL)
        Error("Redefinition of mpLayer (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
      if (CheckKwd(keyword, KID_NumBlocks)) 
      {
        nblocks = GetInt(fp);
      }  
      else if (CheckKwd(keyword, KID_BlockInfo)) 
      {
        nblocks = GetInt(fp);
        for (j = 0; j < nblocks; j++) 
          GetInt(fp); //Blocks' output sizes are not needed
      }  
      else 
      {
        nblocks = 1;
        UngetString();
      }
  
      //***
      // old malloc
      //       if ((block = (Xform **)  malloc(sizeof(Xform*) * nblocks)) == NULL) 
      //         Error("Insufficient memory");
      //   
      //       ret->mpLayer[layer_id-1].mpBlock   = block;
      //       ret->mpLayer[layer_id-1].mNBlocks = nblocks;
      //       
      //       for (j = 0; j < nblocks; j++) 
      //         block[j] = NULL;
      //   
      
      // initialize the blocks within the layer (returns pointer to the begining
      // of block array
      block = ret->mpLayer[layer_id-1].InitBlocks(nblocks);
      
      layer_delay = 0;
      
      for (j = 0; j < nblocks; j++) 
      {
        keyword = GetString(fp, 1);
        
        if (!CheckKwd(keyword, KID_Block)) 
        {
          if (nblocks > 1) {
            Error("Keyword <Block> expected (%s:%d)", mpCurrentMmfName, 
                mCurrentMmfLine);
          }
          
          UngetString();
          block_id = 1;
        } 
        else 
        {
          block_id = GetInt(fp);
        }
  
        if (block_id < 1 || static_cast<size_t>(block_id) > nblocks) {
          Error("Block number out of the range (%s:%d)", mpCurrentMmfName, 
              mCurrentMmfLine);
        }
  
        if (block[block_id-1] != NULL)
          Error("Redefinition of block (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
        block[block_id-1] = ReadXform(fp, NULL);
        ret->mMemorySize += block[block_id-1]->mMemorySize;
        layer_delay = HIGHER_OF(layer_delay ,block[block_id-1]->mDelay);
      }
      ret->mDelay += layer_delay;
    }
  
    for (i=0; i<nlayers; i++) 
    {
      size_t    layer_in_size  = 0;
      size_t    layer_out_size = 0;
  
      for (j=0; j < ret->mpLayer[i].mNBlocks; j++) 
      {
        layer_in_size  += ret->mpLayer[i].mpBlock[j]->mInSize;
        layer_out_size += ret->mpLayer[i].mpBlock[j]->mOutSize;
      }
  
      if (i == nlayers-1) 
        ret->mOutSize = layer_out_size;
      
      if (i == 0)
      {
        ret->mInSize  = layer_in_size;
      }
      else 
      {
        if (prev_out_size < layer_in_size) 
        {
          Error("Output size of mpLayer %d (%d) is smaller then input size of mpLayer %d (%d) (%s:%d)",
                (int) i, (int) prev_out_size, (int) i+1, (int) layer_in_size, mpCurrentMmfName, mCurrentMmfLine);
        }
  
        //***
        // old malloc
        //
        //if ((ret->mpLayer[i-1].mpOutputVector = (FLOAT *) malloc(prev_out_size * sizeof(FLOAT))) == NULL) 
        //{
        //  Error("Insufficient memory");
        //}
        
        //ret->mpLayer[i-1].mpOutputVector = new FLOAT[prev_out_size];
        ret->mpLayer[i-1].mOutputVector.Init(prev_out_size);
      }
  
      prev_out_size = layer_out_size;
    }
  
    
    ret->mpMacro = macro;
    return ret;
  }; // ReadCompositeXform(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************
  //***************************************************************************
  RegionDependentXform *
  ModelSet::
  ReadRegionDependentXform(FILE *fp, Macro *macro) 
  {
    RegionDependentXform* ret;
    char *              keyword;
    size_t              j;
    int                 block_id;
    size_t              nblocks;

    keyword = GetString(fp, 1);
    if (CheckKwd(keyword, KID_NumBlocks)) {
      nblocks = GetInt(fp);
    } else {
      nblocks = 1;
      UngetString();
    }
  
    // create new RegionDependent Xform
    ret = new RegionDependentXform(nblocks);    

    for (j = 0; j < nblocks; j++) {
      keyword = GetString(fp, 1);
      
      if (!CheckKwd(keyword, KID_Block)) {
        if (nblocks > 1) {
          Error("Keyword <Block> expected (%s:%d)", mpCurrentMmfName, 
              mCurrentMmfLine);
        }
        UngetString();
        block_id = 1;
      } else {
        block_id = GetInt(fp);
      }
      if (block_id < 1 || static_cast<size_t>(block_id) > nblocks) {
        Error("Block number out of the range (%s:%d)", mpCurrentMmfName, 
            mCurrentMmfLine);
      }
      if (ret->mpBlock[block_id-1] != NULL) {
        Error("Redefinition of block (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      }
      ret->mpBlock[block_id-1] = ReadXform(fp, NULL);
      ret->mMemorySize += ret->mpBlock[block_id-1]->mMemorySize;
      ret->mDelay = HIGHER_OF(ret->mDelay ,ret->mpBlock[block_id-1]->mDelay);
    }
    ret->mInSize = ret->mOutSize = 0;
    for (j=0; j < ret->mNBlocks; j++) {
      ret->mInSize  = HIGHER_OF(ret->mInSize,  ret->mpBlock[j]->mInSize + nblocks);
      ret->mOutSize = HIGHER_OF(ret->mOutSize, ret->mpBlock[j]->mOutSize);
    }
    ret->mBlockOutputVector.Init(ret->mOutSize);
    ret->mpMacro = macro;
    return ret;
  }; // ReadCompositeXform(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************    
  //***************************************************************************
  LinearXform *
  ModelSet::
  ReadLinearXform(FILE *fp, Macro *macro, bool predefined)
  {
    LinearXform *   ret;
    int             in_size;
    int             out_size;
    size_t          i;
    int             r;
    int             c;
    FLOAT           tmp;
  
    // read the parameters
    out_size = GetInt(fp);  // Rows in Xform matrix
    in_size =  GetInt(fp);   // Cols in Xform matrix
  
    // create new object
    ret = new LinearXform(in_size, out_size);
    ret->mPredefinedID = PLTID_NONE;

    // generating of predefined transforms 
    // - predefined transforms are suported in text format only
    if(predefined)
    {
      char *p_tr_name = GetString(fp, false);
      for(i = 0; i < strlen(p_tr_name); i++) 
      {
        p_tr_name[i] = toupper(p_tr_name[i]);
      }
    
      if(!strcmp(p_tr_name, "DCT")) 
      {
        ret->mPredefinedID = PLTID_DCT;            
        ret->mIncludeC0 = (GetInt(fp) != 0);
  
        if(out_size > in_size || (out_size == in_size && !ret->mIncludeC0))
        {
          Error("Number of base components (%d) of DCT can not be higher then input vector size (%d) (%s:%d)", out_size, in_size, mpCurrentMmfName, mCurrentMmfLine);
        }
  
        GenerateDCTMatrix(ret->mMatrix, out_size, in_size, 1, ret->mIncludeC0);	    
      } 
      else if(!strcmp(p_tr_name, "CONST"))
      {
         ret->mPredefinedID = PLTID_CONST;
         ret->mConstant = GetFloat(fp);

         GenerateConstantMatrix(ret->mMatrix, out_size, in_size, ret->mConstant);
      }
      else if(!strcmp(p_tr_name, "DIAG"))
      {
         ret->mPredefinedID = PLTID_DIAG;
         ret->mConstant = GetFloat(fp);
         
         GenerateDiagMatrix(ret->mMatrix, out_size, in_size, ret->mConstant);         
      }
      else if(!strcmp(p_tr_name, "RANDOM"))
      {
        ret->mPredefinedID = PLTID_RANDOM;
        ret->mMinValue = GetFloat(fp);
        ret->mMaxValue = GetFloat(fp);
        ret->mSeed = static_cast<unsigned int>(GetInt(fp));
            
        GenerateRandomMatrix(ret->mMatrix, out_size, in_size, ret->mMinValue, ret->mMaxValue, ret->mSeed);
      }
      else 
      {
        Error("Unknown predefined Xform %s (%s:%d)", p_tr_name, mpCurrentMmfName, mCurrentMmfLine);	  
      }
    } 
    else
    {
      i = 0;
      for (c=0; c < out_size; c++)
      {
        for (r = 0; r < in_size; r++)
        {
          tmp = GetFloat(fp);
          ret->mMatrix[c][r] = tmp;                
          i++;
        }
      }
    }
  
    ret->mpMacro      = macro;
    return ret;
  }; // ReadLinearXform(FILE *fp, Macro *macro)
  

  //***************************************************************************  
  //***************************************************************************  
  BiasXform *
  ModelSet::
  ReadBiasXform(FILE *fp, Macro *macro, bool predefined)
  {
    BiasXform *   ret;
    size_t        size;
    size_t        i;
    
    size = GetInt(fp);
    ret  = new BiasXform(size);
    
    if(predefined)
    {
      ret->mUsePredefVector = true;
      ret->mPredefVector.Read(fp, size, this);
      
      // This is a hack, BiasXform should have vector instead of matrix.
      // !!! It is hard to fix this in SNet
      BasicVector<FLOAT> vct;
      ret->mPredefVector.Generate(vct);
      
      for(i = 0; i < size; i++)
      {
        ret->mVector[0][i] = vct[i];
      }
    } 
    else
    {
      // fill the object with data    
      for (i=0; i < size; i++) 
      {
        ret->mVector[0][i] = GetFloat(fp);
      }
    }
      
    ret->mpMacro      = macro;
    return ret;
  }; // ReadBiasXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  FuncXform *
  ModelSet::
  ReadFuncXform(FILE *fp, Macro *macro, int funcId)
  {
    FuncXform * ret;
    size_t      size;
  
    //***
    // old malloc
    //ret = (FuncXform *) malloc(sizeof(FuncXform));
    //if (ret == NULL) Error("Insufficient memory");
    
    
    size          = GetInt(fp);
    ret           = new FuncXform(size, funcId);
    ret->mpMacro  = macro;
    return ret;
  }; // ReadFuncXform(FILE *fp, Macro *macro, int funcId)
  
  
  //***************************************************************************
  //***************************************************************************
  CopyXform *
  ModelSet::
  ReadCopyXform(FILE *fp, Macro *macro)
  {
    CopyXform *   ret;
    int           in_size;
    int           out_size;
    int           i=0;
    int           n;
    int           from;
    int           step;
    int           to;
  
    out_size = GetInt(fp);
    in_size = GetInt(fp);
  
    //***
    // old malloc
    //ret = (CopyXform *) malloc(sizeof(CopyXform)+(out_size-1)*sizeof(int));
    //if (ret == NULL) Error("Insufficient memory");
    
    // create new object
    ret = new CopyXform(in_size, out_size);
  
    // fill it with data
    while (i < out_size) 
    {
      RemoveSpaces(fp);
      if ((n = fscanf(fp, "%d:%d:%d", &from, &step, &to)) < 1) {
        if (ferror(fp)) {
          Error("Cannot read input file %s", mpCurrentMmfName);
        }
        Error("Integral number expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
      }
  
      if (n == 2)      { to = step; step = 1; }
      else if (n == 1) { to = from; step = 1; }
  
      if (to < 1 || to > in_size) {
        Error("Copy index %d out of range (%s:%d)",
              to, mpCurrentMmfName, mCurrentMmfLine);
      }
  
      for (n = 0; n < (to-from)/step + 1; n++, i++) {
        ret->mpIndices[i] = from + n * step - 1;
      }
    }
  
    ret->mpMacro      = macro;
    return ret;
  }; //ReadCopyXform(FILE *fp, Macro *macro)

  //***************************************************************************
  //***************************************************************************
  BlockCopyXform *
  ModelSet::
  ReadBlockCopyXform(FILE *fp, Macro *macro)
  {
    BlockCopyXform *   ret;
    int           n;
    int           from;
    int           step;
    int           to;

    int           out_size = 0;
      
    int n_blocks = GetInt(fp);
    int n_rows = GetInt(fp);
    int n_cols = GetInt(fp);
  
    BlockCopyXform::Block *p_blocks = new BlockCopyXform::Block[n_blocks];

    int i;        
    int j;
    int k;
    for(i = 0; i < n_blocks; i++)
    {
      if(!mHmmReadBinary)
      {
        // first dimension
        RemoveSpaces(fp);
        if ((n = fscanf(fp, "%d:%d:%d", &from, &step, &to)) < 1) {
          if (ferror(fp)) {
            Error("Cannot read input file %s", mpCurrentMmfName);
          }
          Error("Integral number expected (%s:%d)", mpCurrentMmfName, 
              mCurrentMmfLine);
        }
      
        if (n == 2)      { to = step; step = 1; }
        else if (n == 1) { to = from; step = 1; }
  
        p_blocks[i].mFromR = from;
        p_blocks[i].mStepR = step;
        p_blocks[i].mToR = to;
      
        // second dimension
        char buff[256];
        if(fscanf(fp, "%255[,]", buff) == 1)
        {
          if ((n = fscanf(fp, "%d:%d:%d", &from, &step, &to)) < 1) {
            if (ferror(fp)) {
              Error("Cannot read input file %s", mpCurrentMmfName);
            }
            Error("Integral number expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
          }
      
          if (n == 2)      { to = step; step = 1; }
          else if (n == 1) { to = from; step = 1; }
  
          p_blocks[i].mFromC = from;
          p_blocks[i].mStepC = step;
          p_blocks[i].mToC = to;                  
        }
        else
        {
          p_blocks[i].mFromC = p_blocks[i].mFromR;
          p_blocks[i].mStepC = p_blocks[i].mStepR;
          p_blocks[i].mToC = p_blocks[i].mToR;      

          p_blocks[i].mFromR = 1;
          p_blocks[i].mStepR = 1;
          p_blocks[i].mToR = 1;      
        }
      }
      else
      {
          p_blocks[i].mFromR = GetInt(fp);
          p_blocks[i].mStepR = GetInt(fp);
          p_blocks[i].mToR = GetInt(fp);
          
          p_blocks[i].mFromC = GetInt(fp);
          p_blocks[i].mStepC = GetInt(fp);
          p_blocks[i].mToC = GetInt(fp);
      }
                  
      if (p_blocks[i].mToR < 1 || p_blocks[i].mToR > n_rows) {
        Error("BlockCopy index %d out of range (%s:%d)",
              p_blocks[i].mToR, mpCurrentMmfName, mCurrentMmfLine);
      }
      
      if (p_blocks[i].mToC < 1 || p_blocks[i].mToC > n_cols) {
        Error("BlockCopy index %d out of range (%s:%d)",
              p_blocks[i].mToC, mpCurrentMmfName, mCurrentMmfLine);
      }
        
      int n_vals_r = 0;
      for(j = p_blocks[i].mFromR; j <= p_blocks[i].mToR; j += p_blocks[i].mStepR)
      {
        n_vals_r++;
      }
      int n_vals_c = 0;
      for(j = p_blocks[i].mFromC; j <= p_blocks[i].mToC; j += p_blocks[i].mStepC)
      {
        n_vals_c++;
      }
      
      out_size += (n_vals_r * n_vals_c);        
      
//      printf(" %d:%d:%d,%d:%d:%d\n", p_blocks[i].mFromR, p_blocks[i].mStepR, p_blocks[i].mToR,
//                                     p_blocks[i].mFromC, p_blocks[i].mStepC, p_blocks[i].mToC);
                         
    }
        
    // create new object
    ret = new BlockCopyXform(n_blocks);
    ret->mpBlocks  = p_blocks;
    ret->mNBlocks  = n_blocks;
    ret->mNRows    = n_rows;
    ret->mInSize  = n_rows * n_cols;
    ret->mOutSize = out_size;
    ret->mpIndices = new int [out_size];
    
    // fill in indices
    int idx = 0;
    for(i = 0; i < n_blocks; i++)
    {
      for(j = p_blocks[i].mFromR; j <= p_blocks[i].mToR; j += p_blocks[i].mStepR)
      {
        for(k = p_blocks[i].mFromC; k <= p_blocks[i].mToC; k += p_blocks[i].mStepC)
        {
          ret->mpIndices[idx] = (j - 1) * n_cols + (k - 1);
          idx++;
        }
      }
    }  
    
    ret->mpMacro      = macro;
    return ret;
  }; //ReadCopyXform(FILE *fp, Macro *macro)

  
  //***************************************************************************
  //***************************************************************************
  TransposeXform *
  ModelSet::
  ReadTransposeXform(FILE *fp, Macro *macro)
  {
    TransposeXform *  ret;
  
    size_t in_rows = GetInt(fp);
    size_t in_cols = GetInt(fp);
      
    // create new object
    ret = new TransposeXform(in_rows, in_cols);
    
    ret->mpMacro      = macro;
    return ret;
  }; //ReadTransposeXform(FILE *fp, Macro *macro)

  //***************************************************************************
  //***************************************************************************
  ConstantXform *
  ModelSet::
  ReadConstantXform(FILE *fp, Macro *macro)
  {
    int        out_size = GetInt(fp);
    ConstantXform * ret = new ConstantXform(out_size);
  
    // load values
    for (int i=0; i < out_size; i++) {
      ret->mVector[0][i]  = GetFloat(fp);
    }
    ret->mpMacro      = macro;
    return ret;
  }; //ReadConstantXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************  
  //***************************************************************************  
  WindowXform *
  ModelSet::
  ReadWindowXform(FILE *fp, Macro *macro, bool predefined)
  {
    WindowXform *   ret;
    size_t        size;
    size_t        i;
    
    size = GetInt(fp);
    ret  = new WindowXform(size);
    
    if(predefined)
    {
      ret->mUsePredefVector = true;
      ret->mPredefVector.Read(fp, size, this);
      ret->mPredefVector.Generate(ret->mVector);
    } 
    else
    {
      // fill the object with data    
      for (i=0; i < size; i++) 
      {
        ret->mVector[i] = GetFloat(fp);
      }
    }
  
    ret->mpMacro      = macro;
    return ret;
  }; // ReadWindowXform(FILE *fp, Macro *macro)

  
  //***************************************************************************
  //***************************************************************************
  FeatureMappingXform*
  ModelSet::
  ReadFeatureMappingXform(FILE* fp, Macro* macro)
  {
    FeatureMappingXform*    ret=NULL;
    int           in_size;
    int           in_size_orig;
    char*         keyword;
    int           state_id;
    
    in_size = GetInt(fp);
  
    // create new object
    ret = new FeatureMappingXform(in_size);
    
    // we need to tell the model set to ignore the original input vector size
    // as the transform can have different dimensiion. The value needs to be 
    // later restored
    in_size_orig     = mInputVectorSize;
    mInputVectorSize = in_size;
    
    // read working states
    keyword = GetString(fp, 1);
    
    if (!CheckKwd(keyword, KID_State)) 
      Error("Keyword <State> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    state_id = GetInt(fp);
 
    if (state_id !=1) 
      Error("State number should be 1 (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    ret->mpStateFrom = ReadState(fp, NULL);
    
    
    keyword = GetString(fp, 1);
    
    if (!CheckKwd(keyword, KID_State)) 
      Error("Keyword <State> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
  
    state_id = GetInt(fp);
 
    if (state_id != 2) 
      Error("State number should be 2 (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    
    ret->mpStateTo   = ReadState(fp, NULL);
  
    
    // check for alike number of mixtures in both states
    if (ret->mpStateFrom->mNMixtures != ret->mpStateTo->mNMixtures)
    {
      Error("Number of mixtures must match in the Xform states (%s:%d)", 
        mpCurrentMmfName, mCurrentMmfLine);
    }
    
    // restore the original input vector size value
    mInputVectorSize  = in_size_orig;
    
    ret->mpMacro      = macro;
    return ret;
  }; //ReadFeatureMappingXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  GmmPosteriorsXform*
  ModelSet::
  ReadGmmPosteriorsXform(FILE* fp, Macro* macro)
  {
    GmmPosteriorsXform* ret=NULL;
    int           in_size_orig;
    //char*         keyword;
    //int           state_id;
    
    int   in_size = GetInt(fp);
    size_t n_best = GetInt(fp);
    FLOAT scale   = GetFloat(fp);

    // create new object
    ret = new GmmPosteriorsXform(in_size);

    // we need to tell the model set to ignore the original input vector size
    // as the transform can have different dimensiion. The value needs to be 
    // later restored
    in_size_orig     = mInputVectorSize;
    mInputVectorSize = in_size;
    
    ret->mpState = ReadState(fp, NULL);    
    ret->mOutSize = ret->mpState->mNMixtures;
    n_best = LOWER_OF(n_best, ret->mOutSize);
    ret->mNBestLikesVector.resize(n_best);
    ret->mScale = scale;    
    
    // restore the original input vector size value
    mInputVectorSize  = in_size_orig;
    
    ret->mpMacro      = macro;
    return ret;
  }; //ReadGmmPosteriorsXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  FrantaProductXform*
  ModelSet::
  ReadFrantaProductXform(FILE* fp, Macro* macro)
  {
    FrantaProductXform*   ret;
    int                   in_size;
    int                   n_parts;
      
    in_size = GetInt(fp);
    n_parts = GetInt(fp);
    
    if (in_size % n_parts)
    {
      Error("Vector size has to be divisible by %d (%s:%d)", 
        n_parts, mpCurrentMmfName, mCurrentMmfLine);
    }  
    
    // create new object
    ret = new FrantaProductXform(in_size, n_parts);
    
    ret->mpMacro      = macro;
    return ret;
  }; //ReadFrantaProductXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  MatlabXform*
  ModelSet::
  ReadMatlabXform(FILE* fp, Macro* macro)
  {
#ifndef MATLAB_ENGINE
    Warning("Matlab engine not available, transform will have unexpected results");    
#endif    
    
    MatlabXform*          ret;
    int                   in_rows;
    int                   in_cols;
    int                   out_size;
    ReadlineData          rld = {0};
    char*                 line;
          
    in_rows   = GetInt(fp);
    in_cols   = GetInt(fp);
    out_size  = GetInt(fp);
    
    // create new object
    ret = new MatlabXform(in_rows, in_cols, out_size);
    
    while ((line = readline(fp, &rld)) != NULL)
    {
      if (!strcmp(line, "%%%% END %%%%"))
        break;
        
      ret->mProgram += line;
      ret->mProgram += '\n';
    }
  
    free(rld.buffer);
    
    ret->mpMacro      = macro;
    return ret;
  }; //MatlabXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  StackingXform *
  ModelSet::
  ReadStackingXform(FILE *fp, Macro *macro)
  {
    StackingXform *   ret;
    size_t            stack_size;
    size_t            in_size;
  
    stack_size = GetInt(fp);
    in_size    = GetInt(fp);
  
    ret = new StackingXform(stack_size, in_size);
  
    ret->mpMacro      = macro;
    return ret;
  }; //ReadStackingXform(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  ModelSet::
  ReadGlobalOptions(FILE *fp, bool readOnly)
  {
    int           i; 
    int           ret = 0;
    char *        keyword;
  
  //puts("ReadGlobalOptions");
    for (;;) 
    {
      if ((keyword = GetString(fp, 0)) == NULL) 
      {
        return ret;
  //      Error("Unexpected end of file %s", mpCurrentMmfName);
      }
  
      if (CheckKwd(keyword, KID_VecSize)) 
      {
        i = GetInt(fp);
        if (mInputVectorSize != -1 && mInputVectorSize != i) 
        {
          Error("Mismatch in <VecSize> redefinition (%s:%d)",
                mpCurrentMmfName, mCurrentMmfLine);
        }
  
        mInputVectorSize = i;
        
        ret = 1;
      } 
      else if (CheckKwd(keyword, KID_StreamInfo)) 
      {
        if (GetInt(fp) != 1) 
        {
          Error("Unsupported definition of multistream (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        }
  
        i = GetInt(fp);
        if (mInputVectorSize != -1 && mInputVectorSize != i) 
        {
          Error("Mismatch in <VecSize> redefinition (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        }
  
        mInputVectorSize = i;
        ret = 1;
      } 
      else if ((i = ReadParmKind(keyword, true)) != -1) 
      {
        if (mParamKind != -1 && mParamKind != i) 
          Error("Mismatch in paramKind redefinition (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        
        mParamKind = i;
        ret = 1;
      } 
      else if ((i = ReadOutPDFKind(keyword)) != -1) 
      {
        if (mOutPdfKind != -1 && mOutPdfKind != i) {
          Error("Mismatch in outPDFKind redefinition (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        }
  
        if (i != KID_PDFObsVec && i != KID_DiagC) {
          Error("Unsupported option '%s' (%s:%d)", keyword, mpCurrentMmfName, mCurrentMmfLine);
        }
  
        mOutPdfKind = static_cast<KeywordID> (i);
        ret = 1;
      } 
      else if ((i = ReadDurKind(keyword)) != -1) 
      {
        if (mDurKind != -1 && mDurKind != i) {
          Error("Mismatch in durKind redefinition (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
        }
        if (i != KID_NullD) {
          Error("Unsupported option '%s' (%s:%d)", keyword, mpCurrentMmfName, mCurrentMmfLine);
        }
        mDurKind = static_cast<KeywordID> (i);
        ret = 1;
      } else if (CheckKwd(keyword, KID_HMMSetID)) {
        ret = 1;
      } else if (CheckKwd(keyword, KID_InputXform)) {
        XformInstance *inputXform;
        Macro *macro = pAddMacro(mt_XformInstance, DEFAULT_XFORM_NAME);

        macro->mReadOnly = readOnly;
        
        if (macro->mpData != NULL) 
        {
          if (mHmmsIgnoreMacroRedefinition == 0)
            Error("Redefinition of <InputXform> (%s:%d)",
                  mpCurrentMmfName, mCurrentMmfLine);
        }
        
        inputXform = ReadXformInstance(fp, macro);
  
        if (macro->mpData != NULL) {
          Warning("Redefinition of <InputXform> (%s:%d)",
                  mpCurrentMmfName, mCurrentMmfLine);
  
          // Macro is redefined. New item must be checked for compatibility with
          // the old one (vector size) All references to old
          // item must be replaced and old item must be released
  
          ReplaceItemUserData ud;
          ud.mpOldData = macro->mpData;
          ud.mpNewData = inputXform;
          ud.mType     = 'j';
  
          this         ->Scan(MTM_XFORM_INSTANCE|MTM_MIXTURE,NULL,ReplaceItem, &ud);
          ud.mpOldData ->Scan(MTM_REVERSE_PASS|MTM_ALL, NULL, ReleaseItem, NULL);
  
          for (unsigned int i = 0; i < mXformInstanceHash.mNEntries; i++) 
          {
            if (mXformInstanceHash.mpEntry[i]->data == ud.mpOldData) 
            {
              mXformInstanceHash.mpEntry[i]->data = ud.mpNewData;
            }
          }
        } 
        else 
        {
          mpInputXform = inputXform;
          macro->mpData = inputXform;
        }
        
        ret = 1;
  //    } else if (CheckKwd(keyword, KID_LinXform)) {
  //      linXform = ReadXformInstance(fp, hmm_set, NULL);
        ret = 1;
      } else {
        UngetString();
        return ret;
      }
    }
  }; //ReadGlobalOptions(FILE *fp)
  
  
  
  //***************************************************************************
  //***************************************************************************
  Transition *
  ModelSet::
  ReadTransition(FILE *fp, Macro *macro)
  {
    Transition *ret;
    char *keyword;
    int nstates, i;
  
  //  puts("ReadTransition");
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~t")) {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mTransitionHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~t %s (%s:%d)", keyword, mpCurrentMmfName, mCurrentMmfLine);
      }
      return (Transition *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_TransP)) {
      Error("Keyword <TransP> expected (%s:%d)", mpCurrentMmfName, mCurrentMmfLine);
    }
  
    nstates = GetInt(fp);
  
    //***
    // old malloc
    //i = mAllocAccums ? 2 * SQR(nstates) + nstates: SQR(nstates);
    //
    //if ((ret = (Transition *) malloc(sizeof(Transition) + i * sizeof(ret->mpMatrixO[0]))) == NULL) {
    //  Error("Insufficient memory");
    //}
    //
    //ret->mNStates = nstates;
    //
    
    ret = new Transition(nstates, mAllocAccums);
    
    for (i=0; i < SQR(nstates); i++) {
      ret->mpMatrixO[i] = GetFloat(fp);
      ret->mpMatrixO[i] = (float) ret->mpMatrixO[i] != 0.0 ? my_log(ret->mpMatrixO[i]) : LOG_0;
    }
  
    ret->mpMacro = macro;
    return ret;
  }
  

    
  //############################################################################
  
  //############################################################################
  //############################################################################
  // MMF OUTPUT
  //############################################################################
  //############################################################################
    
  //****************************************************************************
  //****************************************************************************
  void 
  ModelSet::
  WriteMmf(const char* pFileName, const char* pOutputDir,
           const char* pOutputExt, bool binary)
  {
    FILE*       fp = NULL;
    Macro*      macro;
    char        mmfile[1024];
    char*       lastFileName = NULL;
    int         waitingForNonXform = 1;
    OStkStream  output_stream;
  
    for (macro = mpFirstMacro; macro != NULL; macro = macro->nextAll) 
    {
      // we don't store the read-only macros
      // TODO: check whether this is a good idea.
      if (macro->mReadOnly) {
        continue;
      }

      if (macro->mpFileName == NULL) {
        continue; // Artificial macro not read from file
      }
  
      if (lastFileName == NULL 
      || (!pFileName && strcmp(lastFileName, macro->mpFileName))) 
      {
        // New macro file
        lastFileName = macro->mpFileName;
        
        if (output_stream.is_open()) {
          output_stream.close();
        }
  
        MakeFileName(mmfile, pFileName ? pFileName : macro->mpFileName, pOutputDir, pOutputExt);
          
        output_stream.open(mmfile, std::ios::binary, gpMmfOFilter);
        if (!output_stream.good()) { 
          Error("Cannot open output MMF %s", mmfile);
        }

        fp = output_stream.file();
          
        waitingForNonXform = 1;
        WriteGlobalOptions(fp, binary);
      }
  
      if (macro->mpData == mpInputXform &&
        (!strcmp(macro->mpName, DEFAULT_XFORM_NAME))) 
      {
        fputs("~o ", fp);
        PutKwd(fp, binary, KID_InputXform);
        PutNLn(fp, binary);
      } 
      else 
      {
        fprintf(fp, "~%c \"%s\"", macro->mType, macro->mpName);
        PutNLn(fp, binary);
      }
  
      if (macro->mpData->mpMacro != macro) 
      {
        fprintf(fp, " ~%c \"%s\"", macro->mType, (macro->mpData->mpMacro)->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        switch (macro->mType) 
        {
          case 'x': WriteXform        (fp, binary, static_cast <Xform*>         (macro->mpData)); break;
          case 'j': WriteXformInstance(fp, binary, static_cast <XformInstance*> (macro->mpData)); break;
          case 'u': WriteMean         (fp, binary, static_cast <Mean*>          (macro->mpData)); break;
          case 'v': WriteVariance     (fp, binary, static_cast <Variance*>      (macro->mpData)); break;
          case 't': WriteTransition   (fp, binary, static_cast <Transition*>    (macro->mpData)); break;
          case 'm': WriteMixture      (fp, binary, static_cast <Mixture*>       (macro->mpData)); break;
          case 's': WriteState        (fp, binary, static_cast <State*>         (macro->mpData)); break;
          case 'h': WriteHMM          (fp, binary, static_cast <Hmm*>           (macro->mpData)); break;
        }
      }
    }

    if (output_stream.is_open()) {
      output_stream.close();
    }
  }
  
  //****************************************************************************
  //****************************************************************************
  void 
  ModelSet::
  WriteMmfRaw(const char* pFileName, const char* pOutputDir,
                 const char* pOutputExt, bool binary)
  {
    FILE*       fp = NULL;
    Macro*      macro;
    char        mmfile[1024];
    char*       lastFileName = NULL;
    int         waitingForNonXform = 1;
    OStkStream  output_stream;
  
    for (macro = mpFirstMacro; macro != NULL; macro = macro->nextAll) 
    {
      // we don't store the read-only macros
      // TODO: check whether this is a good idea.
      if (macro->mReadOnly) {
        continue;
      }

      if (macro->mpFileName == NULL) {
        continue; // Artificial macro not read from file
      }
  
      if (lastFileName == NULL 
      || (!pFileName && strcmp(lastFileName, macro->mpFileName))) 
      {
        // New macro file
        lastFileName = macro->mpFileName;
        
        if (output_stream.is_open()) {
          output_stream.close();
        }
  
        MakeFileName(mmfile, pFileName ? pFileName : macro->mpFileName, pOutputDir, pOutputExt);
          
        output_stream.open(mmfile, std::ios::binary, gpMmfOFilter);
        if (!output_stream.good()) { 
          Error("Cannot open output MMF %s", mmfile);
        }

        fp = output_stream.file();
          
        waitingForNonXform = 1;
        WriteGlobalOptions(fp, binary);
      }
  
      if (macro->mpData == mpInputXform &&
        (!strcmp(macro->mpName, DEFAULT_XFORM_NAME))) 
      {
        fputs("~o ", fp);
        PutKwd(fp, binary, KID_InputXform);
        PutNLn(fp, binary);
      } 
      else 
      {
        fprintf(fp, "~%c \"%s\"", macro->mType, macro->mpName);
        PutNLn(fp, binary);
      }
  
      if (macro->mpData->mpMacro != macro) 
      {
        fprintf(fp, " ~%c \"%s\"", macro->mType, (macro->mpData->mpMacro)->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        switch (macro->mType) 
        {
          case 'u': WriteMeanRaw    (fp, binary, static_cast <Mean*>     (macro->mpData)); break;
          case 'v': WriteVarianceRaw(fp, binary, static_cast <Variance*> (macro->mpData)); break;
          case 'm': WriteMixtureRaw (fp, binary, static_cast <Mixture*>  (macro->mpData)); break;
          case 's': WriteStateRaw   (fp, binary, static_cast <State*>    (macro->mpData)); break;
          case 'h': WriteHMMRaw     (fp, binary, static_cast <Hmm*>      (macro->mpData)); break;
          default: break;
        }
      }
    }

    if (output_stream.is_open()) {
      output_stream.close();
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  WriteGlobalOptions(FILE *fp, bool binary)
  {
    char parmkindstr[64];

    if(mSaveGlobOpts == false)
    {
      return;
    }

    fputs("~o ", fp);
    PutKwd(fp, binary, KID_VecSize);
    PutInt(fp, binary, mInputVectorSize);
  
    if (ParmKind2Str(mParamKind, parmkindstr))
    {
      fprintf(fp, "<%s> ", parmkindstr);
    }
    
    if (mOutPdfKind != -1)
    {
      PutKwd(fp, binary, mOutPdfKind);
    }
    
    if (mOutPdfKind != -1)
    {
      PutKwd(fp, binary, mDurKind);
    }
    
    PutNLn(fp, binary);
  }
  
    
  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  WriteHMM(FILE *fp, bool binary, Hmm *hmm)
  {
    size_t i;
  
    PutKwd(fp, binary, KID_BeginHMM);
    PutNLn(fp, binary);
    PutKwd(fp, binary, KID_NumStates);
    PutInt(fp, binary, hmm->mNStates);
    PutNLn(fp, binary);
  
    for (i=0; i < hmm->mNStates-2; i++) 
    {
      PutKwd(fp, binary, KID_State);
      PutInt(fp, binary, i+2);
      PutNLn(fp, binary);
  
      if (hmm->mpState[i]->mpMacro) 
      {
        fprintf(fp, "~s \"%s\"", hmm->mpState[i]->mpMacro->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        WriteState(fp, binary, hmm->mpState[i]);
      }
    }
    
    if (hmm->mpTransition->mpMacro) 
    {
      fprintf(fp, "~t \"%s\"", hmm->mpTransition->mpMacro->mpName);
      PutNLn(fp, binary);
    } 
    else 
    {
      WriteTransition(fp, binary, hmm->mpTransition);
    }
    PutKwd(fp, binary, KID_EndHMM);
    PutNLn(fp, binary);
  } // WriteHMM(FILE *fp, bool binary, Hmm *hmm)
  

  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  WriteHMMRaw(FILE *fp, bool binary, Hmm *hmm)
  {
    size_t i;
  
    for (i=0; i < hmm->mNStates-2; i++) {
      if (!hmm->mpState[i]->mpMacro) {
        WriteState(fp, binary, hmm->mpState[i]);
      }
    }
  } // WriteHMM(FILE *fp, bool binary, Hmm *hmm)
  

  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteState(FILE* fp, bool binary, State* state)
  {
    size_t i;
  
    if (mOutPdfKind == KID_PDFObsVec) 
    {
      PutKwd(fp, binary, KID_ObsCoef);
      PutInt(fp, binary, state->PDF_obs_coef);
      PutNLn(fp, binary);
    } 
    else 
    {
      if (state->mNMixtures > 1) 
      {
        PutKwd(fp, binary, KID_NumMixes);
        PutInt(fp, binary, state->mNMixtures);
        PutNLn(fp, binary);
      }
  
      for (i=0; i < state->mNMixtures; i++) 
      {
        if (state->mNMixtures > 1) 
        {
          PutKwd(fp, binary, KID_Mixture);
          PutInt(fp, binary, i+1);
          PutFlt(fp, binary, my_exp(state->mpMixture[i].mWeight));
          PutNLn(fp, binary);
        }
  
        if (state->mpMixture[i].mpEstimates->mpMacro) 
        {
          fprintf(fp, "~m \"%s\"", state->mpMixture[i].mpEstimates->mpMacro->mpName);
          PutNLn(fp, binary);
        } 
        else 
        {
          WriteMixture(fp, binary, state->mpMixture[i].mpEstimates);
        }
      }
    }
  }
  
  
  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteStateRaw(FILE* fp, bool binary, State* state)
  {
    size_t i;
  
    if (mOutPdfKind != KID_PDFObsVec) 
    {
      for (i=0; i < state->mNMixtures; i++) 
      {
        PutFlt(fp, binary, my_exp(state->mpMixture[i].mWeight));
        PutNLn(fp, binary);
  
        if (!state->mpMixture[i].mpEstimates->mpMacro) {
          WriteMixtureRaw(fp, binary, state->mpMixture[i].mpEstimates);
        }
      }
    }
  }
  
  
  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteMixture(FILE *fp, bool binary, Mixture *mixture)
  {
    //:BUG:
    // The mixture->mpInputXform should never be null, however with feature
    // mapping, it is. Temporary solution is to test it for NULL...
    if (mixture->mpInputXform != mpInputXform
    &&  NULL != mixture->mpInputXform) 
    {
      PutKwd(fp, binary, KID_InputXform);
      
      if (mixture->mpInputXform->mpMacro) 
      {
        fprintf(fp, "~j \"%s\"", mixture->mpInputXform->mpMacro->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        WriteXformInstance(fp, binary, mixture->mpInputXform);
      }
    }
    
    if (mixture->mpMean->mpMacro) 
    {
      fprintf(fp, "~u \"%s\"", mixture->mpMean->mpMacro->mpName);
      PutNLn(fp, binary);
    } 
    else 
    {
      WriteMean(fp, binary, mixture->mpMean);
    }
    
    if (mixture->mpVariance->mpMacro) 
    {
      fprintf(fp, "~v \"%s\"", mixture->mpVariance->mpMacro->mpName);
      PutNLn(fp, binary);
    } 
    else 
    {
      WriteVariance(fp, binary, mixture->mpVariance);
    }
    
    PutKwd(fp, binary, KID_GConst);
    PutFlt(fp, binary, mixture->GConst());
    PutNLn(fp, binary);
  }


  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteMixtureRaw(FILE *fp, bool binary, Mixture *mixture)
  {
    if (!mixture->mpMean->mpMacro) {
      WriteMeanRaw(fp, binary, mixture->mpMean);
    }
    
    if (!mixture->mpVariance->mpMacro) {
      WriteVarianceRaw(fp, binary, mixture->mpVariance);
    }
  }


  //***************************************************************************  
  //*****************************************************************************  
  void
  ModelSet::
  WriteMean(FILE *fp, bool binary, Mean *mean)
  {
    size_t    i;
  
    // Mean specified by a set of vectors
    if (NULL != mean->mpClusterWeightVectors)
    //if(false)
    {
      PutKwd(fp, binary, KID_Weights);
      PutInt(fp, binary, mean->mNClusterWeightVectors);
      PutNLn(fp, binary);
      for (i = 0; i < mean->mNClusterWeightVectors; i++)
      {
        fprintf(fp, "~x \"%s\"", mean->mpClusterWeightVectors[i]->mpMacro->mpName);
        PutNLn(fp, binary);
      }
      
      // go through the rows in the vector
      for (i = 0; i < mean->mClusterMatrixT.Rows(); i++)
      {
        PutKwd(fp, binary, KID_Mean);
        PutInt(fp, binary, mean->VectorSize());
        PutNLn(fp, binary);
      
        for (size_t j=0; j < mean->VectorSize(); j++) 
        {
          PutFlt(fp, binary, mean->mClusterMatrixT(i, j));
        }
      
        PutNLn(fp, binary);
      }
    }
    
    // Mean specified regularly by a single vector
    else
    {
      PutKwd(fp, binary, KID_Mean);
      PutInt(fp, binary, mean->VectorSize());
      PutNLn(fp, binary);
    
      for (i=0; i < mean->VectorSize(); i++) 
      {
        PutFlt(fp, binary, mean->mVector[i]);
      }
    
      PutNLn(fp, binary);
    }
  } //WriteMean(FILE *fp, bool binary, Mean *mean)

  
  //***************************************************************************  
  //*****************************************************************************  
  void
  ModelSet::
  WriteMeanRaw(FILE *fp, bool binary, Mean *mean)
  {
    size_t    i;
  
    for (i=0; i < mean->VectorSize(); i++) 
    {
      PutFlt(fp, binary, mean->mVector[i]);
    }
  
    PutNLn(fp, binary);
  } //WriteMeanRaw(FILE *fp, bool binary, Mean *mean)

  
  //***************************************************************************  
  //*****************************************************************************  
  void
  ModelSet::
  WriteVariance(FILE *fp, bool binary, Variance *variance)
  {
    size_t   i;
  
    PutKwd(fp, binary, KID_Variance);
    PutInt(fp, binary, variance->VectorSize());
    PutNLn(fp, binary);
  
    for (i=0; i < variance->VectorSize(); i++) 
    {
      PutFlt(fp, binary, 1/variance->mVector[i]);
    }
  
    PutNLn(fp, binary);
  }

  
  //***************************************************************************  
  //*****************************************************************************  
  void
  ModelSet::
  WriteVarianceRaw(FILE *fp, bool binary, Variance *variance)
  {
    size_t   i;
  
    PutKwd(fp, binary, KID_Variance);
    PutInt(fp, binary, variance->VectorSize());
    PutNLn(fp, binary);
  
    for (i=0; i < variance->VectorSize(); i++) 
    {
      PutFlt(fp, binary, 1/variance->mVector[i]);
    }
  
    PutNLn(fp, binary);
  }

  
  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteTransition(FILE *fp, bool binary, Transition *transition)
  {
    size_t  i;
    size_t  j;
  
    PutKwd(fp, binary, KID_TransP);
    PutInt(fp, binary, transition->mNStates);
    PutNLn(fp, binary);
  
    for (i=0; i < transition->mNStates; i++) 
    {
      for (j=0; j < transition->mNStates; j++) 
      {
        FLOAT logtp = transition->mpMatrixO[i * transition->mNStates + j];
        PutFlt(fp, binary, logtp > LOG_MIN ? my_exp(logtp) : 0.0);
      }
  
      PutNLn(fp, binary);
    }
  }

  
  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteXformInstance(FILE *fp, bool binary, XformInstance *xformInstance)
  {
    size_t            i;
    bool              isHTKCompatible = true;
    char              parmkindstr[64];
  
    CompositeXform *cxf = (CompositeXform *) xformInstance->mpXform;
    
    if (xformInstance->mpInput != NULL || cxf == NULL || cxf->mpMacro ||
      cxf->mXformType != XT_COMPOSITE || cxf->mNLayers != 1) 
    {
      isHTKCompatible = 0;
    } 
    else
    {
      for (i = 0; i < cxf->mpLayer[0].mNBlocks; i++) 
      {
        if (cxf->mpLayer[0].mpBlock[i]->mXformType != XT_LINEAR) 
        {
          isHTKCompatible = 0;
          break;
        }
      }
    }
  
    if (xformInstance->mpInput != NULL) 
    {
      PutKwd(fp, binary, KID_Input);
      PutNLn(fp, binary);
  
      if (xformInstance->mpInput->mpMacro) 
      {
        fprintf(fp, "~j \"%s\"", xformInstance->mpInput->mpMacro->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        WriteXformInstance(fp, binary, xformInstance->mpInput);
      }
    }
  
    if (isHTKCompatible) {
      PutKwd(fp, binary, KID_MMFIDMask);
      fputs("* ", fp);
  
      if (ParmKind2Str(mParamKind, parmkindstr)) 
      {
        fprintf(fp, "<%s> ", parmkindstr);
      }
  
      PutKwd(fp, binary, KID_LinXform);
    }
  
    PutKwd(fp, binary, KID_VecSize);
    PutInt(fp, binary, xformInstance->OutSize());
    PutNLn(fp, binary);
  
    if (xformInstance->mpXform->mpMacro) 
    {
      fprintf(fp, "~x \"%s\"", xformInstance->mpXform->mpMacro->mpName);
      PutNLn(fp, binary);
    } 
    else 
    {
      WriteXform(fp, binary, xformInstance->mpXform);
    }
  }

  
  //***************************************************************************
  //*****************************************************************************  
  void 
  ModelSet::
  WriteXform(FILE *fp, bool binary, Xform *xform)
  {
    XformType type = xform->mXformType;
    
    switch (type)
    {
      case XT_LINEAR    : WriteLinearXform   (fp, binary, static_cast<LinearXform *>(xform)); break;
      case XT_COPY      : WriteCopyXform     (fp, binary, static_cast<CopyXform *>(xform)); break;      
      case XT_BLOCKCOPY : WriteBlockCopyXform(fp, binary, static_cast<BlockCopyXform *>(xform)); break;      
      case XT_TRANSPOSE : WriteTransposeXform(fp, binary, static_cast<TransposeXform *>(xform)); break;
      case XT_WINDOW    : WriteWindowXform   (fp, binary, static_cast<WindowXform *>(xform)); break;
      case XT_GMM_POSTERIORS:  WriteGmmPosteriorsXform (fp, binary, static_cast<GmmPosteriorsXform  *>(xform)); break;
      case XT_FEATURE_MAPPING: WriteFeatureMappingXform(fp, binary, static_cast<FeatureMappingXform *>(xform)); break;
      case XT_FRANTA_PRODUCT : WriteFrantaProductXform (fp, binary, static_cast<FrantaProductXform  *>(xform)); break;
      case XT_FUNC      : WriteFuncXform     (fp, binary, static_cast<FuncXform *>(xform)); break;
      case XT_BIAS      : WriteBiasXform     (fp, binary, static_cast<BiasXform *>(xform)); break;
      case XT_STACKING  : WriteStackingXform (fp, binary, static_cast<StackingXform *>(xform)); break;
      case XT_CONSTANT  : WriteConstantXform (fp, binary, static_cast<ConstantXform *>(xform)); break;
      case XT_REGIONDEPENDENT : WriteRegionDependentXform(fp, binary, static_cast<RegionDependentXform *>(xform)); break;
      case XT_COMPOSITE : WriteCompositeXform(fp, binary, static_cast<CompositeXform *>(xform)); break;
      default:  break;    
    }
  } //WriteXform(FILE *fp, bool binary, Xform *xform)

  
  //***************************************************************************
  //*****************************************************************************  
  void 
  ModelSet::
  WriteCompositeXform(FILE *fp, bool binary, CompositeXform *xform)
  {
    size_t          i;
    size_t          j;
    bool            isHTKCompatible = true;
  
    if (xform->mpMacro || xform->mNLayers != 1) 
    {
      isHTKCompatible = 0;
    } 
    else
    {
      for (i = 0; i < xform->mpLayer[0].mNBlocks; i++) 
      {
        if (xform->mpLayer[0].mpBlock[i]->mXformType != XT_LINEAR) 
        {
          isHTKCompatible = 0;
          break;
        }
      }
    }
    
    if (xform->mNLayers > 1) 
    {
      PutKwd(fp, binary, KID_NumLayers);
      PutInt(fp, binary, xform->mNLayers);
      PutNLn(fp, binary);
    }
  
    for (i=0; i < xform->mNLayers; i++) 
    {
      if (xform->mNLayers > 1) 
      {
        PutKwd(fp, binary, KID_Layer);
        PutInt(fp, binary, i+1);
        PutNLn(fp, binary);
      }
  
      if (isHTKCompatible) 
      {
        PutKwd(fp, binary, KID_BlockInfo);
        PutInt(fp, binary, xform->mpLayer[i].mNBlocks);
        PutNLn(fp, binary);
  
        for (j = 0; j < xform->mpLayer[i].mNBlocks; j++) 
        {
          PutInt(fp, binary, xform->mpLayer[i].mpBlock[j]->mOutSize);
        }
  
        PutNLn(fp, binary);
      } 
      else if (xform->mpLayer[i].mNBlocks > 1) 
      {
        PutKwd(fp, binary, KID_NumBlocks);
        PutInt(fp, binary, xform->mpLayer[i].mNBlocks);
        PutNLn(fp, binary);
      }
  
      for (j = 0; j < xform->mpLayer[i].mNBlocks; j++) {
        if (isHTKCompatible || xform->mpLayer[i].mNBlocks > 1) {
          PutKwd(fp, binary, KID_Block);
          PutInt(fp, binary, j+1);
          PutNLn(fp, binary);
        }
        if (xform->mpLayer[i].mpBlock[j]->mpMacro) {
          fprintf(fp, "~x \"%s\"", xform->mpLayer[i].mpBlock[j]->mpMacro->mpName);
          PutNLn(fp, binary);
        } else {
          WriteXform(fp, binary, xform->mpLayer[i].mpBlock[j]);
        }
      }
    }
  } //WriteCompositeXform(FILE *fp, bool binary, CompositeXform *xform)

  
  //***************************************************************************
  //*****************************************************************************  
  void 
  ModelSet::
  WriteRegionDependentXform(FILE *fp, bool binary, RegionDependentXform *xform)
  {
    size_t          j;
  
    PutKwd(fp, binary, KID_RegionDependent);
    PutNLn(fp, binary);
  
    if (xform->mNBlocks > 1) 
    {
      PutKwd(fp, binary, KID_NumBlocks);
      PutInt(fp, binary, xform->mNBlocks);
      PutNLn(fp, binary);
    }
  
    for (j = 0; j < xform->mNBlocks; j++) {
      if (xform->mNBlocks > 1) {
        PutKwd(fp, binary, KID_Block);
        PutInt(fp, binary, j+1);
        PutNLn(fp, binary);
      }
      if (xform->mpBlock[j]->mpMacro) {
        fprintf(fp, "~x \"%s\"", xform->mpBlock[j]->mpMacro->mpName);
        PutNLn(fp, binary);
      } else {
        WriteXform(fp, binary, xform->mpBlock[j]);
      }
    }
  } //WriteCompositeXform(FILE *fp, bool binary, CompositeXform *xform)


  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteFuncXform(FILE *fp, bool binary, FuncXform *xform)
  {
    PutKwd(fp, binary, gFuncTable[xform->mFuncId].KID);
    PutInt(fp, binary, xform->mOutSize);
    PutNLn(fp, binary);
  } //WriteFuncXform(FILE *fp, bool binary, FuncXform *xform)


  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteBiasXform(FILE *fp, bool binary, BiasXform* xform)
  {
    size_t  i;
    
    if(xform->mUsePredefVector)  
    { 
      PutKwd(fp, binary, KID_BiasPredef);
      PutInt(fp, binary, xform->mOutSize);
      
      xform->mPredefVector.Write(fp, binary, this);               
    }
    else            
    {      
      PutKwd(fp, binary, KID_Bias);
      PutInt(fp, binary, xform->mOutSize);
      PutNLn(fp, binary);
       
      // save raw data
      for (i=0; i < xform->mOutSize; i++) 
      {
        PutFlt(fp, binary, xform->mVector[0][i]);
      }
      PutNLn(fp, binary);
    } 
  } //WriteBiasXform(FILE *fp, bool binary, BiasXform *xform)

  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteWindowXform(FILE *fp, bool binary, WindowXform* xform)
  {
    size_t  i;
    
    if(xform->mUsePredefVector)  
    { 
      PutKwd(fp, binary, KID_WindowPredef);
      PutInt(fp, binary, xform->mOutSize);
      
      xform->mPredefVector.Write(fp, binary, this);
    }
    else            
    {      
      PutKwd(fp, binary, KID_Window);
      PutInt(fp, binary, xform->mOutSize);
      PutNLn(fp, binary);
      
      // save raw data
      for (i=0; i < xform->mOutSize; i++) 
      {
        PutFlt(fp, binary, xform->mVector[i]);
      }
      PutNLn(fp, binary);
    } 
  } //WriteWindowXform(FILE *fp, bool binary, WindowXform *xform)

  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteLinearXform(FILE *fp, bool binary, LinearXform *xform)
  {
    size_t  i;
    size_t  j;
    
    if(xform->mPredefinedID != PLTID_NONE)
    {    
      PutKwd(fp, binary, KID_XformPredef);
      PutInt(fp, binary, xform->mOutSize);
      PutInt(fp, binary, xform->mInSize);
      // PutNLn(fp, binary);

      switch(xform->mPredefinedID)
      {
        case PLTID_DCT:
          PutString(fp, binary, "Dct");
          PutSpace(fp, binary);
          PutInt(fp, binary, (xform->mIncludeC0 ? 1 : 0));
          PutNLn(fp, binary);
          break;

        case PLTID_CONST:
          PutString(fp, binary, "Const");
          PutSpace(fp, binary);
          PutFlt(fp, binary, xform->mConstant);
          PutNLn(fp, binary);
          break;

        case PLTID_DIAG:
          PutString(fp, binary, "Diag");
          PutSpace(fp, binary);
          PutFlt(fp, binary, xform->mConstant);
          PutNLn(fp, binary);
          break;
  
        case PLTID_RANDOM:
          PutString(fp, binary, "Random");
          PutSpace(fp, binary);  
          PutFlt(fp, binary, xform->mMinValue);
          PutFlt(fp, binary, xform->mMaxValue);
          PutInt(fp, binary, static_cast<int>(xform->mSeed));
          PutNLn(fp, binary);
          break;
      }       
    }
    else
    {
      PutKwd(fp, binary, KID_Xform);
      PutInt(fp, binary, xform->mOutSize);
      PutInt(fp, binary, xform->mInSize);
      PutNLn(fp, binary);    
       
      // save raw data
      for (i=0; i < xform->mOutSize; i++) 
      {
        for (j=0; j < xform->mInSize; j++) 
        {
//        if (mUseNewMatrix)
//        {
            PutFlt(fp, binary, xform->mMatrix[i][j]);
//        }
//        else
//        {
//          PutFlt(fp, binary, xform->mpMatrixO[i * xform->mInSize + j]);
//        }
        }
        PutNLn(fp, binary);
      }       
    }   
  } //WriteLinearXform(FILE *fp, bool binary, LinearXform *xform)
  

  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteStackingXform(FILE *fp, bool binary, StackingXform *xform)
  {
    PutKwd(fp, binary, KID_Stacking);
    PutInt(fp, binary, xform->mOutSize / xform->mInSize);
    PutInt(fp, binary, xform->mInSize);
    PutNLn(fp, binary);
  } //WriteStackingXform(FILE *fp, bool binary, StackingXform *xform)
  
  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteCopyXform(FILE *fp, bool binary, CopyXform *xform)
  {
    size_t      i;
    size_t      j;
    int         step = 0;
    int *       ids = xform->mpIndices;
  
    fprintf(fp, "<Copy> %d %d\n", (int) xform->mOutSize, (int) xform->mInSize);
  
    for (i=0; i < xform->mOutSize; i++) 
    {
      if (i + 1 < xform->mOutSize) {
        step = ids[i+1] - ids[i];
      }
  
      for (j = i + 2; j < xform->mOutSize && ids[j] - ids[j - 1] == step; j++)
        ;
  
      if (step == 1 && j > i + 2) 
      {
        fprintf(fp, " %d:%d",    ids[i] + 1, ids[j - 1]+1);
        i = j - 1;
      } 
      else if (j > i + 3) 
      {
        fprintf(fp, " %d:%d:%d", ids[i]+1, step, ids[j-1]+1);
        i = j-1;
      } 
      else 
      {
        fprintf(fp, " %d", ids[i]+1);
      }      
    }
    
    fputs("\n", fp);
  } //WriteCopyXform(FILE *fp, bool binary, CopyXform *xform)

  void 
  ModelSet::
  WriteConstantXform(FILE *fp, bool binary, ConstantXform* xform)
  {
    size_t  i;
  
    PutKwd(fp, binary, KID_Constant);
    PutInt(fp, binary, xform->mOutSize);
    PutNLn(fp, binary);
    for (i=0; i < xform->mOutSize; i++) 
    {
      PutFlt(fp, binary, xform->mVector[0][i]);
    }
    PutNLn(fp, binary);
  } //WriteConstantXform(FILE *fp, bool binary, BiasXform *xform)
  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteBlockCopyXform(FILE *fp, bool binary, BlockCopyXform *xform)
  {
    size_t      i;
    //size_t      j;
    //int         step = 0;
    //int *       ids = xform->mpIndices;
  
    BlockCopyXform::Block *p_blocks = xform->mpBlocks;    
    
    if(!binary)
    {
      fprintf(fp, "<BlockCopy> %d %d %d\n", (int) xform->mNBlocks, (int)xform->mNRows, (int)(xform->mInSize / xform->mNRows));
    
      for (i=0; i < xform->mNBlocks; i++) 
      {
        fputs(" ", fp);
        if(p_blocks[i].mFromR != 1 || p_blocks[i].mToR != 1)  // ignore row range if 1:1:1
        {
          if(p_blocks[i].mStepR == 1)                         // ignore step if step == 1
          {
            if(p_blocks[i].mFromR == p_blocks[i].mToR)        // ignore range end if if from == to
            {
              fprintf(fp, "%d,", p_blocks[i].mFromR);
            }
            else
            {
              fprintf(fp, "%d:%d,", p_blocks[i].mFromR, p_blocks[i].mToR);          
            }
          }
          else
          {
            fprintf(fp, "%d:%d:%d,", p_blocks[i].mFromR, p_blocks[i].mStepR, p_blocks[i].mToR);      
          }
        }
        if(p_blocks[i].mStepC == 1)                          // ignore step if step == 1
        {
          if(p_blocks[i].mFromC == p_blocks[i].mToC)         // ignore range end if if from == to
          {
            fprintf(fp, "%d", p_blocks[i].mFromC);
          }
          else
          {
            fprintf(fp, "%d:%d", p_blocks[i].mFromC, p_blocks[i].mToC);        
          }
        }
        else
        {
          fprintf(fp, "%d:%d:%d", p_blocks[i].mFromC, p_blocks[i].mStepC, p_blocks[i].mToC);      
        }
      }
    
      fputs("\n", fp);
    }
    else
    {
      PutKwd(fp, binary, KID_BlockCopy);
      
      PutInt(fp, binary, xform->mNBlocks);
      PutInt(fp, binary, xform->mNRows);
      PutInt(fp, binary, xform->mInSize / xform->mNRows);
      
      for (i=0; i < xform->mNBlocks; i++)
      {
        PutInt(fp, binary, p_blocks[i].mFromR);
        PutInt(fp, binary, p_blocks[i].mStepR);
        PutInt(fp, binary, p_blocks[i].mToR);
        PutInt(fp, binary, p_blocks[i].mFromC);
        PutInt(fp, binary, p_blocks[i].mStepC);
        PutInt(fp, binary, p_blocks[i].mToC);        
      }      
    }
  } //WriteCopyXform(FILE *fp, bool binary, CopyXform *xform)

  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteTransposeXform(FILE *fp, bool binary, TransposeXform *xform)
  {
    PutKwd(fp, binary, KID_Transpose);
    PutInt(fp, binary, xform->mInRows);
    PutInt(fp, binary, xform->mOutSize / xform->mInRows);
    PutNLn(fp, binary);
  } //WriteTransposeXform(FILE *fp, bool binary, CopyXform *xform)

  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteFeatureMappingXform(FILE* fp, bool binary, FeatureMappingXform* xform)
  {
    fprintf(fp, "<ExtendedXform> FeatureMapping %d\n", (int) xform->mInSize);
    
    fprintf(fp, "<State> 1\n");
    WriteState(fp, binary, xform->mpStateFrom);
    fputs("\n", fp);
    
    fprintf(fp, "<State> 2\n");
    WriteState(fp, binary, xform->mpStateTo);
    fputs("\n", fp);
  } //WriteFeatureMappingXform(FILE *fp, bool binary, FeatureMappingXform *xform)

  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteGmmPosteriorsXform(FILE* fp, bool binary, GmmPosteriorsXform* xform)
  {
    PutKwd(fp, binary, KID_GmmPosteriors);
    PutInt(fp, binary, xform->mInSize);
    PutInt(fp, binary, xform->mNBestLikesVector.size());
    PutFlt(fp, binary, xform->mScale);
    PutNLn(fp, binary);
    WriteState(fp, binary, xform->mpState);
  } //WriteGmmPosteriorsXform(FILE *fp, bool binary, GmmPosteriorsXform *xform)

  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteFrantaProductXform(FILE* fp, bool binary, FrantaProductXform* xform)
  {
    fprintf(fp, "<ExtendedXform> FrantaProduct %d %d\n", 
      (int) xform->mInSize, (int) xform->NParts());
  
    fputs("\n", fp);
  } //WriteFrantaProductXform(FILE *fp, bool binary, FrantaProductXform *xform)

  
  //###########################################################################################################
  
  //###########################################################################
  //###########################################################################
  // XFORM LIST INPUT
  //###########################################################################
  //###########################################################################
    
  
  //***************************************************************************
  //*****************************************************************************  
  void 
  ModelSet::
  ReadXformList(const char * pFileName)
  {
    char      line[1024];
    FILE *    fp;
    size_t    nlines=0;
  
    if ((fp = fopen(pFileName, "rt")) == NULL)
        Error("ReadXformList: Cannot open file: '%s'", pFileName);
  
    while (fgets(line, sizeof(line), fp)) 
    {
      char*             xformName;
      char*             makeXformShellCommand = NULL;
      char              termChar = '\0';
      size_t            i = strlen(line);
      Macro*            macro;
      MakeXformCommand* mxfc;
  
      nlines++;
      
      if (line[i-1] != '\n' && getc(fp) != EOF) 
      {
        Error("ReadXformList: Line %d is too long in file: %s",
              (int) nlines, pFileName);
      }
  
      for (; i > 0 && isspace(line[i-1]); i--) 
        line[i-1] = '\0';
  
      if (i == 0) 
        continue;
  
      for (i = 0; isspace(line[i]); i++)
        ;
  
      if (line[i] == '\"' || line[i] == '\'')
        termChar = line[i++];
  
      xformName = &line[i];
      
      while (line[i] != '\0' && line[i] != termChar &&
            (termChar != '\0' || !isspace(line[i]))) 
        i++;
  
      if (termChar != '\0') 
      {
        if (line[i] != termChar) 
        {
          Error("ReadXformList: Terminanting %c expected at line %d in file %s",
                termChar, (int) nlines, pFileName);
        }
        
        line[i++] = '\0';
      }
  
      if (line[i] != '\0') 
      { // shell command follows
        for (; isspace(line[i]); i++) 
          line[i] = '\0';
        makeXformShellCommand = &line[i];
      }
  
      macro = FindMacro(&mXformHash, xformName);
      
      if (macro == NULL) 
      {
        Error("ReadXformList: Undefined Xform '%s' at line %d in file %s",
              xformName, (int) nlines, pFileName);
      }
  
      mpXformToUpdate = (MakeXformCommand*)
        realloc(mpXformToUpdate,
                sizeof(MakeXformCommand) * ++mNumberOfXformsToUpdate);
  
      if (mpXformToUpdate == NULL) 
      {
        Error("ReadXformList: Insufficient memory");
      }
  
      mxfc = &mpXformToUpdate[mNumberOfXformsToUpdate-1];
      mxfc->mpXform = (Xform *) macro->mpData;
      mxfc->mpShellCommand = NULL;
  
      if (makeXformShellCommand) 
      {
        if ((mxfc->mpShellCommand = strdup(makeXformShellCommand)) == NULL) 
        {
          Error("ReadXformList: Insufficient memory");
        }
      }
    }
    fclose(fp);
  }; //ReadXformList(const std::string & rFileName);
  

  //***************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteClusterWeightsVector(size_t i)
  {
    char buff[1024];
    
    MakeFileName(buff, mpClusterWeightVectors[i]->mpMacro->mpFileName, 
      mClusterWeightsOutPath.c_str(), NULL);
    
    if (!mClusterWeightsStream.is_open())
    {
      mClusterWeightsStream.open(buff);
    }
    else if (mClusterWeightsStream.name() != buff)
    {
      mClusterWeightsStream.close();
      mClusterWeightsStream.open(buff);
    }
    
    // check if nothing's wrong with the stream
    if (!mClusterWeightsStream.good())
    {
      Error("Cannot open MMF output file %s", buff);
    }
    
    fprintf(mClusterWeightsStream.file(), "~x \"%s\"\n", 
      mpClusterWeightVectors[i]->mpMacro->mpName);
    WriteBiasXform(mClusterWeightsStream.file(), false, mpClusterWeightVectors[i]);

  }

  
  //###########################################################################################################
  //###########################################################################
  //###########################################################################
  // Accumulator INPUT
  //###########################################################################
  //###########################################################################
  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ReadAccum(int macro_type, HMMSetNodeName nodeName,
                  MacroData * pData, void *pUserData) 
  {
    unsigned int        i;
    unsigned int        j;
    int                 c;
    unsigned int        size   = 0;
    FLOAT*              vector = NULL;
    Macro*              macro;
    ReadAccumUserData*  ud = (ReadAccumUserData *) pUserData;
    //  FILE *fp =        ((ReadAccumUserData *) pUserData)->fp;
    //  char *fn =        ((ReadAccumUserData *) pUserData)->fn;
    //  ModelSet *hmm_set = ((ReadAccumUserData *) pUserData)->hmm_set;
    //  float weight =    ((ReadAccumUserData *) pUserData)->weight;
    char                xfName[128];
  
    xfName[sizeof(xfName)-1] = '\0';

    if (macro_type == mt_mixture)  
    {
      Mixture* mixture = static_cast<Mixture*>(pData);
      
      //************************************************************************
      // cluster adaptive training - cluster parameters accumulators 
      // we expect all accums to be initialized if mAccumG is initialized
      if (mixture->mAccumG.IsInitialized())
      {
        size_t r;     // row counter
        size_t rows;
        
        rows = mixture->mAccumG.Rows();
        for (r = 0; r < rows; r++)
        {                                                              
          if (faddfloat(mixture->mAccumG[r], mixture->mAccumG.Cols(), 1, ud->mpFp) 
              != mixture->mAccumG.Cols())
          {
            Error("G: Incompatible accumulator file: '%s'", ud->mpFileName);
          }
        }

        rows = mixture->mAccumK.Rows();
        for (r = 0; r < rows; r++)
        {
          if (faddfloat(mixture->mAccumK[r], mixture->mAccumK.Cols(), 1, ud->mpFp) 
              != mixture->mAccumK.Cols())
          {
            Error("K: Incompatible accumulator file: '%s'", ud->mpFileName);
          }
        }

        if (faddfloat(mixture->mAccumL.pData(), mixture->mAccumL.Length(), 1, ud->mpFp) 
            != mixture->mAccumL.Length())
        {
          Error("L: Incompatible accumulator file: '%s'", ud->mpFileName);
        }
        
      } // if (mixture->mAccumG.IsInitialized())
      //************************************************************************
      
    } // if (macro_type == mt_mixture)  
    
    else if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      XformStatAccum *  xfsa = NULL;
      //UINT_32           size_inf;
      size_t            size_inf;
      UINT_32           nxfsa_inf;
      size_t            nxfsa = 0;
  
      if (macro_type == mt_mean) 
      {
        xfsa   = ((Mean*)pData)->mpXformStatAccum;
        nxfsa  = ((Mean*)pData)->mNumberOfXformStatAccums;
        size   = ((Mean*)pData)->VectorSize();
        //vector = ((Mean *)pData)->mpVectorO+size;
        vector = ((Mean*)pData)->mpAccums;
        size   = size + 1;
      } 
      else if (macro_type == mt_variance) 
      {
        xfsa   = ((Variance*)pData)->mpXformStatAccum;
        nxfsa  = ((Variance*)pData)->mNumberOfXformStatAccums;
        size   = ((Variance*)pData)->VectorSize();
        //vector = ((Variance *)pData)->mpVectorO+size;
        vector = ((Variance*)pData)->mpAccums;
        size   = size * 2 + 1;
      }
  
      if (ud->mMmi) 
        vector += size;
  
      if (faddfloat(vector, size, ud->mWeight,    ud->mpFp) != size ||
          fread(&nxfsa_inf, sizeof(nxfsa_inf), 1, ud->mpFp) != 1) 
      {
        Error("Incompatible accumulator file: '%s'", ud->mpFileName);
      }
  
      if (!ud->mMmi && !ud->mpModelSet->mCmllrStats ) 
      { // MMI estimation of Xform statistics has not been implemented yet
        for (i = 0; i < nxfsa_inf; i++) 
        {
          if (getc(ud->mpFp) != '"') 
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
          for (j=0; (c=getc(ud->mpFp)) != EOF && c != '"' && j < sizeof(xfName)-1; j++) 
            xfName[j] = c;
  
          xfName[j] = '\0';
          
          if (c == EOF)
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
          macro = FindMacro(&ud->mpModelSet->mXformHash, xfName);
  
          if (fread(&size_inf, sizeof(size_inf), 1, ud->mpFp) != 1) 
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
          if (macro != NULL) 
          {
            size = ((LinearXform *) macro->mpData)->mInSize;
            size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
  
            if (size != size_inf)
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
            for (j = 0; j < nxfsa && xfsa[j].mpXform != macro->mpData; j++)
            {}
            
            if (j < nxfsa) 
            {
              if (faddfloat(xfsa[j].mpStats, size, ud->mWeight, ud->mpFp) != size  ||
                  faddfloat(&xfsa[j].mNorm,     1, ud->mWeight, ud->mpFp) != 1) 
              {
                Error("Invalid accumulator file: '%s'", ud->mpFileName);
              }
            } 
            else 
            {
              macro = NULL;
            }
          }
  
          // Skip Xform accumulator
          if (macro == NULL) 
          { 
            FLOAT f;
            for (j = 0; j < size_inf+1; j++) 
              fread(&f, sizeof(f), 1, ud->mpFp);
          }
        }
      }
    } 
    
    else if (macro_type == mt_state) 
    {
      State *state = (State *) pData;
      if (state->mOutPdfKind == KID_DiagC) {
        FLOAT junk;
        for (i = 0; i < state->mNMixtures; i++) {
          if (ud->mMmi == 1) {
            if (faddfloat(&state->mpMixture[i].mWeightAccumDen, 1, ud->mWeight, ud->mpFp) != 1 ||
                faddfloat(&junk,                                1, ud->mWeight, ud->mpFp) != 1) {
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
            }
          } else if (ud->mMmi == 2) {
            if (faddfloat(&junk,                                1, ud->mWeight, ud->mpFp) != 1 ||
                faddfloat(&state->mpMixture[i].mWeightAccumDen, 1, ud->mWeight, ud->mpFp) != 1) {
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
            }
          } else {
            if (faddfloat(&state->mpMixture[i].mWeightAccum,     1, ud->mWeight, ud->mpFp) != 1 ||
                faddfloat(&state->mpMixture[i].mWeightAccumDen,  1, ud->mWeight, ud->mpFp) != 1) {
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
            }
          }
        }
      }
    } else if (macro_type == mt_transition) {
      FLOAT f;
      size   = SQR(((Transition *) pData)->mNStates);
      vector =     ((Transition *) pData)->mpMatrixO + size;
  
      for (i = 0; i < size; i++) {
        if (fread(&f, sizeof(FLOAT), 1, ud->mpFp) != 1) {
        Error("Incompatible accumulator file: '%s'", ud->mpFileName);
        }
        if (!ud->mMmi) { // MMI estimation of transition probabilities has not been implemented yet
          f += my_log(ud->mWeight);
          LOG_INC(vector[i], f);
        }
      }
    } else if (macro_type == mt_Xform) {
      if (ud->mpModelSet->mCmllrStats 
      && ((Xform *) pData)->mpCmllrStats != NULL) {
        size = ((Xform *) pData)->mInSize;
        size = (size + size*(size+1)/2) * (size -1) + 1;
        vector = ((Xform *) pData)->mpCmllrStats;
        
        if (faddfloat(vector, size, ud->mWeight,    ud->mpFp) != size) {
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
        }    
      }
    }
  } // void ReadAccum(int macro_type, HMMSetNodeName nodeName...)
  

  //****************************************************************************  
  //****************************************************************************  
  void
  ModelSet::
  ReadAccums(const FileListElem& rFile, long* totFrames, FLOAT* totLogLike, 
      int mmiDenominatorAccums)
  {
    IStkStream                in;
    FILE*                     fp;
    char                      macro_name[128];
    MyHSearchData*            hash;
    unsigned int              i;
    int                       t = 0;
    int                       c;
    int                       skip_accum = 0;
    INT_32                    occurances;
    ReadAccumUserData         ud;
    Macro*                    macro;
    int                       mtm = MTM_PRESCAN | 
                                    MTM_STATE |
                                    MTM_MIXTURE |
                                    MTM_MEAN |   
                                    MTM_VARIANCE | 
                                    MTM_TRANSITION;
  
    macro_name[sizeof(macro_name)-1] = '\0';
  
    // open the file
    in.open(rFile.Physical().c_str(), std::ios::binary);
    if (!in.good())
    {
      Error("Cannot open input accumulator file: '%s'", rFile.Physical().c_str());
    }
    
    fp = in.file();
    
    INT_32 i32;
    if (fread(&i32,       sizeof(i32),  1, fp) != 1 ||
        fread(totLogLike, sizeof(FLOAT), 1, fp) != 1) 
    {
      Error("Invalid accumulator file: '%s'", rFile.Physical().c_str());
    }

    *totFrames = i32;
    
//    *totFrames  *= weight; // Not sure whether we should report weighted quantities or not
//    *totLogLike *= weight;
  
    ud.mpFileName   = rFile.Physical().c_str();
    ud.mpFp         = fp;
    ud.mpModelSet   = this;
    ud.mWeight      = rFile.Weight();;
    ud.mMmi         = mmiDenominatorAccums;
  
    for (;;) {
      // Skip to the begining of the next macro accumulator
      if (skip_accum) { 
        for (;;) {
          // read until '~'
          while ((c = getc(fp)) != '~' && c != EOF)
          { }
          
          if (c == EOF) {
            break;
          }
            
          if (strchr("hsmuvtx", t = c = getc(fp)) &&
             (c = getc(fp)) == ' ' && (c = getc(fp)) == '"')
          {  
            break;
          }
          
          ungetc(c, fp);
        }
        
        if (c == EOF) {
          break;
        }
      } 
      else {
        if ((c = getc(fp)) == EOF) {
          break;
        }
        
        if (c != '~'      || !strchr("hsmuvtx", t = getc(fp)) ||
          getc(fp) != ' ' || getc(fp) != '"') 
        {
          Error("Incompatible accumulator file: '%s'", rFile.Physical().c_str());
        }
      }
  
      for (i=0; (c = getc(fp))!=EOF && c!='"' && i<sizeof(macro_name)-1; i++) {
        macro_name[i] = c;
      }

      macro_name[i] = '\0';
  
      hash = t == 'h' ? &mHmmHash :
             t == 's' ? &mStateHash :
             t == 'm' ? &mMixtureHash :
             t == 'u' ? &mMeanHash :
             t == 'v' ? &mVarianceHash :
             t == 't' ? &mTransitionHash : 
             t == 'x' ? &mXformHash : NULL;
  
      assert(hash);
      if ((macro = FindMacro(hash, macro_name)) == NULL) 
      {
        skip_accum = 1;
        continue;
      }
  
      skip_accum = 0;
      if (fread(&occurances, sizeof(occurances), 1, fp) != 1) 
      {
        Error("Invalid accumulator file: '%s'", rFile.Physical().c_str());
      }
      
      if (!mmiDenominatorAccums) macro->mOccurances += occurances;
      switch (t) 
      {
        case 'h': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
        case 's': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
        case 'm': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
        case 'u': ReadAccum(mt_mean,       NULL, macro->mpData, &ud);          break;
        case 'v': ReadAccum(mt_variance,   NULL, macro->mpData, &ud);          break;
        case 't': ReadAccum(mt_transition, NULL, macro->mpData, &ud);          break;
        case 'x': ReadAccum(mt_Xform,      NULL, macro->mpData, &ud);          break;
        default:  assert(0);
      }
    }
    
    in.close();
    //delete [] ud.mpFileName;
    //free(ud.mpFileName);    
  }; // ReadAccums(...)



  //****************************************************************************  
  //****************************************************************************  
  //  void
  //  ModelSet::
  //  ReadAccums(const char * pFileName, 
  //             float        weight,
  //             long *       totFrames, 
  //             FLOAT *      totLogLike, 
  //             int          mmiDenominatorAccums)
  //  {
  //    IStkStream                in;
  //    FILE*                     fp;
  //    char                      macro_name[128];
  //    MyHSearchData*            hash;
  //    unsigned int              i;
  //    int                       t = 0;
  //    int                       c;
  //    int                       skip_accum = 0;
  //    INT_32                    occurances;
  //    ReadAccumUserData         ud;
  //    Macro*                    macro;
  //    int                       mtm = MTM_PRESCAN | 
  //                                    MTM_STATE |
  //                                    MTM_MIXTURE |
  //                                    MTM_MEAN |   
  //                                    MTM_VARIANCE | 
  //                                    MTM_TRANSITION;
  //  
  //    macro_name[sizeof(macro_name)-1] = '\0';
  //  
  //    // open the file
  //    in.open(pFileName, ios::binary);
  //    if (!in.good())
  //    {
  //      Error("Cannot open input accumulator file: '%s'", pFileName);
  //    }
  //    
  //    fp = in.file();
  //    
  //    INT_32 i32;
  //    if (fread(&i32,       sizeof(i32),  1, fp) != 1 ||
  //        fread(totLogLike, sizeof(FLOAT), 1, fp) != 1) 
  //    {
  //      Error("Invalid accumulator file: '%s'", pFileName);
  //    }
  //    *totFrames = i32;
  //    
////      *totFrames  *= weight; // Not sure whether we should report weighted quantities or not
////      *totLogLike *= weight;
  //  
  //    ud.mpFileName   = pFileName;
  //    ud.mpFp         = fp;
  //    ud.mpModelSet   = this;
  //    ud.mWeight      = weight;
  //    ud.mMmi         = mmiDenominatorAccums;
  //  
  //    for (;;) 
  //    {
  //      if (skip_accum) 
  //      { // Skip to the begining of the next macro accumulator
  //        for (;;) 
  //        {
  //          while ((c = getc(fp)) != '~' && c != EOF)
  //            ;
  //          
  //          if (c == EOF) 
  //            break;
  //            
  //          if (strchr("hsmuvtx", t = c = getc(fp)) &&
  //            (c = getc(fp)) == ' ' && (c = getc(fp)) == '"')
  //          {  
  //            break;
  //          }
  //          
  //          ungetc(c, fp);
  //        }
  //        
  //        if (c == EOF) 
  //          break;
  //      } 
  //      else 
  //      {
  //        if ((c = getc(fp)) == EOF) break;
  //        
  //        if (c != '~'      || !strchr("hsmuvtx", t = getc(fp)) ||
  //          getc(fp) != ' ' || getc(fp) != '"') 
  //        {
  //          Error("Incomatible accumulator file: '%s'", pFileName);
  //        }
  //      }
  //  
  //      for (i=0; (c = getc(fp))!=EOF && c!='"' && i<sizeof(macro_name)-1; i++) 
  //      {
  //        macro_name[i] = c;
  //      }
  //      macro_name[i] = '\0';
  //  
  //      hash = t == 'h' ? &mHmmHash :
  //             t == 's' ? &mStateHash :
  //             t == 'm' ? &mMixtureHash :
  //             t == 'u' ? &mMeanHash :
  //             t == 'v' ? &mVarianceHash :
  //             t == 't' ? &mTransitionHash : 
  //             t == 'x' ? &mXformHash : NULL;
  //  
  //      assert(hash);
  //      if ((macro = FindMacro(hash, macro_name)) == NULL) 
  //      {
  //        skip_accum = 1;
  //        continue;
  //      }
  //  
  //      skip_accum = 0;
  //      if (fread(&occurances, sizeof(occurances), 1, fp) != 1) 
  //      {
  //        Error("Invalid accumulator file: '%s'", pFileName);
  //      }
  //      
  //      if (!mmiDenominatorAccums) macro->mOccurances += occurances;
  //      switch (t) 
  //      {
  //        case 'h': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
  //        case 's': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
  //        case 'm': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
  //        case 'u': ReadAccum(mt_mean,       NULL, macro->mpData, &ud); break;
  //        case 'v': ReadAccum(mt_variance,   NULL, macro->mpData, &ud); break;
  //        case 't': ReadAccum(mt_transition, NULL, macro->mpData, &ud); break;
  //        case 'x': ReadAccum(mt_Xform,      NULL, macro->mpData, &ud); break;
  //        default:  assert(0);
  //      }
  //    }
  //    
  //    in.close();
  //    //delete [] ud.mpFileName;
  //    //free(ud.mpFileName);    
  //  }; // ReadAccums(...)
  
    
  //############################################################################
  //############################################################################
  //############################################################################
  // Accumulator OUTPUT
  //############################################################################
  //############################################################################
  
  //****************************************************************************
  //****************************************************************************
  void 
  WriteAccum(int macro_type, HMMSetNodeName nodeName,
             MacroData* pData, void *pUserData) 
  {
    size_t                i;
    size_t                size;
    FLOAT*                vector = NULL;
  //  FILE *fp = (FILE *) pUserData;
    WriteAccumUserData*   ud = static_cast<WriteAccumUserData*>(pUserData);
    Macro*                macro = pData->mpMacro;
  
    if (macro) {
      INT_32 occurances = macro->mOccurances;

      if ( fprintf(ud->mpFp, "~%c \"%s\"", macro->mType, macro->mpName) < 0 
      ||   fwrite(&occurances, sizeof(INT_32), 1, ud->mpFp) != 1)
       //fwrite(&macro->mOccurances, sizeof(macro->mOccurances), 1, ud->mpFp) != 1)) 
      {
        Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
      }
    }
  
    if (macro_type == mt_mixture)
    {
      Mixture* mixture = static_cast<Mixture*>(pData);
      
      // save cluster parameters accumulators
      if (mixture->mAccumG.IsInitialized())
      {
        size_t r;     // row counter
        size_t rows;
        
        // we need to update the accumulators from the per-speaker partial
        // accumulators
        mixture->UpdateClusterParametersAccums();
        
        // the G accumulator
        rows = mixture->mAccumG.Rows();
        for (r = 0; r < rows; r++)
        {
          fwrite(mixture->mAccumG[r], sizeof(FLOAT), mixture->mAccumG.Cols(), 
            ud->mpFp);
        }
        
        // the K accumulator
        rows = mixture->mAccumK.Rows();
        for (r = 0; r < rows; r++)
        {
          fwrite(mixture->mAccumK[r], sizeof(FLOAT), mixture->mAccumK.Cols(), 
            ud->mpFp);
        }
        
        // the L accumulator
        fwrite(mixture->mAccumL.cpData(), sizeof(FLOAT), mixture->mAccumL.Length(), 
          ud->mpFp);
      }
    }
    
    else if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      XformStatAccum*     xfsa = NULL;
      UINT_32             nxfsa = 0;
  
      if (macro_type == mt_mean) 
      {
        xfsa   = ((Mean*)pData)->mpXformStatAccum;
        nxfsa  = ((Mean*)pData)->mNumberOfXformStatAccums;
        size   = ((Mean*)pData)->VectorSize();
        //vector = ((Mean*)pData)->mpVectorO+size;
        vector = ((Mean*)pData)->mpAccums;
        size   = size + 1;
      } 
      else if (macro_type == mt_variance) 
      {
        xfsa   = ((Variance*)pData)->mpXformStatAccum;
        nxfsa  = ((Variance*)pData)->mNumberOfXformStatAccums;
        size   = ((Variance*)pData)->VectorSize();
        //vector = ((Variance*)pData)->mpVectorO+size;
        vector = ((Variance*)pData)->mpAccums;
        size   = size * 2 + 1;
      }
  
  //    if (ud->mMmi) vector += size; // Move to MMI accums, which follows ML accums
  
      if (fwrite(vector, sizeof(FLOAT), size, ud->mpFp) != size ||
          fwrite(&nxfsa, sizeof(nxfsa),    1, ud->mpFp) != 1) 
      {
        Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
      }
  
  //    if (!ud->mMmi) { // MMI estimation of Xform statistics has not been implemented yet
      if(!ud->mpModelSet->mCmllrStats) {
        for (i = 0; i < nxfsa; i++) 
        {
          size = xfsa[i].mpXform->mInSize;
          size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
          assert(xfsa[i].mpXform->mpMacro != NULL);
          if (fprintf(ud->mpFp, "\"%s\"", xfsa[i].mpXform->mpMacro->mpName) < 0 ||
            fwrite(&size,           sizeof(size),        1, ud->mpFp) != 1    ||
            fwrite(xfsa[i].mpStats, sizeof(FLOAT), size, ud->mpFp) != size ||
            fwrite(&xfsa[i].mNorm,  sizeof(FLOAT),    1, ud->mpFp) != 1) 
          {
            Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
          }
        }
      }
  //    }
    }
    
    else if (macro_type == mt_state) 
    {
      State *state = (State *) pData;
      if (state->mOutPdfKind == KID_DiagC) 
      {
        for (i = 0; i < state->mNMixtures; i++) 
        {
          if (fwrite(&state->mpMixture[i].mWeightAccum,
                    sizeof(FLOAT), 1, ud->mpFp) != 1 ||
            fwrite(&state->mpMixture[i].mWeightAccumDen,
                    sizeof(FLOAT), 1, ud->mpFp) != 1) 
          {
            Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
          }
        }
      }
    } 
    
    else if (macro_type == mt_transition) 
    {
      size   = SQR(((Transition *) pData)->mNStates);
      vector = ((Transition *) pData)->mpMatrixO + size;
  
      if (fwrite(vector, sizeof(FLOAT), size, ud->mpFp) != size) 
      {
        Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
      }
    }
    else if (macro_type == mt_Xform) 
    {
      if (ud->mpModelSet->mCmllrStats 
      && ((Xform *) pData)->mpCmllrStats != NULL) {
        size = ((Xform *) pData)->mInSize;
        size = (size + size*(size+1)/2) * (size -1) + 1;
        vector = ((Xform *) pData)->mpCmllrStats;
  
        if (fwrite(vector, sizeof(FLOAT), size, ud->mpFp) != size) {
          Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
        }
      }
    }
  } // WriteAccum(int macro_type, HMMSetNodeName nodeName,..)
  

  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  WriteAccums(const char*  pFileName, 
              const char*  pOutputDir,
              long         totFrames, 
              FLOAT        totLogLike)
  {
    FILE*                 fp;
    OStkStream            output_stream;
    char                  file_name[1024];
    WriteAccumUserData    ud;  
    
    // synthesize new name for the accumulator
    MakeFileName(file_name, pFileName, pOutputDir, NULL);
  
    // open the stream
    output_stream.open(file_name, std::ios::binary);
    if (!output_stream.good()) {
      Error("Cannot open output file: '%s'", file_name);
    }
  
    // TODO: make full use of streams
    fp = output_stream.file();

    // if ((fp = fopen(file_name, "wb")) == NULL) 
    // {
    //   Error("Cannot open output file: '%s'", file_name);
    // }
  
    // write some header info
    INT_32 i32 = totFrames;
    if (fwrite(&i32,  sizeof(i32),  1, fp) != 1 ||
        fwrite(&totLogLike, sizeof(FLOAT), 1, fp) != 1) 
    {
      Error("Cannot write accumulators to file: '%s'", file_name);
    }
  
    ud.mpFp         = fp;
    ud.mpFileName   = file_name;
    ud.mpModelSet   = this;
  //  ud.mMmi = MMI_denominator_accums;
  
    // recursively scan the structure and write accumulator for each element
    Scan(MTM_PRESCAN | (MTM_ALL & ~(MTM_XFORM_INSTANCE)), //|MTM_XFORM
              NULL, WriteAccum, &ud);
  
    // close the stream
    // fclose(fp);
    output_stream.close();
  }; // WriteAccums(...)
  
  
  //###########################################################################
  //###########################################################################
  // PREDEFINED VECTOR
  //###########################################################################
  //##########################################################################
  void PredefinedVector::Read(FILE *fp, size_t vctSize, ModelSet *pModelSet)
  {
    mSize = vctSize;
  
    char *p_win_name = pModelSet->GetString(fp, false);
    size_t i;
    for(i = 0; i < strlen(p_win_name); i++)
    {
      p_win_name[i] = toupper(p_win_name[i]);
    }
    
    if(!strcmp(p_win_name, "HAMMING")) 
    {
      mID = PWID_HAMMING;
      mNRepetitions = pModelSet->GetInt(fp);
            
      size_t win_len = vctSize / mNRepetitions;
      if(win_len * mNRepetitions != vctSize)
      {
        Error("Hamming windows do not fit to the vector size (%s:%d)", pModelSet->mpCurrentMmfName, pModelSet->mCurrentMmfLine);
      }
    }
    else if(!strcmp(p_win_name, "TRIANG")) 
    {
      mID = PWID_TRIANG;
      mNRepetitions = pModelSet->GetInt(fp);

      size_t win_len = vctSize / mNRepetitions;
      if(win_len * mNRepetitions != vctSize)
      {
        Error("Triangular windows do not fit to the vector size (%s:%d)", pModelSet->mpCurrentMmfName, pModelSet->mCurrentMmfLine);
      }
    }
    else if(!strcmp(p_win_name, "CONST")) 
    {
      mID = PWID_CONST;
      mConstant = pModelSet->GetFloat(fp);
    }
    else if(!strcmp(p_win_name, "RANDOM")) 
    {
      mID = PWID_RANDOM;
      mMinValue = pModelSet->GetFloat(fp);
      mMaxValue = pModelSet->GetFloat(fp);
      mSeed = static_cast<unsigned int>(pModelSet->GetInt(fp));
    }
    else if(!strcmp(p_win_name, "LINSPACE")) 
    {
      mID = PWID_LINSPACE;
      mStartValue = pModelSet->GetFloat(fp);
      mEndValue = pModelSet->GetFloat(fp);
      mNRepetitions = pModelSet->GetInt(fp);

      size_t win_len = vctSize / mNRepetitions;
      if(win_len * mNRepetitions != vctSize)
      {
        Error("LinSpace windows do not fit to the vector size (%s:%d)", pModelSet->mpCurrentMmfName, pModelSet->mCurrentMmfLine);
      }
    }
    else 
    {
      Error("Unknown predefined window %s (%s:%d)", p_win_name, pModelSet->mpCurrentMmfName, pModelSet->mCurrentMmfLine);  
    }  
  }
  
  void PredefinedVector::Write(FILE *fp, bool binary, ModelSet *pModelSet)
  {
    switch(mID)
    {
      case PWID_HAMMING:
        pModelSet->PutString(fp, binary, "Hamming");
        pModelSet->PutSpace(fp, binary);
        pModelSet->PutInt(fp, binary, mNRepetitions);
        pModelSet->PutNLn(fp, binary);
        break;
      case PWID_TRIANG:
        pModelSet->PutString(fp, binary, "Triang");
        pModelSet->PutSpace(fp, binary);
        pModelSet->PutInt(fp, binary, mNRepetitions);
        pModelSet->PutNLn(fp, binary);
        break;
      case PWID_CONST:
        pModelSet->PutString(fp, binary, "Const");
        pModelSet->PutSpace(fp, binary);
        pModelSet->PutFlt(fp, binary, mConstant);
        pModelSet->PutNLn(fp, binary);
        break;
      case PWID_RANDOM:
        pModelSet->PutString(fp, binary, "Random");
        pModelSet->PutSpace(fp, binary);
        pModelSet->PutFlt(fp, binary, mMinValue);
        pModelSet->PutFlt(fp, binary, mMaxValue);
        pModelSet->PutInt(fp, binary, static_cast<int>(mSeed));
        pModelSet->PutNLn(fp, binary);
        break;
      case PWID_LINSPACE:
        pModelSet->PutString(fp, binary, "LinSpace");
        pModelSet->PutSpace(fp, binary);
        pModelSet->PutFlt(fp, binary, mStartValue);
        pModelSet->PutFlt(fp, binary, mEndValue);
        pModelSet->PutInt(fp, binary, mNRepetitions);
        pModelSet->PutNLn(fp, binary);
        break;
    }  
  }
        
  void PredefinedVector::Generate(BasicVector<FLOAT> &rVct)
  {    
    size_t win_len;
    
    switch(mID)
    {
      case PWID_HAMMING:
        win_len = mSize / mNRepetitions;    
        GenerateHammingWindow(rVct, win_len, mNRepetitions); 
        break;
      case PWID_TRIANG:
        win_len = mSize / mNRepetitions;          
        GenerateTriangWindow(rVct, win_len, mNRepetitions);
        break;
      case PWID_CONST:
        GenerateConstantWindow(rVct, mSize, mConstant);
        break;
      case PWID_RANDOM:
        GenerateRandomWindow(rVct, mSize, mMinValue, mMaxValue, mSeed);
        break;
      case PWID_LINSPACE:
        win_len = mSize / mNRepetitions;
        GenerateLinSpaceWindow(rVct, win_len, mStartValue, mEndValue, mNRepetitions);
        break;
    }
  }    
    
}; // namespace STK
