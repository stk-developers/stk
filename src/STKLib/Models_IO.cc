#include "Models.h"

namespace STK
{
  int             gCurrentMmfLine       = 1;
  const char *    gpCurrentMmfName;
  bool            gHmmReadBinary;
  bool            gStringUnget          = false;
  char *          gpKwds[KID_MaxKwdID]  = {0};
  
  
  
  //###########################################################################
  //###########################################################################
  // GENERAL FUNCTIONS
  //###########################################################################
  //###########################################################################
  
  //***************************************************************************
  //***************************************************************************
  void InitKwdTable()
  {
    gpKwds[KID_BeginHMM   ] = "BeginHMM";    gpKwds[KID_Use        ] = "Use";
    gpKwds[KID_EndHMM     ] = "EndHMM";      gpKwds[KID_NumMixes   ] = "NumMixes";
    gpKwds[KID_NumStates  ] = "NumStates";   gpKwds[KID_StreamInfo ] = "StreamInfo";
    gpKwds[KID_VecSize    ] = "VecSize";     gpKwds[KID_NullD      ] = "NullD";
    gpKwds[KID_PoissonD   ] = "PoissonD";    gpKwds[KID_GammaD     ] = "GammaD";
    gpKwds[KID_RelD       ] = "RelD";        gpKwds[KID_GenD       ] = "GenD";
    gpKwds[KID_DiagC      ] = "DiagC";       gpKwds[KID_FullC      ] = "FullC";
    gpKwds[KID_XFormC     ] = "XFormC";      gpKwds[KID_State      ] = "State";
    gpKwds[KID_TMix       ] = "TMix";        gpKwds[KID_Mixture    ] = "Mixture";
    gpKwds[KID_Stream     ] = "Stream";      gpKwds[KID_SWeights   ] = "SWeights";
    gpKwds[KID_Mean       ] = "Mean";        gpKwds[KID_Variance   ] = "Variance";
    gpKwds[KID_InvCovar   ] = "InvCovar";    gpKwds[KID_XForm      ] = "XForm";
    gpKwds[KID_GConst     ] = "GConst";      gpKwds[KID_Duration   ] = "Duration";
    gpKwds[KID_InvDiagC   ] = "InvDiagC";    gpKwds[KID_TransP     ] = "TransP";
    gpKwds[KID_DProb      ] = "DProb";       gpKwds[KID_LLTC       ] = "LLTC";
    gpKwds[KID_LLTCovar   ] = "LLTCovar";    gpKwds[KID_XFormKind  ] = "XFormKind";
    gpKwds[KID_ParentXForm] = "ParentXForm"; gpKwds[KID_NumXForms  ] = "NumXForms";
    gpKwds[KID_XFormSet   ] = "XFormSet";    gpKwds[KID_LinXForm   ] = "LinXForm";
    gpKwds[KID_Offset     ] = "Offset";      gpKwds[KID_Bias       ] = "Bias";
    gpKwds[KID_BlockInfo  ] = "BlockInfo";   gpKwds[KID_Block      ] = "Block";
    gpKwds[KID_BaseClass  ] = "BaseClass";   gpKwds[KID_Class      ] = "Class";
    gpKwds[KID_XFormWgtSet] = "XFormWgtSet"; gpKwds[KID_ClassXForm ] = "ClassXForm";
    gpKwds[KID_MMFIDMask  ] = "MMFIDMask";   gpKwds[KID_Parameters ] = "Parameters";
    gpKwds[KID_NumClasses ] = "NumClasses";  gpKwds[KID_AdaptKind  ] = "AdaptKind";
    gpKwds[KID_Prequal    ] = "Prequal";     gpKwds[KID_InputXForm ] = "InputXForm";
    gpKwds[KID_RClass     ] = "RClass";      gpKwds[KID_RegTree    ] = "RegTree";
    gpKwds[KID_Node       ] = "Node";        gpKwds[KID_TNode      ] = "TNode";
    gpKwds[KID_HMMSetID   ] = "HMMSetID";    gpKwds[KID_ParmKind   ] = "ParmKind";
  
    /* Non-HTK keywords */
    gpKwds[KID_FrmExt     ] = "FrmExt";      gpKwds[KID_PDFObsVec  ] = "PDFObsVec";
    gpKwds[KID_ObsCoef    ] = "ObsCoef";     gpKwds[KID_Input      ] = "Input";
    gpKwds[KID_NumLayers  ] = "NumLayers";   gpKwds[KID_NumBlocks  ] = "NumBlocks";
    gpKwds[KID_Layer      ] = "Layer";       gpKwds[KID_Copy       ] = "Copy";
    gpKwds[KID_Stacking   ] = "Stacking";
    /* Numeric functions - FuncXForm*/
    gpKwds[KID_Sigmoid    ] = "Sigmoid";     gpKwds[KID_Log        ] = "Log";
    gpKwds[KID_Exp        ] = "Exp";         gpKwds[KID_Sqrt       ] = "Sqrt";
    gpKwds[KID_SoftMax    ] = "SoftMax";  
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int CheckKwd(const char *str, KeywordID kwdID)
  {
    const char *chptr;
    if (str[0] == ':') {
      gHmmReadBinary = true;
      return str[1] == kwdID;
    }
  
    if (str[0] != '<') return 0;
    for (chptr = gpKwds[kwdID], str++; *chptr; chptr++, str++) {
      if (toupper(*chptr) != toupper(*str)) return 0;
    }
  
    if (str[0] != '>') return 0;
  
    assert(str[1] == '\0');
    gHmmReadBinary = false;
    return 1;
  }
  
  //***************************************************************************
  //***************************************************************************
  char *GetString(FILE *fp, int eofNotExpected)
  {
    static char buffer[1024];
    char ch, *chptr = buffer;
    int lines = 0;
  
  //  fputs("GetString: ", stdout);
    if (gStringUnget) {
      gStringUnget = false;
  
  //    puts(buffer);
      return buffer;
    }
  
    RemoveSpaces(fp);
  
    ch = getc(fp);
    if (ch == '\"' || ch == '\'' ) {
      char termChar = ch;
  
      while (((ch = getc(fp)) != EOF) && 
            (ch != termChar) && 
            ((chptr-buffer) < (sizeof(buffer)-1))) 
      {
        if (ch == '\n') {
          ++lines;
        }
        *chptr++ = ch;
      }
  
      if (ch == EOF && ferror(fp)) {
        Error("Cannot read input file %s", gpCurrentMmfName);
      }
  
      if (ch != termChar) {
        Error("Unterminated string constant (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
      }
      gCurrentMmfLine += lines;
    } else if (ch == '<') {
      *chptr++ = '<';
      while ((ch = getc(fp)) != EOF && !isspace(ch) && ch != '>' && chptr-buffer < sizeof(buffer)-1) {
        *chptr++ = ch;
      }
  
      if (ch == EOF && ferror(fp)) {
        Error("Cannot read input file %s", gpCurrentMmfName);
      }
  
      if (ch != '>') {
        Error("Unterminated keyword %s (%s:%d)", buffer, gpCurrentMmfName, gCurrentMmfLine);
      }
  
      *chptr++ = '>';
    } else if (ch == ':') {
      *chptr++ = ':';
      *chptr++ = ch = getc(fp);
  
      if (ch == EOF){
      if (ferror(fp)) Error("Cannot read input file %s", gpCurrentMmfName);
      else           Error("Unexpected end of file %s", gpCurrentMmfName);
      }
    } else {
      while (ch != EOF && !isspace(ch) && chptr-buffer < sizeof(buffer)-1) {
        *chptr++ = ch;
        ch = getc(fp);
      }
  
      if (ch != EOF) {
        ungetc(ch, fp);
      } else if (ferror(fp)) {
        Error("Cannot read input file %s", gpCurrentMmfName);
      }
  
      if (chptr == buffer) {
        if (eofNotExpected) {
          Error("Unexpected end of file %s", gpCurrentMmfName);
        }
        return NULL;
      }
    }
  
    *chptr = '\0';
  //  puts(buffer);
    return buffer;
  
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void UngetString(void)
  {
    gStringUnget = true;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int GetInt(FILE *fp)
  {
    int   cc;
    short ret;
  //puts("GetInt");
    
    if (gHmmReadBinary) {
      cc = fread(&ret, sizeof(short), 1, fp);
      if (!isBigEndian()) swap2(ret);
    } else {
      RemoveSpaces(fp);
      cc = fscanf(fp, "%hd", &ret);
    }
    
    if (cc != 1) {
      if (ferror(fp)) {
        Error("Cannot read input file %s", gpCurrentMmfName);
      }
      Error("Integral number expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT GetFloat(FILE *fp)
  {
    int cc;
    float ret;
  //puts("GetFloat");
    
    if (gHmmReadBinary) {
      cc = fread(&ret, sizeof(float), 1, fp);
      if (!isBigEndian()) swap4(ret);  
    } else {
      RemoveSpaces(fp);
      cc = fscanf(fp, "%f", &ret);
    }
    
    if (cc != 1) {
      if (ferror(fp)) {
        Error("Cannot read input file %s", gpCurrentMmfName);
      }
      Error("Float number expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void RemoveSpaces(FILE *fp)
  {
    char ch;
  
  //  puts("RemoveSpaces");
    while (isspace(ch = getc(fp))) {
      if (ch == '\n') {
        ++gCurrentMmfLine;
      }
    }
    if (ch != EOF) {
      ungetc(ch, fp);
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  KeywordID ReadOutPDFKind(char *str)
  {
    if (     CheckKwd(str, KID_DiagC))    return KID_DiagC;
    else if (CheckKwd(str, KID_InvDiagC)) return KID_InvDiagC;
    else if (CheckKwd(str, KID_FullC))    return KID_FullC;
    else if (CheckKwd(str, KID_XFormC))   return KID_XFormC;
    else if (CheckKwd(str, KID_LLTC))     return KID_LLTC;
    else if (CheckKwd(str, KID_PDFObsVec))return KID_PDFObsVec;
    else return KID_UNSET;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  KeywordID ReadDurKind(char *str)
  {
    if (     CheckKwd(str, KID_NullD))    return KID_NullD;
    else if (CheckKwd(str, KID_PoissonD)) return KID_PoissonD;
    else if (CheckKwd(str, KID_GammaD))   return KID_GammaD;
    else if (CheckKwd(str, KID_GenD))     return KID_GenD;
    else return KID_UNSET;
  }

  
  //***************************************************************************
  //***************************************************************************
  void PutKwd(FILE *fp, bool binary, KeywordID kwdID)
  {
    if (binary) {
      putc(':', fp);
      putc(kwdID, fp);
    } else {
      putc('<', fp);
      fputs(gpKwds[kwdID], fp);
      fputs("> ", fp);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void PutInt(FILE *fp, bool binary, int i)
  {
    if (binary) {
      short b = i;
      if (!isBigEndian()) swap2(b);
      fwrite(&b, sizeof(short), 1, fp);
    } else {
      fprintf(fp, "%d ", i);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void PutFlt(FILE *fp, bool binary, FLOAT f)
  {
    if (binary) {
      float b = f;
      if (!isBigEndian()) swap4(b);
      fwrite(&b, sizeof(float), 1, fp);
    } else {
      fprintf(fp, FLOAT_FMT" ", f);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void PutNLn(FILE *fp, bool binary)
  {
    if (!binary) putc('\n', fp);
  }
  
  
  //###########################################################################################################
  
  //###########################################################################
  //###########################################################################
  // HMM LIST INPUT
  //###########################################################################
  //###########################################################################
  
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  ReadHMMList(const char * pFileName, 
              const char * pInMmfDir, 
              const char * pInMmfExt)
  {
    struct ReadlineData      rld = {0};
    char *                    lhmm;
    char *                    fhmm;
    char *                    chptr;
    Macro *                   macro;
    Macro *                   macro2;
    int                       line_no = 0;
    char                      mmfile[1024];
    FILE *                    fp;
  
    if ((fp = my_fopen(pFileName, "rt", gpHListFilter)) == NULL) 
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
        
        ParseMmf(mmfile, fhmm);
        
        if ((macro = FindMacro(&mHmmHash, fhmm)) == NULL)
          Error("Definition of model '%s' not found in file '%s'", fhmm, mmfile);
      }
      
      if (lhmm != fhmm) 
      {
        gpCurrentMmfName = NULL; // Global variable; macro will not be written to any output MMF
        macro2 = this->pAddMacro('h', lhmm);
        assert(macro2 != NULL);
        
        if (macro2->mpData != NULL) 
        {
          if (gHmmsIgnoreMacroRedefinition == 0) 
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
  ParseMmf(const char * pFileName, char * expectHMM)
  {
    //ReadHMMSet(rName.c_str(), this, expectHMM);
    istkstream    input_stream;
    FILE *        fp;
    char *        keyword;
    Macro *       macro;
    MacroData *   data;
  
    data = NULL;
    
    gCurrentMmfLine = 1;
    gpCurrentMmfName = pFileName;//mmFileName;
  
    /*
    if ((fp = fopen(mmFileName, "rb")) == NULL) 
    {
      Error("Cannot open input MMF %s", mmFileName);
    }
    */
    
    // try to open the stream
    input_stream.open(pFileName, ios::in|ios::binary);
    
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
        fclose(fp);
        return;
      }
  
      
      if (keyword[0] == '~' && keyword[2] == '\0' ) 
      {
        char type = keyword[1];
  
        if (type == 'o') {
          if (!ReadGlobalOptions(fp)) {
            Error("No global option defined (%s:%d)", pFileName, gCurrentMmfLine);
          }
        } else {
          keyword = GetString(fp, 1);
          if ((macro = pAddMacro(type, keyword)) == NULL) {
            Error("Unrecognized macro type ~%c (%s:%d)", type, pFileName, gCurrentMmfLine);
          }
  
          if (macro->mpData != NULL) {
            if (gHmmsIgnoreMacroRedefinition == 0) {
              Error("Redefinition of macro ~%c %s (%s:%d)", type, keyword, pFileName, gCurrentMmfLine);
            } else {
  //            Warning("Redefinition of macro ~%c %s (%s:%d) is ignored", type, keyword, mmFileName, gCurrentMmfLine);
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
            case 'j': data =  ReadXFormInstance(fp, macro); break;
            case 'x': data =  ReadXForm        (fp, macro); break;
            default : data = NULL; break;
          }  
  
          assert(data != NULL);
          
          if (macro->mpData == NULL) {
            macro->mpData = data;
          } else {
            Warning("Redefinition of macro ~%c %s (%s:%d)",
                    type, keyword, pFileName, gCurrentMmfLine);
  
            // Macro is redefined. New item must be checked for compatibility with
            // the old one (vector sizes, delays, memory sizes) All references to
            // old item must be replaced and old item must be released
            // !!! How about AllocateAccumulatorsForXFormStats() and ResetAccumsForHMMSet()
  
  
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
              //free(ud.mpOldData);
              hash = &mMeanHash;
              break;
            case 'v':
              this->Scan(MTM_MIXTURE, NULL,ReplaceItem, &ud);
              delete ud.mpOldData;
              //free(ud.mpOldData);
              hash = &mVarianceHash;
              break;
            case 't':
              this->Scan(MTM_HMM, NULL,ReplaceItem, &ud);
              delete ud.mpOldData;
              //free(ud.mpOldData);
              hash = &mTransitionHash;
              break;
            case 'j':
              this->Scan(MTM_XFORM_INSTANCE | MTM_MIXTURE, NULL, ReplaceItem, &ud);
              ud.mpOldData->Scan(MTM_REVERSE_PASS|MTM_ALL, NULL, ReleaseItem, NULL);
              hash = &mXFormInstanceHash;
              break;
            case 'x':
              this->Scan(MTM_XFORM | MTM_XFORM_INSTANCE,NULL,ReplaceItem, &ud);
              ud.mpOldData->Scan(MTM_REVERSE_PASS | MTM_ALL, NULL,ReleaseItem,NULL);
              hash = &mXFormHash;
              break;
            }
            
            for (i = 0; i < hash->mNEntries; i++) {
              if (hash->mpEntry[i]->data == ud.mpOldData) {
                hash->mpEntry[i]->data = ud.mpNewData;
              }
            }
          }
        }
      } else if (CheckKwd(keyword, KID_BeginHMM)) {
        UngetString();
        if (expectHMM == NULL) {
          Error("Macro definition expected (%s:%d)",pFileName,gCurrentMmfLine);
        }
        //macro = AddMacroToHMMSet('h', expectHMM, hmm_set);
        macro = pAddMacro('h', expectHMM);
        
        macro->mpData = ReadHMM(fp, macro);
      } else {
        Error("Unexpected keyword %s (%s:%d)",keyword,pFileName,gCurrentMmfLine);
      }
    }
  }
  

  //***************************************************************************
  //***************************************************************************
  Hmm *
  ModelSet::
  ReadHMM(FILE * fp, Macro * macro)
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
        keyword, gpCurrentMmfName, gCurrentMmfLine);
      }
      return (Hmm *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_BeginHMM))
      Error("Keyword <BeginHMM> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
    ReadGlobalOptions(fp);
  
    if (mInputVectorSize == -1)
      Error("<VecSize> is not defined yet (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
    if (mDurKind == -1) 
      mDurKind = KID_NullD;
  
    keyword = GetString(fp, 1);
    
    
    if (!CheckKwd(keyword, KID_NumStates))
      Error("Keyword <NumStates> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
    nstates = GetInt(fp);
    //***
    // old malloc
    // if ((ret = (Hmm *) malloc(sizeof(Hmm) + (nstates-3) * sizeof(State *))) == NULL) 
    //   Error("Insufficient memory");
    // ret->mNStates = nstates;
    ret = new Hmm(nstates);  
    
  
    for (i=0; i<nstates-2; i++) 
      ret->mpState[i] = NULL;
  
    for (i=0; i<nstates-2; i++) 
    {
      keyword = GetString(fp, 1);
      
      if (!CheckKwd(keyword, KID_State)) 
        Error("Keyword <State> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
      state_id = GetInt(fp);
  
  //    printf("%d\n", state_id);
      if (state_id < 2 || state_id >= nstates) 
        Error("State number out of the range (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
      if (ret->mpState[state_id-2] != NULL) 
        Error("Redefinition of state (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
      ret->mpState[state_id-2] = ReadState(fp, NULL);
  //    printf("\n%d: %x\n", state_id-2, ret->mpState[state_id-2]);
    }
  
    ret->mpTransition    = ReadTransition(fp, NULL);
    ret->mpMacro         = macro;
  
    if (ret->mpTransition->mNStates != nstates) 
      Error("Invalid transition matrix size (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
    keyword = GetString(fp, 1);
    
    if (!CheckKwd(keyword, KID_EndHMM))
      Error("Keyword <EndHMM> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    
    return ret;
  }

  
  //***************************************************************************
  //***************************************************************************
  State *
  ModelSet::
  ReadState(FILE * fp, Macro * macro)
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
        Error("Undefined reference to macro ~s %s (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
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
              gpCurrentMmfName, gCurrentMmfLine);
      }
      ret->PDF_obs_coef = GetInt(fp) - 1;
      range = mpInputXForm ? mpInputXForm->mOutSize
                                  : mInputVectorSize;
      if (ret->PDF_obs_coef < 0 || ret->PDF_obs_coef >= range) 
      {
        Error("Parameter <ObsCoef> is out of the range 1:%d (%s:%d)",
              range, gpCurrentMmfName, gCurrentMmfLine);
      }
    } 
    else 
    {
      ret->mNumberOfMixtures = num_mixes;
    //  printf("ptr: %x num_mixes: %d\n", ret, num_mixes);
  
      for (i=0; i<num_mixes; i++) 
        ret->mpMixture[i].mpEstimates = NULL;
  
      if (CheckKwd(keyword, KID_Stream)) 
      {
        if (GetInt(fp) != 1) 
        {
          Error("Stream number out of the range (%s:%d)",
                gpCurrentMmfName, gCurrentMmfLine);
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
                  gpCurrentMmfName, gCurrentMmfLine);
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
                gpCurrentMmfName, gCurrentMmfLine);
        }
  
        if (ret->mpMixture[mixture_id-1].mpEstimates != NULL) 
        {
          Error("Redefinition of mixture %d (%s:%d)",
                mixture_id, gpCurrentMmfName, gCurrentMmfLine);
        }
  
        ret->mpMixture[mixture_id-1].mpEstimates  = ReadMixture(fp, NULL);
        ret->mpMixture[mixture_id-1].mWeight      = log(mixture_weight);
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
  Mixture *
  ModelSet::
  ReadMixture(FILE *fp, Macro *macro)
  {
    Mixture * ret;
    char *    keyword;
    int       size;
  
  //  puts("ReadMixture");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~m")) 
    {
      keyword = GetString(fp, 1);
      
      if ((macro = FindMacro(&mMixtureHash, keyword)) == NULL)
        Error("Undefined reference to macro ~m %s (%s:%d)",
              keyword, gpCurrentMmfName, gCurrentMmfLine);
      
      return (Mixture *) macro->mpData;
    }
  
    //***
    // old free
    //
    //if ((ret = (Mixture *) malloc(sizeof(Mixture))) == NULL)
    //  Error("Insufficient memory");
    ret = new Mixture;  
  
    if (CheckKwd(keyword, KID_InputXForm)) 
    {
      ret->mpInputXForm = ReadXFormInstance(fp, NULL);
    } 
    else 
    {
      ret->mpInputXForm = mpInputXForm;
      UngetString();
    }
  
    ret->mpMean = ReadMean(fp, NULL);
    ret->mpVariance = ReadVariance(fp, NULL);
  
    size = ret->mpInputXForm ? ret->mpInputXForm->mOutSize : mInputVectorSize;
    
    if (size == -1) 
    {
      Error("<VecSize> is not defined yet (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    } 
    else if (ret->mpMean->mVectorSize != size || ret->mpVariance->mVectorSize != size) 
    {
      Error("Invalid mean or variance vector size (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    if ((keyword = GetString(fp, 0)) != NULL && CheckKwd(keyword, KID_GConst)) 
    {
      ret->mGConst = GetFloat(fp);
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
  Mean *
  ModelSet::
  ReadMean(FILE *fp, Macro *macro)
  {
    Mean *    ret;
    char *    keyword;
    int       vec_size;
    int       i;
    int       accum_size = 0;
  
  //  puts("ReadMean");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~u")) 
    {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mMeanHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~u %s (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
      }
      return  (Mean *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_Mean))
      Error("Keyword <Mean> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
    vec_size = GetInt(fp);
  
    //***
    // old malloc
    //  
    //if (mAllocAccums) 
    //  accum_size = (vec_size + 1) * 2; // * 2 for MMI accums
    //
    //ret = (Mean *) malloc(sizeof(Mean) + (vec_size+accum_size-1) * sizeof(FLOAT));
    //if (ret == NULL)  Error("Insufficient memory");
    //ret->mVectorSize = vec_size;
    ret = new Mean(vec_size, mAllocAccums);
  
//  printf("vec_size: %d\n", vec_size);
    for (i=0; i<vec_size; i++) {
      ret->mpVectorO[i] = GetFloat(fp);
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
    int vec_size, i, accum_size = 0;
  
  //  puts("ReadVariance");
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~v")) {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mVarianceHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~v %s (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
  
      }
      return (Variance *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_Variance)) {
      Error("Keyword <Variance> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    vec_size = GetInt(fp);
  
    //***
    // old malloc
    //if (mAllocAccums) accum_size = (2 * vec_size + 1) * 2; // * 2 for MMI accums
    //
    //ret = (Variance *) malloc(sizeof(Variance) + (vec_size+accum_size-1) * sizeof(FLOAT));
    //if (ret == NULL) Error("Insufficient memory");
    //
    //ret->mVectorSize = vec_size;
    
    ret = new Variance(vec_size, mAllocAccums);
  
    for (i=0; i<vec_size; i++) {
      ret->mpVectorO[i] = 1.0 / GetFloat(fp);
    }
  
    ret->mpMacro = macro;
    return ret;
  }
    
  
  //***************************************************************************
  //***************************************************************************
  XFormInstance *
  ModelSet::
  ReadXFormInstance(FILE *fp, Macro *macro) 
  {
    XFormInstance * ret;
    XFormInstance * input = NULL;
    char *          keyword;
    int             out_vec_size = -1;
    int             i;
  
  //  puts("ReadXFormInstance");
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~j")) 
    {
      keyword = GetString(fp, 1);
      
      if ((macro = FindMacro(&mXFormInstanceHash, keyword)) == NULL)
        Error("Undefined reference to macro ~j %s (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
      
      return (XFormInstance *) macro->mpData;
    }
  
    if (CheckKwd(keyword, KID_Input)) 
    {
      input = ReadXFormInstance(fp, NULL);
      keyword = GetString(fp, 1);
    }
  
    if (CheckKwd(keyword, KID_MMFIDMask)) 
    {
      keyword = GetString(fp, 1);
      if (strcmp(keyword, "*"))
        Error("<MMFIdMask> different than '*' is not supported (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
      keyword = GetString(fp, 1);
    }
  
    if ((i = ReadParmKind(keyword, TRUE)) != -1) {
      if (mParamKind != -1 && mParamKind != i) 
        Error("ParamKind mismatch (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
      keyword = GetString(fp, 1);
    }
  
    if (CheckKwd(keyword, KID_LinXForm)) 
    {
      keyword = GetString(fp, 1);
  //    Error("Keyword <LinXForm> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    if (!CheckKwd(keyword, KID_VecSize))
      Error("Keyword <VecSize> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
    out_vec_size = GetInt(fp);
  
    //***
    // old malloc
    // ret = (XFormInstance *) malloc(sizeof(*ret)+(out_vec_size-1)*sizeof(FLOAT));
    // if (ret == NULL) Error("Insufficient memory");
    //   
    // 
    
    // create new instance
    ret = new XFormInstance(out_vec_size);

    // read instance from file    
    ret->mpXForm = ReadXForm(fp, NULL);
    
    
    if (input == NULL && mInputVectorSize == -1) {
      Error("<VecSize> has not been defined yet (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    if (out_vec_size != ret->mpXForm->mOutSize /* * ret->stackSize*/) {
      Error("XFormInstance <VecSize> must equal to XForm "
            "output size (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    }
  
    if (input == NULL) {
      if (ret->mpXForm->mInSize != mInputVectorSize) {
        Error("XForm input size must equal to ~o <VecSize> (%s:%d)",
              gpCurrentMmfName, gCurrentMmfLine);
      }
    } else {
      if (ret->mpXForm->mInSize != input->mOutSize) {
        Error("XForm input size must equal to <Input> <VecSize> (%s:%d)",
              gpCurrentMmfName, gCurrentMmfLine);
      }
    }
  
    ret->mpMemory = NULL;
    
    //***
    // old malloc
    //if (ret->mpXForm->mMemorySize > 0 &&
    //  ((ret->mpMemory = (char *) calloc(1, ret->mpXForm->mMemorySize)) == NULL)) {
    //  Error("Insufficient memory");
    //}
    if (ret->mpXForm->mMemorySize > 0)
      ret->mpMemory = new char[ret->mpXForm->mMemorySize];
  
    ret->mpNext = mpXFormInstances;
    mpXFormInstances = ret;
    ret->mpInput = input;
    ret->mpMacro = macro;
  //  puts("ReadXFormInstance exit");
  
    ret->mNumberOfXFormStatCaches = 0;
    ret->mpXFormStatCache   = NULL;
    ret->mTotalDelay = ret->mpXForm->mDelay + (input ? input->mTotalDelay : 0);
  
    mTotalDelay = HIGHER_OF(mTotalDelay, ret->mTotalDelay);
    return ret;
  }; //ReadXFormInstance(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************
  XForm *
  ModelSet::
  ReadXForm(FILE *fp, Macro *macro) 
  {
    char *keyword;
    unsigned int i;
  
  //  puts("ReadXForm");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~x")) {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mXFormHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~x %s (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
      }
      return (XForm *) macro->mpData;
    }
  
    if (CheckKwd(keyword, KID_XForm)) {
      return (XForm *) ReadLinearXForm(fp, macro);
    }
  
    if (CheckKwd(keyword, KID_Bias)) {
      return (XForm *) ReadBiasXForm(fp, macro);
    }
  
    if (CheckKwd(keyword, KID_Copy)) {
      return (XForm *) ReadCopyXForm(fp, macro);
    }
  
    if (CheckKwd(keyword, KID_Stacking)) {
      return (XForm *) ReadStackingXForm(fp, macro);
    }
  
    if (CheckKwd(keyword, KID_NumLayers) ||
      CheckKwd(keyword, KID_NumBlocks) ||
      CheckKwd(keyword, KID_BlockInfo))
    {
      UngetString();
      return (XForm *) ReadCompositeXForm(fp, macro);
    }
  
    for (i=0; i < gFuncTableSize/sizeof(*gFuncTable); i++) 
    {
      if (CheckKwd(keyword, gFuncTable[i].KID)) {
        return (XForm *) ReadFuncXForm(fp, macro, i);
      }
    }
  
    Error("Invalid XForm definition (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
    return NULL;
  }; //ReadXForm(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************
  //***************************************************************************
  CompositeXForm *
  ModelSet::
  ReadCompositeXForm(FILE *fp, Macro *macro) 
  {
    CompositeXForm *    ret;
    XForm **            block;
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
  
    //if ((ret = (CompositeXForm *) malloc(sizeof(CompositeXForm) +
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
  
    // create new Composite XForm
    ret = new CompositeXForm(nlayers);    
    
    // load each layer
    for (i=0; i<nlayers; i++) 
    {
      keyword = GetString(fp, 1);
      
      if (!CheckKwd(keyword, KID_Layer)) 
      {
        if (nlayers > 1) 
          Error("Keyword <Layer> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
        
        layer_id = 1;
      } 
      else 
      {
        layer_id = GetInt(fp);
        keyword = GetString(fp, 1);
      }
  
      if (layer_id < 1 || layer_id > nlayers)
        Error("Layer number out of the range (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
      
      if (ret->mpLayer[layer_id-1].mpBlock != NULL)
        Error("Redefinition of mpLayer (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
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
      //       if ((block = (XForm **)  malloc(sizeof(XForm*) * nblocks)) == NULL) 
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
          if (nblocks > 1) 
            Error("Keyword <Block> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
          
          UngetString();
          block_id = 1;
        } 
        else 
        {
          block_id = GetInt(fp);
        }
  
        if (block_id < 1 || block_id > nblocks)
          Error("Block number out of the range (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
        if (block[block_id-1] != NULL)
          Error("Redefinition of block (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
  
        block[block_id-1] = ReadXForm(fp, NULL);
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
              i, prev_out_size, i+1, layer_in_size, gpCurrentMmfName, gCurrentMmfLine);
        }
  
        //***
        // old malloc
        //
        //if ((ret->mpLayer[i-1].mpOutputVector = (FLOAT *) malloc(prev_out_size * sizeof(FLOAT))) == NULL) 
        //{
        //  Error("Insufficient memory");
        //}
        
        ret->mpLayer[i-1].mpOutputVector = new FLOAT[prev_out_size];
      }
  
      prev_out_size = layer_out_size;
    }
  
    
    ret->mpMacro = macro;
    return ret;
  }; // ReadCompositeXForm(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************    
  //***************************************************************************
  LinearXForm *
  ModelSet::
  ReadLinearXForm(FILE *fp, Macro *macro)
  {
    LinearXForm *   ret;
    int             in_size;
    int             out_size;
    int             i;
    int             r;
    int             c;
    FLOAT           tmp;
  
    // read the parameters
    out_size = GetInt(fp);
    in_size = GetInt(fp);
  
    //*** 
    // old malloc
    //ret = (LinearXForm *) malloc(sizeof(LinearXForm)+(out_size*in_size-1)*sizeof(ret->mpMatrixO[0]));
    //if (ret == NULL) Error("Insufficient memory");
    
    // create new object
    ret = new LinearXForm(in_size, out_size);
    
    // fill the object with data    
    i = 0;
    for (r=0; r < in_size; r++)
    {
      for (c = 0; c < out_size; c++)
      {
        tmp = GetFloat(fp);
        ret->mMatrix(r, c) = tmp;
        //:TODO: Get rid of this old stuff
        ret->mpMatrixO[i] = tmp;
        i++;
      }
    }
  
    ret->mpMacro      = macro;
    return ret;
  }; // ReadLinearXForm(FILE *fp, Macro *macro)
  

  //***************************************************************************  
  BiasXForm *
  ModelSet::
  ReadBiasXForm(FILE *fp, Macro *macro)
  {
    BiasXForm *   ret;
    size_t        size;
    size_t        i;
    FLOAT         tmp;
    
    size = GetInt(fp);
  
    //***
    // old malloc
    //     ret = (BiasXForm *) malloc(sizeof(BiasXForm)+(size-1)*sizeof(ret->mpVectorO[0]));
    //     if (ret == NULL) Error("Insufficient memory");
    //     
    
    ret = new BiasXForm(size);
    
    // load values
    for (i=0; i < size; i++) 
    {
      tmp = GetFloat(fp);
      ret->mVector(0, i);
      
      //:TODO: Get rid of the obsolete thing
      ret->mpVectorO[i] = tmp;      
    }
  
    ret->mpMacro      = macro;
    return ret;
  }; // ReadBiasXForm(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  //***************************************************************************
  FuncXForm *
  ModelSet::
  ReadFuncXForm(FILE *fp, Macro *macro, int funcId)
  {
    FuncXForm * ret;
    size_t      size;
  
    //***
    // old malloc
    //ret = (FuncXForm *) malloc(sizeof(FuncXForm));
    //if (ret == NULL) Error("Insufficient memory");
    
    
    size          = GetInt(fp);
    ret           = new FuncXForm(size, funcId);
    ret->mpMacro  = macro;
    return ret;
  }; // ReadFuncXForm(FILE *fp, Macro *macro, int funcId)
  
  //***************************************************************************
  //***************************************************************************
  CopyXForm *
  ModelSet::
  ReadCopyXForm(FILE *fp, Macro *macro)
  {
    CopyXForm *   ret;
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
    //ret = (CopyXForm *) malloc(sizeof(CopyXForm)+(out_size-1)*sizeof(int));
    //if (ret == NULL) Error("Insufficient memory");
    
    // create new object
    ret = new CopyXForm(in_size, out_size);
  
    // fill it with data
    while (i < out_size) 
    {
      RemoveSpaces(fp);
      if ((n = fscanf(fp, "%d:%d:%d", &from, &step, &to)) < 1) {
        if (ferror(fp)) {
          Error("Cannot read input file %s", gpCurrentMmfName);
        }
        Error("Integral number expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
      }
  
      if (n == 2)      { to = step; step = 1; }
      else if (n == 1) { to = from; step = 1; }
  
      if (to < 1 || to > in_size) {
        Error("Copy index %d out of range (%s:%d)",
              to, gpCurrentMmfName, gCurrentMmfLine);
      }
  
      for (n = 0; n < (to-from)/step + 1; n++, i++) {
        ret->mpIndices[i] = from + n * step - 1;
      }
    }
  
    ret->mpMacro      = macro;
    return ret;
  }; //ReadCopyXForm(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  StackingXForm *
  ModelSet::
  ReadStackingXForm(FILE *fp, Macro *macro)
  {
    StackingXForm *   ret;
    size_t            stack_size;
    size_t            in_size;
  
    stack_size = GetInt(fp);
    in_size    = GetInt(fp);
  
    //***
    // old malloc
    //ret = (StackingXForm *) malloc(sizeof(StackingXForm));
    //if (ret == NULL) Error("Insufficient memory");
    ret = new StackingXForm(stack_size, in_size);
  
    ret->mpMacro      = macro;
    return ret;
  }; //ReadStackingXForm(FILE *fp, Macro *macro)
  
  //***************************************************************************
  int 
  ModelSet::
  ReadGlobalOptions(FILE *fp)
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
  //      Error("Unexpected end of file %s", gpCurrentMmfName);
      }
  
      if (CheckKwd(keyword, KID_VecSize)) 
      {
        i = GetInt(fp);
        if (mInputVectorSize != -1 && mInputVectorSize != i) 
        {
          Error("Mismatch in <VecSize> redefinition (%s:%d)",
                gpCurrentMmfName, gCurrentMmfLine);
        }
  
        mInputVectorSize = i;
        ret = 1;
      } 
      else if (CheckKwd(keyword, KID_StreamInfo)) 
      {
        if (GetInt(fp) != 1) 
        {
          Error("Unsupported definition of multistream (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
        }
  
        i = GetInt(fp);
        if (mInputVectorSize != -1 && mInputVectorSize != i) 
        {
          Error("Mismatch in <VecSize> redefinition (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
        }
  
        mInputVectorSize = i;
        ret = 1;
      } 
      else if ((i = ReadParmKind(keyword, TRUE)) != -1) 
      {
        if (mParamKind != -1 && mParamKind != i) 
          Error("Mismatch in paramKind redefinition (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
        
        mParamKind = i;
        ret = 1;
      } 
      else if ((i = ReadOutPDFKind(keyword)) != -1) 
      {
        if (mOutPdfKind != -1 && mOutPdfKind != i) {
          Error("Mismatch in outPDFKind redefinition (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
        }
  
        if (i != KID_PDFObsVec && i != KID_DiagC) {
          Error("Unsupported option '%s' (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
        }
  
        mOutPdfKind = static_cast<KeywordID> (i);
        ret = 1;
      } 
      else if ((i = ReadDurKind(keyword)) != -1) 
      {
        if (mDurKind != -1 && mDurKind != i) {
          Error("Mismatch in durKind redefinition (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
        }
        if (i != KID_NullD) {
          Error("Unsupported option '%s' (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
        }
        mDurKind = static_cast<KeywordID> (i);
        ret = 1;
      } else if (CheckKwd(keyword, KID_HMMSetID)) {
        ret = 1;
      } else if (CheckKwd(keyword, KID_InputXForm)) {
        XFormInstance *inputXForm;
        Macro *macro = pAddMacro(mt_XFormInstance, DEFAULT_XFORM_NAME);
        
        if (macro->mpData != NULL) 
        {
          if (gHmmsIgnoreMacroRedefinition == 0)
            Error("Redefinition of <InputXForm> (%s:%d)",
                  gpCurrentMmfName, gCurrentMmfLine);
        }
        
        inputXForm = ReadXFormInstance(fp, macro);
  
        if (macro->mpData != NULL) {
          Warning("Redefinition of <InputXForm> (%s:%d)",
                  gpCurrentMmfName, gCurrentMmfLine);
  
          // Macro is redefined. New item must be checked for compatibility with
          // the old one (vector size) All references to old
          // item must be replaced and old item must be released
  
          ReplaceItemUserData ud;
          ud.mpOldData = macro->mpData;
          ud.mpNewData = inputXForm;
          ud.mType     = 'j';
  
          this         ->Scan(MTM_XFORM_INSTANCE|MTM_MIXTURE,NULL,ReplaceItem, &ud);
          ud.mpOldData ->Scan(MTM_REVERSE_PASS|MTM_ALL, NULL, ReleaseItem, NULL);
  
          for (unsigned int i = 0; i < mXFormInstanceHash.mNEntries; i++) 
          {
            if (mXFormInstanceHash.mpEntry[i]->data == ud.mpOldData) 
            {
              mXFormInstanceHash.mpEntry[i]->data = ud.mpNewData;
            }
          }
        } 
        else 
        {
          mpInputXForm = inputXForm;
          macro->mpData = inputXForm;
        }
        
        ret = 1;
  //    } else if (CheckKwd(keyword, KID_LinXForm)) {
  //      linXForm = ReadXFormInstance(fp, hmm_set, NULL);
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
        Error("Undefined reference to macro ~t %s (%s:%d)", keyword, gpCurrentMmfName, gCurrentMmfLine);
      }
      return (Transition *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_TransP)) {
      Error("Keyword <TransP> expected (%s:%d)", gpCurrentMmfName, gCurrentMmfLine);
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
      ret->mpMatrixO[i] = (float) ret->mpMatrixO[i] != 0.0 ? log(ret->mpMatrixO[i]) : LOG_0;
    }
  
    ret->mpMacro = macro;
    return ret;
  }
  

    
  //###########################################################################################################
  
  //###########################################################################
  //###########################################################################
  // MMF OUTPUT
  //###########################################################################
  //###########################################################################
    
  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  WriteMmf(const char * pFileName, const char * pOutputDir,
           const char * pOutputExt, bool binary)
  {
    FILE *    fp = NULL;
    Macro *   macro;
    char      mmfile[1024];
    char *    lastFileName = NULL;
    int       waitingForNonXForm = 1;
  
    for (macro = mpFirstMacro; macro != NULL; macro = macro->nextAll) 
    {
      if (macro->mpFileName == NULL) 
        continue; // Artificial macro not read from file
  
      if (lastFileName == NULL || (!pFileName && strcmp(lastFileName, macro->mpFileName))) 
      {
        // New macro file
        lastFileName = macro->mpFileName;
        
        if (fp && fp != stdout) 
          fclose(fp);
  
        if (!strcmp(pFileName ? pFileName : macro->mpFileName, "-")) 
        {
          fp = stdout;
        }
        else 
        {
          MakeFileName(mmfile, pFileName ? pFileName : macro->mpFileName, pOutputDir, pOutputExt);
          
          if ((fp  = fopen(mmfile, "wb")) == NULL) 
            Error("Cannot open output MMF %s", mmfile);
        }
        
        waitingForNonXForm = 1;
        WriteGlobalOptions(fp, binary);
      }
  
      if (macro->mpData == mpInputXForm &&
        (!strcmp(macro->mpName, DEFAULT_XFORM_NAME))) 
      {
        fputs("~o ", fp);
        PutKwd(fp, binary, KID_InputXForm);
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
          case 'x': WriteXForm        (fp, binary, static_cast <XForm*>         (macro->mpData)); break;
          case 'j': WriteXFormInstance(fp, binary, static_cast <XFormInstance*> (macro->mpData)); break;
          case 'u': WriteMean         (fp, binary, static_cast <Mean*>          (macro->mpData)); break;
          case 'v': WriteVariance     (fp, binary, static_cast <Variance*>      (macro->mpData)); break;
          case 't': WriteTransition   (fp, binary, static_cast <Transition*>    (macro->mpData)); break;
          case 'm': WriteMixture      (fp, binary, static_cast <Mixture*>       (macro->mpData)); break;
          case 's': WriteState        (fp, binary, static_cast <State*>         (macro->mpData)); break;
          case 'h': WriteHMM          (fp, binary, static_cast <Hmm*>           (macro->mpData)); break;
        }
      }
    }
    if (fp && fp != stdout) fclose(fp);
  }


  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  WriteGlobalOptions(FILE *fp, bool binary)
  {
    char parmkindstr[64];
  
    fputs("~o ", fp);
    PutKwd(fp, binary, KID_VecSize);
    PutInt(fp, binary, mInputVectorSize);
  
    if (ParmKind2Str(mParamKind, parmkindstr)) {
      fprintf(fp, "<%s> ", parmkindstr);
    }
    if (mOutPdfKind != -1) {
      PutKwd(fp, binary, mOutPdfKind);
    }
    if (mOutPdfKind != -1) {
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
  

  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteState(FILE *fp, bool binary, State *state)
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
      if (state->mNumberOfMixtures > 1) 
      {
        PutKwd(fp, binary, KID_NumMixes);
        PutInt(fp, binary, state->mNumberOfMixtures);
        PutNLn(fp, binary);
      }
  
      for (i=0; i < state->mNumberOfMixtures; i++) 
      {
        if (state->mNumberOfMixtures > 1) 
        {
          PutKwd(fp, binary, KID_Mixture);
          PutInt(fp, binary, i+1);
          PutFlt(fp, binary, exp(state->mpMixture[i].mWeight));
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
  WriteMixture(FILE *fp, bool binary, Mixture *mixture)
  {
    if (mixture->mpInputXForm != mpInputXForm) 
    {
      PutKwd(fp, binary, KID_InputXForm);
      if (mixture->mpInputXForm->mpMacro) 
      {
        fprintf(fp, "~j \"%s\"", mixture->mpInputXForm->mpMacro->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        WriteXFormInstance(fp, binary, mixture->mpInputXForm);
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
    PutFlt(fp, binary, mixture->mGConst);
    PutNLn(fp, binary);
  }


  //***************************************************************************  
  //*****************************************************************************  
  void
  ModelSet::
  WriteMean(FILE *fp, bool binary, Mean *mean)
  {
    size_t    i;
  
    PutKwd(fp, binary, KID_Mean);
    PutInt(fp, binary, mean->mVectorSize);
    PutNLn(fp, binary);
  
    for (i=0; i < mean->mVectorSize; i++) {
      PutFlt(fp, binary, mean->mpVectorO[i]);
    }
  
    PutNLn(fp, binary);
  } //WriteMean(FILE *fp, bool binary, Mean *mean)

  
  //***************************************************************************  
  //*****************************************************************************  
  void
  ModelSet::
  WriteVariance(FILE *fp, bool binary, Variance *variance)
  {
    size_t   i;
  
    PutKwd(fp, binary, KID_Variance);
    PutInt(fp, binary, variance->mVectorSize);
    PutNLn(fp, binary);
  
    for (i=0; i < variance->mVectorSize; i++) {
      PutFlt(fp, binary, 1/variance->mpVectorO[i]);
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
        PutFlt(fp, binary, logtp > LOG_MIN ? exp(logtp) : 0.0);
      }
  
      PutNLn(fp, binary);
    }
  }

  
  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteXFormInstance(FILE *fp, bool binary, XFormInstance *xformInstance)
  {
    size_t            i;
    bool              isHTKCompatible = true;
    char              parmkindstr[64];
  
    CompositeXForm *cxf = (CompositeXForm *) xformInstance->mpXForm;
    
    if (xformInstance->mpInput != NULL || cxf == NULL || cxf->mpMacro ||
      cxf->mXFormType != XT_COMPOSITE || cxf->mNLayers != 1) 
    {
      isHTKCompatible = 0;
    } 
    else
    {
      for (i = 0; i < cxf->mpLayer[0].mNBlocks; i++) 
      {
        if (cxf->mpLayer[0].mpBlock[i]->mXFormType != XT_LINEAR) 
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
        WriteXFormInstance(fp, binary, xformInstance->mpInput);
      }
    }
  
    if (isHTKCompatible) {
      PutKwd(fp, binary, KID_MMFIDMask);
      fputs("* ", fp);
  
      if (ParmKind2Str(mParamKind, parmkindstr)) 
      {
        fprintf(fp, "<%s> ", parmkindstr);
      }
  
      PutKwd(fp, binary, KID_LinXForm);
    }
  
    PutKwd(fp, binary, KID_VecSize);
    PutInt(fp, binary, xformInstance->mOutSize);
    PutNLn(fp, binary);
  
    if (xformInstance->mpXForm->mpMacro) 
    {
      fprintf(fp, "~x \"%s\"", xformInstance->mpXForm->mpMacro->mpName);
      PutNLn(fp, binary);
    } 
    else 
    {
      WriteXForm(fp, binary, xformInstance->mpXForm);
    }
  }

  
  //***************************************************************************
  //*****************************************************************************  
  void 
  ModelSet::
  WriteXForm(FILE *fp, bool binary, XForm *xform)
  {
    XFormType type = xform->mXFormType;
    
    switch (type)
    {
      case XT_LINEAR    : WriteLinearXForm   (fp, binary, static_cast<LinearXForm *>(xform)); break;
      case XT_COPY      : WriteCopyXForm     (fp, binary, static_cast<CopyXForm *>(xform)); break;
      case XT_FUNC      : WriteFuncXForm     (fp, binary, static_cast<FuncXForm *>(xform)); break;
      case XT_BIAS      : WriteBiasXForm     (fp, binary, static_cast<BiasXForm *>(xform)); break;
      case XT_STACKING  : WriteStackingXForm (fp, binary, static_cast<StackingXForm *>(xform)); break;
      case XT_COMPOSITE : WriteCompositeXForm(fp, binary, static_cast<CompositeXForm *>(xform)); break;
      default:  break;    
    }
  } //WriteXForm(FILE *fp, bool binary, XForm *xform)

  
  //***************************************************************************
  //*****************************************************************************  
  void 
  ModelSet::
  WriteCompositeXForm(FILE *fp, bool binary, CompositeXForm *xform)
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
        if (xform->mpLayer[0].mpBlock[i]->mXFormType != XT_LINEAR) 
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
          WriteXForm(fp, binary, xform->mpLayer[i].mpBlock[j]);
        }
      }
    }
  } //WriteCompositeXForm(FILE *fp, bool binary, CompositeXForm *xform)

  
  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteFuncXForm(FILE *fp, bool binary, FuncXForm *xform)
  {
    PutKwd(fp, binary, gFuncTable[xform->mFuncId].KID);
    PutInt(fp, binary, xform->mOutSize);
    PutNLn(fp, binary);
  } //WriteFuncXForm(FILE *fp, bool binary, FuncXForm *xform)


  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteBiasXForm(FILE *fp, bool binary, BiasXForm *xform)
  {
    size_t  i;
  
    PutKwd(fp, binary, KID_Bias);
    PutInt(fp, binary, xform->mOutSize);
    PutNLn(fp, binary);
    for (i=0; i < xform->mOutSize; i++) 
    {
      PutFlt(fp, binary, xform->mpVectorO[i]);
    }
    PutNLn(fp, binary);
  } //WriteBiasXForm(FILE *fp, bool binary, BiasXForm *xform)


  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteLinearXForm(FILE *fp, bool binary, LinearXForm *xform)
  {
    size_t  i;
    size_t  j;
  
    PutKwd(fp, binary, KID_XForm);
    PutInt(fp, binary, xform->mOutSize);
    PutInt(fp, binary, xform->mInSize);
    PutNLn(fp, binary);
    for (i=0; i < xform->mOutSize; i++) 
    {
      for (j=0; j < xform->mInSize; j++) 
      {
        PutFlt(fp, binary, xform->mpMatrixO[i * xform->mInSize + j]);
      }
      PutNLn(fp, binary);
    }
  } //WriteLinearXForm(FILE *fp, bool binary, LinearXForm *xform)
  

  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteStackingXForm(FILE *fp, bool binary, StackingXForm *xform)
  {
    PutKwd(fp, binary, KID_Stacking);
    PutInt(fp, binary, xform->mOutSize / xform->mInSize);
    PutInt(fp, binary, xform->mInSize);
    PutNLn(fp, binary);
  } //WriteStackingXForm(FILE *fp, bool binary, StackingXForm *xform)
  
  
  //*****************************************************************************
  //*****************************************************************************  
  void
  ModelSet::
  WriteCopyXForm(FILE *fp, bool binary, CopyXForm *xform)
  {
    size_t      i;
    size_t      j;
    int         step = 0;
    int *       ids = xform->mpIndices;
  
    fprintf(fp, "<Copy> %d %d\n", xform->mOutSize, xform->mInSize);
  
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
  } //WriteCopyXForm(FILE *fp, bool binary, CopyXForm *xform)

  

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
  ReadXFormList(const char * pFileName)
  {
    char      line[1024];
    FILE *    fp;
    size_t    nlines=0;
  
    if ((fp = fopen(pFileName, "rt")) == NULL)
        Error("ReadXFormList: Cannot open file: '%s'", pFileName);
  
    while (fgets(line, sizeof(line), fp)) 
    {
      char *            xformName;
      char *            makeXFormShellCommand = NULL;
      char              termChar = '\0';
      size_t            i = strlen(line);
      Macro *           macro;
      MakeXFormCommand *mxfc;
  
      nlines++;
      
      if (line[i-1] != '\n' && getc(fp) != EOF) 
      {
        Error("ReadXFormList: Line %d is too long in file: %s",
              nlines, pFileName);
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
          Error("ReadXFormList: Terminanting %c expected at line %d in file %s",
                termChar, nlines, pFileName);
        }
        
        line[i++] = '\0';
      }
  
      if (line[i] != '\0') 
      { // shell command follows
        for (; isspace(line[i]); i++) 
          line[i] = '\0';
        makeXFormShellCommand = &line[i];
      }
  
      macro = FindMacro(&mXFormHash, xformName);
      
      if (macro == NULL) 
      {
        Error("ReadXFormList: Undefined XForm '%s' at line %d in file %s",
              xformName, nlines, pFileName);
      }
  
      mpXFormToUpdate = (MakeXFormCommand*)
        realloc(mpXFormToUpdate,
                sizeof(MakeXFormCommand) * ++mNumberOfXFormsToUpdate);
  
      if (mpXFormToUpdate == NULL) {
        Error("ReadXFormList: Insufficient memory");
      }
  
      mxfc = &mpXFormToUpdate[mNumberOfXFormsToUpdate-1];
      mxfc->mpXForm = (XForm *) macro->mpData;
      mxfc->mpShellCommand = NULL;
  
      if (makeXFormShellCommand) {
        if ((mxfc->mpShellCommand = strdup(makeXFormShellCommand)) == NULL) {
          Error("ReadXFormList: Insufficient memory");
        }
      }
    }
    fclose(fp);
  }; //ReadXFormList(const std::string & rFileName);
  
}; // namespace STK
