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

#include "Models.h"
#include "stkstream.h"
#include "common.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

namespace STK
{
  int             hmm_read_binary;
  int             current_mmf_line = 1;
  int             string_unget = 0;
  const char *    current_mmf_name;
  const char *    gpHListFilter;
  bool            gHmmsIgnoreMacroRedefinition = true;
  char *          gpKwds[KID_MaxKwdID] = {0};
  
  FLOAT weight_accum_den;
  
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

  
  FunctionTable gFuncTable[] = 
  {
    {sigmoid_vec, KID_Sigmoid}, {log_vec,     KID_Log},
    {exp_vec,     KID_Exp},     {sqrt_vec,    KID_Sqrt},
    {softmax_vec, KID_SoftMax}
  };

  
  void ReplaceItem(int macro_type, HMMSetNodeName nodeName, MacroData * pData, void * pUserData)
  {
    ReplaceItemUserData *  ud = (ReplaceItemUserData *) pUserData;
    size_t                  i;
    size_t                  j;
  
    if (macro_type == 'h') 
    {
      Hmm *hmm = (Hmm *) pData;
      
      if (ud->mType == 's') 
      {
        for (i = 0; i < hmm->mNStates-2; i++) 
        {
          if (hmm->state[i] == ud->mpOldData) 
          {
            hmm->state[i] = (State*) ud->mpNewData;
          }
        }
      } 
      else if (hmm->mpTransition == ud->mpOldData) 
      {
        hmm->mpTransition = (Transition*) ud->mpNewData;
      }
    } 
    
    else if (macro_type == 's') 
    {
      State *state = (State *) pData;
      if (state->mOutPdfKind != KID_PDFObsVec) 
      {
        for (i = 0; i < state->mNumberOfMixtures; i++) 
        {
          if (state->mpMixture[i].estimates == ud->mpOldData) 
          {
            state->mpMixture[i].estimates = (Mixture*) ud->mpNewData;
          }
        }
      }
    } 
    
    else if (macro_type == 'm') 
    {
      Mixture *mixture = (Mixture *) pData;
      if (mixture->mpMean      == ud->mpOldData) mixture->mpMean       = (Mean*)          ud->mpNewData;
      if (mixture->mpVariance  == ud->mpOldData) mixture->mpVariance   = (Variance*)      ud->mpNewData;
      if (mixture->mpInputXForm== ud->mpOldData) mixture->mpInputXForm = (XFormInstance*) ud->mpNewData;
    } 
    
    else if (macro_type == 'x') 
    {
      CompositeXForm *cxf = (CompositeXForm *) pData;
      if (cxf->mXFormType == XT_COMPOSITE) 
      {
        for (i = 0; i < cxf->mNLayers; i++) 
        {
          for (j = 0; j < cxf->layer[i].mNBlocks; j++) 
          {
            if (cxf->layer[i].block[j] == ud->mpOldData) 
            {
              cxf->layer[i].block[j] = (XForm*) ud->mpNewData;
            }
          }
        }
      }
    } 
    
    else if (macro_type == 'j') 
    {
      XFormInstance *xformInstance = (XFormInstance *) pData;
      if (xformInstance->mpInput == ud->mpOldData) xformInstance->mpInput = (XFormInstance*) ud->mpNewData;
      if (xformInstance->mpXForm == ud->mpOldData) xformInstance->mpXForm = (XForm*)         ud->mpNewData;
    }
  }
  
    
  bool IsXFormIn1stLayer(XForm *xform, XForm *topXForm)
  {
    size_t      i;
    if (topXForm == NULL)  return false;
    if (topXForm == xform) return true;
  
    if (topXForm->mXFormType == XT_COMPOSITE) 
    {
      CompositeXForm *cxf = static_cast<CompositeXForm *>(topXForm);
      
      for (i=0; i < cxf->layer[0].mNBlocks; i++) 
      {
        if (IsXFormIn1stLayer(xform, cxf->layer[0].block[i])) 
          return true;
      }
    }
    
    return false;
  }
  
  
  bool Is1Layer1BlockLinearXForm(XForm * pXForm)
  {
    CompositeXForm *cxf = static_cast<CompositeXForm *>(pXForm);
    
    if (cxf == NULL)                                        return false;
    if (cxf->mXFormType == XT_LINEAR)                       return true;
    if (cxf->mXFormType != XT_COMPOSITE)                    return false;
    if (cxf->mNLayers > 1 || cxf->layer[0].mNBlocks > 1)    return false;
    
    return Is1Layer1BlockLinearXForm(cxf->layer[0].block[0]);
  }
  
  
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
  
  KeywordID ReadDurKind(char *str)
  {
    if (     CheckKwd(str, KID_NullD))    return KID_NullD;
    else if (CheckKwd(str, KID_PoissonD)) return KID_PoissonD;
    else if (CheckKwd(str, KID_GammaD))   return KID_GammaD;
    else if (CheckKwd(str, KID_GenD))     return KID_GenD;
    else return KID_UNSET;
  }
  
  
  
  Macro *FindMacro(struct my_hsearch_data *macro_hash, const char *name) 
  {
    ENTRY e, *ep;
    e.key = (char *) name;
    my_hsearch_r(e, FIND, &ep, macro_hash);
    return (Macro *) (ep ? ep->data : NULL);
  }
  
  void ReleaseMacroHash(struct my_hsearch_data *macro_hash) 
  {
    unsigned int i;
    for (i = 0; i < macro_hash->nentries; i++) 
    {
      Macro *macro = (Macro *) macro_hash->entry[i]->data;
      free(macro->mpName);
      free(macro->mpFileName);
      free(macro);
      macro_hash->entry[i]->data = NULL;
    }
    my_hdestroy_r(macro_hash, 0);
  }
  
  void ReleaseItem(int macro_type,HMMSetNodeName nodeName,MacroData * pData, void * pUserData)
  {
    unsigned int i;
    
    CompositeXForm *cxf = (CompositeXForm *) pData;
    
    if (macro_type == 'x' && cxf->mXFormType == XT_COMPOSITE) 
    {
      for (i = 0; i < cxf->mNLayers-1; i++) free(cxf->layer[i].out_vec);
      for (i = 0; i < cxf->mNLayers;   i++) free(cxf->layer[i].block);
    }
    if (macro_type == 'j') free(((XFormInstance *) pData)->memory);
  
    free(pData);
  }
  
  
  
  int CheckKwd(const char *str, KeywordID kwdID)
  {
    const char *chptr;
    if (str[0] == ':') {
      hmm_read_binary = 1;
      return str[1] == kwdID;
    }
  
    if (str[0] != '<') return 0;
    for (chptr = gpKwds[kwdID], str++; *chptr; chptr++, str++) {
      if (toupper(*chptr) != toupper(*str)) return 0;
    }
  
    if (str[0] != '>') return 0;
  
    assert(str[1] == '\0');
    hmm_read_binary = 0;
    return 1;
  }
  
  char *GetString(FILE *fp, int eofNotExpected)
  {
    static char buffer[1024];
    char ch, *chptr = buffer;
    int lines = 0;
  
  //  fputs("GetString: ", stdout);
    if (string_unget) {
      string_unget = 0;
  
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
        Error("Cannot read input file %s", current_mmf_name);
      }
  
      if (ch != termChar) {
        Error("Unterminated string constant (%s:%d)", current_mmf_name, current_mmf_line);
      }
      current_mmf_line += lines;
    } else if (ch == '<') {
      *chptr++ = '<';
      while ((ch = getc(fp)) != EOF && !isspace(ch) && ch != '>' && chptr-buffer < sizeof(buffer)-1) {
        *chptr++ = ch;
      }
  
      if (ch == EOF && ferror(fp)) {
        Error("Cannot read input file %s", current_mmf_name);
      }
  
      if (ch != '>') {
        Error("Unterminated keyword %s (%s:%d)", buffer, current_mmf_name, current_mmf_line);
      }
  
      *chptr++ = '>';
    } else if (ch == ':') {
      *chptr++ = ':';
      *chptr++ = ch = getc(fp);
  
      if (ch == EOF){
      if (ferror(fp)) Error("Cannot read input file %s", current_mmf_name);
      else           Error("Unexpected end of file %s", current_mmf_name);
      }
    } else {
      while (ch != EOF && !isspace(ch) && chptr-buffer < sizeof(buffer)-1) {
        *chptr++ = ch;
        ch = getc(fp);
      }
  
      if (ch != EOF) {
        ungetc(ch, fp);
      } else if (ferror(fp)) {
        Error("Cannot read input file %s", current_mmf_name);
      }
  
      if (chptr == buffer) {
        if (eofNotExpected) {
          Error("Unexpected end of file %s", current_mmf_name);
        }
        return NULL;
      }
    }
  
    *chptr = '\0';
  //  puts(buffer);
    return buffer;
  
  }
  
  void UngetString(void)
  {
    string_unget = 1;
  }
  
  int GetInt(FILE *fp)
  {
    int   cc;
    short ret;
  //puts("GetInt");
    
    if (hmm_read_binary) {
      cc = fread(&ret, sizeof(short), 1, fp);
      if (!isBigEndian()) swap2(ret);
    } else {
      RemoveSpaces(fp);
      cc = fscanf(fp, "%hd", &ret);
    }
    
    if (cc != 1) {
      if (ferror(fp)) {
        Error("Cannot read input file %s", current_mmf_name);
      }
      Error("Integral number expected (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    return ret;
  }
  
  FLOAT GetFloat(FILE *fp)
  {
    int cc;
    float ret;
  //puts("GetFloat");
    
    if (hmm_read_binary) {
      cc = fread(&ret, sizeof(float), 1, fp);
      if (!isBigEndian()) swap4(ret);  
    } else {
      RemoveSpaces(fp);
      cc = fscanf(fp, "%f", &ret);
    }
    
    if (cc != 1) {
      if (ferror(fp)) {
        Error("Cannot read input file %s", current_mmf_name);
      }
      Error("Float number expected (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    return ret;
  }
  
  void RemoveSpaces(FILE *fp)
  {
    char ch;
  
  //  puts("RemoveSpaces");
    while (isspace(ch = getc(fp))) {
      if (ch == '\n') {
        ++current_mmf_line;
      }
    }
    if (ch != EOF) {
      ungetc(ch, fp);
    }
  }
  
  int qsmacrocmp(const void *a, const void *b) 
  {
    return strcmp(((Macro *)a)->mpName, ((Macro *)b)->mpName);
  }
  
  
  
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
  
  void PutNLn(FILE *fp, bool binary)
  {
    if (!binary) putc('\n', fp);
  }
  
  
  
  
  //*****************************************************************************
  
  void ResetAccum(int macro_type, HMMSetNodeName nodeName,
                  MacroData * pData, void *pUserData) 
  {
    size_t    i;
    size_t    size;
    FLOAT *   vector = NULL;
  
    if (macro_type == mt_mean || macro_type == mt_variance) {
      if (macro_type == mt_mean) {
        size   = ((Mean *)pData)->mVectorSize;
        vector = ((Mean *)pData)->mVector + size;
        size   = (size + 1) * 2;
      } else if (macro_type == mt_variance) {
        size   = ((Variance *)pData)->mVectorSize;
        vector = ((Variance *)pData)->mVector + size;
        size   = (size * 2 + 1) * 2;
      }
  
      for (i = 0; i < size; i++) vector[i] = 0;
  
    } else if (macro_type == mt_state) {
      State *state = (State *) pData;
      if (state->mOutPdfKind == KID_DiagC) {
        for (i = 0; i < state->mNumberOfMixtures; i++) {
          state->mpMixture[i].weight_accum     = 0;
          state->mpMixture[i].weight_accum_den = 0;
        }
      }
    } else if (macro_type == mt_transition) {
      size   = SQR(((Transition *) pData)->mNStates);
      vector = ((Transition *) pData)->matrix + size;
  
      for (i = 0; i < size; i++) vector[i] = LOG_0;
    }
  }
  
  
  
  
  
  
  
  
  
  
  void GlobalStats(int macro_type, HMMSetNodeName nn, MacroData * pData, void * pUserData)
  {
    size_t        i;
    
    GlobalStatsUserData *ud = (GlobalStatsUserData *) pUserData;
  
    if (macro_type == mt_state) 
    {
      State * state = (State *) pData;
      
      if (state->mOutPdfKind == KID_DiagC) 
      {
        for (i = 0; i < state->mNumberOfMixtures; i++) 
        {
          state->mpMixture[i].weight_accum += 1;
        }
      }
    } 
    else 
    { // macro_type == mt_mixture
      Mixture * mixture   = (Mixture *) pData;
      size_t    vec_size  = mixture->mpMean->mVectorSize;
      FLOAT *   obs       = XFormPass(mixture->mpInputXForm,ud->observation,ud->time,FORWARD);
      
      for (i = 0; i < vec_size; i++) 
      {
        mixture->mpMean->mVector[vec_size + i] += obs[i];
      }
      
      mixture->mpMean->mVector[2 * vec_size] += 1;
  
      for (i = 0; i < vec_size; i++) 
      {
        mixture->mpVariance->mVector[vec_size  +i] += SQR(obs[i]);
        mixture->mpVariance->mVector[2*vec_size+i] += obs[i];
      }
      mixture->mpVariance->mVector[3 * vec_size] += 1;
    }
  }
  
  FLOAT *XFormPass(XFormInstance *xformInst, FLOAT *in_vec, int time, PropagDir dir)
  {
    if (xformInst == NULL) return in_vec;
  
    if (time != UNDEF_TIME && xformInst->time == time) return xformInst->out_vec;
  
    xformInst->time = time;
  
    if (xformInst->mpInput)
      in_vec = XFormPass(xformInst->mpInput, in_vec, time, dir);
  
    xformInst->mpXForm->Evaluate(in_vec,xformInst->out_vec,xformInst->memory,dir);
  
    return xformInst->out_vec;
  }
  
  
  void 
  ModelSet::
  Scan(int mask, HMMSetNodeName nodeName,
       ScanAction action, void *pUserData)
  {
    Macro *macro;
    
    if (nodeName != NULL) 
      strcpy(nodeName+sizeof(HMMSetNodeName)-4, "...");
  
    // walk through the list of macros and decide what to do... 
    // we also decide which direction to walk based on the MTM_REVERSE_PASS flag
    for (macro = mask & MTM_REVERSE_PASS ? mpLastMacro : mpFirstMacro;
        macro != NULL;
        macro = mask & MTM_REVERSE_PASS ? macro->prevAll      : macro->nextAll) 
    {
      if (macro->mpData != NULL && macro->mpData->mpMacro != macro) 
        continue;
      
      if (nodeName != NULL) 
        strncpy(nodeName, macro->mpName, sizeof(HMMSetNodeName)-4);
  
      switch (macro->mType) 
      {
        case mt_XForm:
          if (!(mask & MTM_XFORM)) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_XFormInstance:
          if (!(mask & (MTM_XFORM | MTM_XFORM_INSTANCE))) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_mean:
          if (!(mask & MTM_MEAN)) break;
          action(mt_mean, nodeName, macro->mpData, pUserData);
          break;
        case mt_variance:
          if (!(mask & MTM_VARIANCE)) break;
          action(mt_variance, nodeName, macro->mpData, pUserData);
          break;
        case mt_transition:
          if (!(mask & MTM_TRANSITION)) break;
          action(mt_transition, nodeName, macro->mpData, pUserData);
          break;
        case mt_mixture:
          if (!(mask & (MTM_ALL & ~(MTM_STATE | MTM_HMM | MTM_TRANSITION)))) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_state:
          if (!(mask & (MTM_ALL & ~(MTM_HMM | MTM_TRANSITION)))) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_hmm:
          if (!(mask & MTM_ALL)) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        default: assert(0);
      }
    }
  }
  
  void 
  Hmm::
  //ScanHMM(Hmm *hmm, int mask, HMMSetNodeName nodeName,
  //            ScanAction action, void *pUserData)
  Scan(int mask, HMMSetNodeName nodeName,
       ScanAction action, void *pUserData)
  {
    size_t    i;
    size_t    n = 0;
    char *    chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    
    if (mask & MTM_HMM && mask & MTM_PRESCAN) 
      action(mt_hmm, nodeName, this, pUserData);
  
    if (mask & (MTM_ALL & ~(MTM_HMM | MTM_TRANSITION))) 
    {
      for (i=0; i < mNStates-2; i++) 
      {
        if (!state[i]->mpMacro) 
        {
          if (n > 0 ) snprintf(chptr, n, ".state[%d]", i+2);            
          state[i]->Scan(mask, nodeName, action, pUserData);
        }
      }
    }
  
    if (mask & MTM_TRANSITION && !mpTransition->mpMacro) 
    {
      if (n > 0) strncpy(chptr, ".transP", n);
      action(mt_transition, nodeName, mpTransition, pUserData);
    }
  
    if (mask & MTM_HMM && !(mask & MTM_PRESCAN)) 
    {
      if (n > 0) chptr = '\0';
      action(mt_hmm, nodeName, this, pUserData);
    }
  }
  
  void 
  State::
  Scan(int mask, HMMSetNodeName nodeName,  ScanAction action, void *pUserData)
  {
    size_t    i;
    size_t    n = 0;
    char *    chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_STATE && mask & MTM_PRESCAN) 
    {
      action(mt_state, nodeName, this, pUserData);
    }
  
    if (mOutPdfKind != KID_PDFObsVec &&
      mask & (MTM_ALL & ~(MTM_STATE | MTM_HMM | MTM_TRANSITION))) 
    {
      for (i=0; i < mNumberOfMixtures; i++) 
      {
        if (!mpMixture[i].estimates->mpMacro) 
        {
          if (n > 0 ) snprintf(chptr, n, ".mix[%d]", i+1);
          mpMixture[i].estimates->Scan(mask, nodeName,
                      action, pUserData);
        }
      }
    }
    if (mask & MTM_STATE && !(mask & MTM_PRESCAN)) 
    {
      if (n > 0) chptr = '\0';
      action(mt_state, nodeName, this, pUserData);
    }
  }
  
  void 
  Mixture::
  Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *pUserData)
  {
    int n = 0;
    char *chptr = NULL;
    
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_MIXTURE && mask & MTM_PRESCAN) {
      action(mt_mixture, nodeName, this, pUserData);
    }
  
    if (mask & MTM_MEAN && !mpMean->mpMacro) {
      if (n > 0) strncpy(chptr, ".mean", n);
      action(mt_mean, nodeName, mpMean, pUserData);
    }
  
    if (mask & MTM_VARIANCE && !mpVariance->mpMacro) {
      if (n > 0) strncpy(chptr, ".cov", n);
      action(mt_variance, nodeName, mpVariance, pUserData);
    }
  
    if (mask & MTM_XFORM_INSTANCE && mpInputXForm &&
      !mpInputXForm->mpMacro) {
      if (n > 0) strncpy(chptr, ".input", n);
      mpInputXForm->Scan(mask, nodeName, action, pUserData);
    }
  
    if (mask & MTM_MIXTURE && !(mask & MTM_PRESCAN)) {
      if (n > 0) chptr = '\0';
      action(mt_mixture, nodeName, this, pUserData);
    }
  }
  
  void 
  XFormInstance::
  Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *pUserData)
  {
    int n = 0;
    char *chptr = NULL;
    
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_XFORM_INSTANCE && mask & MTM_PRESCAN) {
      action(mt_XFormInstance, nodeName, this, pUserData);
    }
    
    if (mpInput != NULL) {
      if (!mpInput->mpMacro) {
      if (n > 0) strncpy(chptr, ".input", n);
        mpInput->Scan(mask, nodeName, action, pUserData);
      }
    }
  
    if (mask & MTM_XFORM && mpXForm != NULL && 
      !mpXForm->mpMacro) {
      if (n > 0) strncpy(chptr, ".mpXForm", n);
      mpXForm->Scan(mask, nodeName, action, pUserData);
    }
    
    if (mask & MTM_XFORM_INSTANCE && !(mask & MTM_PRESCAN)) {
      if (n > 0) chptr = '\0';
      action(mt_XFormInstance, nodeName, this, pUserData);
    }
  }
  
  void 
  XForm::
  Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *pUserData)
  {
    size_t      n = 0;
    char *      chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_PRESCAN) 
    {
      action(mt_XForm, nodeName, this, pUserData);
    }
    
    if (mXFormType == XT_COMPOSITE) 
    {
      CompositeXForm *  cxf = dynamic_cast<CompositeXForm *> (this);
      size_t            i;
      size_t            j;
  
      for (i=0; i < cxf->mNLayers; i++) 
      {
        for (j = 0; j < cxf->layer[i].mNBlocks; j++) 
        {
          if (!cxf->layer[i].block[j]->mpMacro) 
          {
            if (n > 0) 
              snprintf(chptr, n, ".part[%d,%d]", i+1, j+1);
            cxf->layer[i].block[j]->Scan(mask, nodeName, action, pUserData);
          }
        }
      }
    }
  
    if (!(mask & MTM_PRESCAN)) 
    {
      if (n > 0) 
        chptr = '\0';
        
      action(mt_XForm, nodeName, this, pUserData);
    }
  }
  
  
  enum StatType {MEAN_STATS, COV_STATS};
  
  void AllocXFormStatAccums(XFormStatAccum **xformStatAccum,
                            size_t *nxformStatAccums,
                            XFormInstance *xformInstance,
                            enum StatType stat_type) 
  {
    size_t    i;
    size_t    j;
    
    if (xformInstance == NULL) 
      return;
  
    for (i = 0; i < xformInstance->mNumberOfXFormStatCaches; i++) 
    {
      XFormStatCache *xfsc = &xformInstance->mpXFormStatCache[i];
      XFormStatAccum *xfsa = *xformStatAccum;
  
      for (j = 0; j < *nxformStatAccums; j++, xfsa++) 
      {
        if (xfsa->mpXForm == xfsc->mpXForm) 
          break;
      }
  
      if (j == *nxformStatAccums) 
      {
        size_t size = xfsc->mpXForm->mInSize; //mean : mean+covariance
        size = (stat_type == MEAN_STATS) ? size : size+size*(size+1)/2;
  
        *xformStatAccum =
          (XFormStatAccum *) realloc(*xformStatAccum,
                                    sizeof(XFormStatAccum) * ++*nxformStatAccums);
  
        if (*xformStatAccum == NULL)
          Error("Insufficient memory");
  
        xfsa = *xformStatAccum + *nxformStatAccums - 1;
  
        if ((xfsa->mpStats = (FLOAT *) malloc(sizeof(FLOAT) * size)) == NULL)
          Error("Insufficient memory");
  
        xfsa->mpXForm    = xfsc->mpXForm;
        xfsa->norm     = 0.0;
        
        for (j = 0; j < size; j++) 
          xfsa->mpStats[j] = 0.0;
      }
    }
  }
  
  void AllocateXFormStatCachesAndAccums(int macro_type, HMMSetNodeName nodeName,
                                        MacroData * pData, void * hmm_set)
  {
    if (macro_type == mt_XFormInstance) 
    {
      //Allocate XForm stat caches for XFormInstance
  
      XFormInstance *   xfi =(XFormInstance *) pData;
      size_t            i;
      size_t            j;
  
      for (i=0; i < ((ModelSet *) hmm_set)->mNumberOfXFormsToUpdate; i++) 
      {
        XForm *xform = ((ModelSet *) hmm_set)->mpXFormToUpdate[i].mpXForm;
        int instanceContainXfrom = IsXFormIn1stLayer(xform, xfi->mpXForm);
  
        //Does instance one level up contain cache for this xform
        XFormStatCache *upperLevelStats = NULL;
        if (xfi->mpInput != NULL) {
          for (j=0; j < xfi->mpInput->mNumberOfXFormStatCaches; j++) {
            if (xfi->mpInput->mpXFormStatCache[j].mpXForm == xform) {
              upperLevelStats = &xfi->mpInput->mpXFormStatCache[j];
              break;
            }
          }
        }
  
        if (instanceContainXfrom || upperLevelStats != NULL) 
        {
          XFormStatCache *xfsc;
  
          xfi->mpXFormStatCache = (XFormStatCache *)
            realloc(xfi->mpXFormStatCache,
                    sizeof(XFormStatCache) * ++xfi->mNumberOfXFormStatCaches);
  
          if (xfi->mpXFormStatCache == NULL) {
            Error("Insufficient memory");
          }
  
          xfsc = &xfi->mpXFormStatCache[xfi->mNumberOfXFormStatCaches-1];
  
          if (instanceContainXfrom) {
            int size = xform->mInSize;
            size = size+size*(size+1)/2;
  
            if ((xfsc->mpStats = (FLOAT *) malloc(sizeof(FLOAT) * size))==NULL) {
              Error("Insufficient memory");
            }
          } else {
            xfsc->mpStats = upperLevelStats->mpStats;
          }
  
          xfsc->norm = 0;
          xfsc->mpXForm = xform;
          xfsc->mpUpperLevelStats = upperLevelStats;
        }
      }
    } 
    
    else if (macro_type == mt_mixture) 
    {
      //Allocate XForm stat accumulators for mean and covariance
  
      Mixture *mix = (Mixture *) pData;
      AllocXFormStatAccums(&mix->mpMean->mpXFormStatAccum,
                          &mix->mpMean->mNumberOfXFormStatAccums,
                          mix->mpInputXForm, MEAN_STATS);
  
      AllocXFormStatAccums(&mix->mpVariance->mpXFormStatAccum,
                          &mix->mpVariance->mNumberOfXFormStatAccums,
                          mix->mpInputXForm, COV_STATS);
  
      if (mix->mpInputXForm == NULL || mix->mpInputXForm->mNumberOfXFormStatCaches == 0)
        return;
        
      if (mix->mpInputXForm->mNumberOfXFormStatCaches != 1 ||
        !Is1Layer1BlockLinearXForm(mix->mpInputXForm->mpXForm) ||
        mix->mpInputXForm->mpXFormStatCache[0].mpUpperLevelStats != NULL) 
      {
        mix->mpVariance->mUpdatableFromStatAccums = false;
        mix->mpMean    ->mUpdatableFromStatAccums = false;
        ((ModelSet *) hmm_set)->mAllMixuresUpdatableFromStatAccums = false;  
      } 
      
      else if (mix->mpMean->mNumberOfXFormStatAccums != 1) 
      {
        assert(mix->mpMean->mNumberOfXFormStatAccums > 1);
        mix->mpMean->mUpdatableFromStatAccums = false;
        ((ModelSet *) hmm_set)->mAllMixuresUpdatableFromStatAccums = false;
      } 
      
      else if (mix->mpVariance->mNumberOfXFormStatAccums != 1) 
      {
        assert(mix->mpVariance->mNumberOfXFormStatAccums > 1);
        mix->mpVariance->mUpdatableFromStatAccums = false;
        ((ModelSet *) hmm_set)->mAllMixuresUpdatableFromStatAccums = false;
      }
    }
  }
  
  
  void NormalizeStatsForXForm(int macro_type, HMMSetNodeName nodeName,
                              MacroData * pData, void * pUserData) 
  {
    XFormStatAccum *xfsa = NULL;
    int i, j, k, nxfsa = 0, size;
    FLOAT *mean, *cov, inorm;
  
    if (macro_type == mt_mean) 
    {
      xfsa  = ((Mean *)pData)->mpXFormStatAccum;
      nxfsa = ((Mean *)pData)->mNumberOfXFormStatAccums;
    } 
    
    else if (macro_type == mt_variance) 
    {
      xfsa  = ((Variance *)pData)->mpXFormStatAccum;
      nxfsa = ((Variance *)pData)->mNumberOfXFormStatAccums;
    }
  
    for (i = 0; i < nxfsa; i++) 
    {
      size = xfsa[i].mpXForm->mInSize;
      mean = xfsa[i].mpStats;
      cov  = xfsa[i].mpStats + size;
      inorm = 1.0 / xfsa[i].norm;
  
      for (j = 0; j < size; j++) 
        mean[j] *= inorm; //normalize means
  
      if (macro_type == mt_variance) 
      {
        for (k=0; k < size; k++) 
        {
          for (j=0; j <= k; j++) 
          {                 //normalize covariances
            cov[k*(k+1)/2+j] = cov[k*(k+1)/2+j] * inorm - mean[k] * mean[j];
          }
        }
      }
    }
  }
  
  
  void WriteStatsForXForm(int macro_type, HMMSetNodeName nodeName,
                          MacroData * pData, void * pUserData) 
  {
    XFormStatsFileNames *    file = NULL;
    XFormStatAccum *              xfsa = NULL;
    size_t                        i;
    size_t                        j;
    size_t                        k;
    size_t                        nxfsa = 0;
    int                           cc = 0;
    size_t                        size;
    FLOAT *                       mean;
    FLOAT *                       cov;
    WriteStatsForXFormUserData * ud = (WriteStatsForXFormUserData *) pUserData;
  
    if (macro_type == mt_mean) 
    {
      file  = &ud->mMeanFile;
      xfsa  = ((Mean *)pData)->mpXFormStatAccum;
      nxfsa = ((Mean *)pData)->mNumberOfXFormStatAccums;
    } 
    else if (macro_type == mt_variance) 
    {
      file  = &ud->mCovFile;
      xfsa  = ((Variance *)pData)->mpXFormStatAccum;
      nxfsa = ((Variance *)pData)->mNumberOfXFormStatAccums;
    }
  
    for (i = 0; i < nxfsa && xfsa[i].mpXForm != (XForm *) ud->mpXForm; i++)
      ;
    
    if (i == nxfsa) 
      return;
  
    if (fprintf(file->mpOccupP, "%s "FLOAT_FMT"\n", nodeName, xfsa[i].norm) < 0) 
    {
      Error("Cannot write to file: %s", file->mpOccupN);
    }
  
    size = xfsa[i].mpXForm->mInSize;
    mean = xfsa[i].mpStats;
    cov  = xfsa[i].mpStats + size;
  
    if (macro_type == mt_mean) 
    {
      if (ud->mBinary) 
      {
        if (!isBigEndian()) 
          for (i = 0; i < size; i++) swapFLOAT(mean[i]);
          
        cc |= (fwrite(mean, sizeof(FLOAT), size, file->mpStatsP) != size);
        
        if (!isBigEndian()) 
          for (i = 0; i < size; i++) swapFLOAT(mean[i]);
      } 
      else 
      {
        for (j=0;j<size;j++) 
        {
          cc |= fprintf(file->mpStatsP, FLOAT_FMT" ", mean[j]) < 0;
        }
        
        cc |= fputs("\n", file->mpStatsP) == EOF;
      }
    } 
    else 
    {
      if (ud->mBinary) 
      {
        size = size*(size+1)/2;
        if (!isBigEndian()) for (i = 0; i < size; i++) swapFLOAT(cov[i]);
        cc |= fwrite(cov, sizeof(FLOAT), size, file->mpStatsP) != size;
        if (!isBigEndian()) for (i = 0; i < size; i++) swapFLOAT(cov[i]);
      } 
      else
      {
        for (k=0; k < size; k++) 
        {
          for (j=0;j<=k;j++) {
            cc |= fprintf(file->mpStatsP, FLOAT_FMT" ", cov[k*(k+1)/2+j]) < 0;
          }
  
          for (;j<size; j++) {
            cc |= fprintf(file->mpStatsP, FLOAT_FMT" ", cov[j*(j+1)/2+k]) < 0;
          }
  
          cc |= fputs("\n", file->mpStatsP) == EOF;
        }
        cc |= fputs("\n", file->mpStatsP) == EOF;
      }
    }
  
    if (cc) {
      Error("Cannot write to file %s", file->mpStatsN);
    }
  }
  
  
  void ReadStatsForXForm(int macro_type, HMMSetNodeName nodeName,
                        MacroData * pData, void *pUserData) 
  {
    char                        buff[128];
    XFormStatsFileNames *  file = NULL;
    XFormStatAccum *            xfsa = NULL;
    size_t                      i;
    size_t                      j;
    size_t                      k;
    size_t                      nxfsa = 0;
    int                         cc = 0;
    size_t                      size;
    FLOAT *                     mean;
    FLOAT *                     cov;
    FLOAT                       f;
    WriteStatsForXFormUserData * ud = (WriteStatsForXFormUserData *) pUserData;
  
    if (macro_type == mt_mean) 
    {
      file  = &ud->mMeanFile;
      xfsa  = ((Mean *)pData)->mpXFormStatAccum;
      nxfsa = ((Mean *)pData)->mNumberOfXFormStatAccums;
    } 
    else if (macro_type == mt_variance) 
    {
      file  = &ud->mCovFile;
      xfsa  = ((Variance *)pData)->mpXFormStatAccum;
      nxfsa = ((Variance *)pData)->mNumberOfXFormStatAccums;
    }
  
    for (i = 0; i < nxfsa && xfsa[i].mpXForm != (XForm *) ud->mpXForm; i++)
      ;
      
    if (i == nxfsa) 
      return;
  
    j = fscanf(file->mpOccupP, "%128s "FLOAT_FMT"\n", buff, &xfsa[i].norm);
    
    if (j < 1) 
      Error("Unexpected end of file: %s", file->mpOccupN);
    
    else if (strcmp(buff, nodeName)) 
      Error("'%s' expected but '%s' found in file: %s",nodeName,buff,file->mpOccupN);
    
    else if (j < 2) 
      Error("Decimal number expected after '%s'in file: %s", buff, file->mpOccupN);
    
    size = xfsa[i].mpXForm->mInSize;
    mean = xfsa[i].mpStats;
    cov  = xfsa[i].mpStats + size;
  
    if (macro_type == mt_mean) 
    {
      if (ud->mBinary) 
      {
        j = fread(mean, sizeof(FLOAT), size, file->mpStatsP);
        cc |= j != size;
        if (!isBigEndian()) for (i = 0; i < size; i++) swapFLOAT(mean[i]);
      } 
      else 
      {
        for (j=0;j<size;j++)
        {
          cc |= (fscanf(file->mpStatsP, FLOAT_FMT" ", &mean[j]) != 1);
        }
      }
    } 
    else 
    {
      if (ud->mBinary) 
      {
        size = size*(size+1)/2;
        cc |= (fread(cov, sizeof(FLOAT), size, file->mpStatsP) != size);
        
        if (!isBigEndian()) 
          for (i = 0; i < size; i++) swapFLOAT(cov[i]);
      } 
      else
      {
        for (k=0; k < size; k++) 
        {
          for (j=0;j<k;j++) 
          {
            cc |= (fscanf(file->mpStatsP, FLOAT_FMT" ", &f) != 1);
            if (f != cov[k*(k+1)/2+j]) 
            {
              Error("Covariance matrix '%s' in file '%s' must be symetric",
                    nodeName, file->mpStatsP);
            }
          }
  
          for (;j<size; j++) {
            cc |= (fscanf(file->mpStatsP, FLOAT_FMT" ", &cov[j*(j+1)/2+k]) != 1);
          }
        }
      }
    }
  
    if (ferror(file->mpStatsP)) 
    {
      Error("Cannot read file '%s'", file->mpStatsN);
    } 
    else if (cc) 
    {
      Error("Invalid file with XForm statistics '%s'", file->mpStatsN);
    }
  }
  
  
  void WriteAccum(int macro_type, HMMSetNodeName nodeName,
                  MacroData * pData, void *pUserData) 
  {
    size_t                i;
    size_t                size;
    FLOAT *               vector = NULL;
  //  FILE *fp = (FILE *) pUserData;
    WriteAccumUserData *  ud = (WriteAccumUserData * ) pUserData;
    Macro *               macro = static_cast <MacroData *> (pData)->mpMacro;
  
    if (macro &&
      (fprintf(ud->fp, "~%c \"%s\"", macro->mType, macro->mpName) < 0 ||
      fwrite(&macro->mOccurances, sizeof(macro->mOccurances), 1, ud->fp) != 1)) 
    {
      Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
    }
  
    if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      XFormStatAccum *    xfsa = NULL;
      size_t              nxfsa = 0;
  
      if (macro_type == mt_mean) 
      {
        xfsa   = ((Mean *)pData)->mpXFormStatAccum;
        nxfsa  = ((Mean *)pData)->mNumberOfXFormStatAccums;
        size   = ((Mean *)pData)->mVectorSize;
        vector = ((Mean *)pData)->mVector+size;
        size   = size + 1;
      } 
      else if (macro_type == mt_variance) 
      {
        xfsa   = ((Variance *)pData)->mpXFormStatAccum;
        nxfsa  = ((Variance *)pData)->mNumberOfXFormStatAccums;
        size   = ((Variance *)pData)->mVectorSize;
        vector = ((Variance *)pData)->mVector+size;
        size   = size * 2 + 1;
      }
  
  //    if (ud->mMmi) vector += size; // Move to MMI accums, which follows ML accums
  
      if (fwrite(vector, sizeof(FLOAT), size, ud->fp) != size ||
        fwrite(&nxfsa, sizeof(nxfsa),    1, ud->fp) != 1) 
      {
        Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
      }
  
  //    if (!ud->mMmi) { // MMI estimation of XForm statistics has not been implemented yet
      for (i = 0; i < nxfsa; i++) 
      {
        size = xfsa[i].mpXForm->mInSize;
        size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
        assert(xfsa[i].mpXForm->mpMacro != NULL);
        if (fprintf(ud->fp, "\"%s\"", xfsa[i].mpXForm->mpMacro->mpName) < 0 ||
          fwrite(&size,         sizeof(int),      1, ud->fp) != 1    ||
          fwrite(xfsa[i].mpStats, sizeof(FLOAT), size, ud->fp) != size ||
          fwrite(&xfsa[i].norm, sizeof(FLOAT),    1, ud->fp) != 1) 
        {
          Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
        }
      }
  //    }
    }
    
    else if (macro_type == mt_state) 
    {
      State *state = (State *) pData;
      if (state->mOutPdfKind == KID_DiagC) 
      {
        for (i = 0; i < state->mNumberOfMixtures; i++) 
        {
          if (fwrite(&state->mpMixture[i].weight_accum,
                    sizeof(FLOAT), 1, ud->fp) != 1 ||
            fwrite(&state->mpMixture[i].weight_accum_den,
                    sizeof(FLOAT), 1, ud->fp) != 1) 
          {
            Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
          }
        }
      }
    } 
    
    else if (macro_type == mt_transition) 
    {
      size   = SQR(((Transition *) pData)->mNStates);
      vector = ((Transition *) pData)->matrix + size;
  
      if (fwrite(vector, sizeof(FLOAT), size, ud->fp) != size) 
      {
        Error("Cannot write accumulators to file: '%s'", ud->mpFileName);
      }
    }
  }
  
  
  void NormalizeAccum(int macro_type, HMMSetNodeName nodeName,
                    MacroData * pData, void *pUserData) 
  {
    size_t      i;
    size_t      j;
    size_t      size;
    FLOAT *     vector = NULL;
  
    if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      XFormStatAccum *  xfsa = NULL;
      size_t            nxfsa = 0;
  
      if (macro_type == mt_mean) 
      {
        xfsa   = ((Mean *)pData)->mpXFormStatAccum;
        nxfsa  = ((Mean *)pData)->mNumberOfXFormStatAccums;
        size   = ((Mean *)pData)->mVectorSize;
        vector = ((Mean *)pData)->mVector+size;
        size   = size + 1;
      } 
      else if (macro_type == mt_variance) 
      {
        xfsa   = ((Variance *)pData)->mpXFormStatAccum;
        nxfsa  = ((Variance *)pData)->mNumberOfXFormStatAccums;
        size   = ((Variance *)pData)->mVectorSize;
        vector = ((Variance *)pData)->mVector+size;
        size   = size * 2 + 1;
      }
  
      for (i=0; i < size; i++) 
        vector[i] /= vector[size-1];
  
      for (i = 0; i < nxfsa; i++) 
      {
        size = xfsa[i].mpXForm->mInSize;
        size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
  
        for (j=0; j < size; j++) 
          xfsa[i].mpStats[j] /= xfsa[i].norm;
        
        xfsa[i].norm = 1.0;
      }
    } 
    else if (macro_type == mt_state) 
    {
      State *state = (State *) pData;
      
      if (state->mOutPdfKind == KID_DiagC) 
      {
        FLOAT accum_sum = 0.0;
  
        for (i = 0; i < state->mNumberOfMixtures; i++)
          accum_sum += state->mpMixture[i].weight_accum;
  
        if (accum_sum > 0.0) 
        {
          for (i = 0; i < state->mNumberOfMixtures; i++)
            state->mpMixture[i].weight_accum /= accum_sum;
        }
      }
    } 
    else if (macro_type == mt_transition) 
    {
      size_t nstates = ((Transition *) pData)->mNStates;
      vector = ((Transition *) pData)->matrix + SQR(nstates);
  
      for (i=0; i < nstates; i++) 
      {
        FLOAT nrm = LOG_0;
        
        for (j=0; j < nstates; j++) 
        {
          LOG_INC(nrm, vector[i * nstates + j]);
        }
        
        if (nrm < LOG_MIN) nrm = 0.0;
        
        for (j=0; j < nstates; j++) 
        {
          vector[i * nstates + j] -= nrm;
        }
      }
    }
  }
  
  
  unsigned int faddfloat(FLOAT *vec, size_t size, float mul_const, FILE *fp) 
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
  
  
  void ReadAccum(int macro_type, HMMSetNodeName nodeName,
                  MacroData * pData, void *pUserData) 
  {
    unsigned int        i;
    unsigned int        j;
    int                 c;
    unsigned int        size;
    FLOAT *             vector = NULL;
    Macro *             macro;
    ReadAccumUserData * ud = (ReadAccumUserData *) pUserData;
  //  FILE *fp =        ((ReadAccumUserData *) pUserData)->fp;
  //  char *fn =        ((ReadAccumUserData *) pUserData)->fn;
  //  ModelSet *hmm_set = ((ReadAccumUserData *) pUserData)->hmm_set;
  //  float weight =    ((ReadAccumUserData *) pUserData)->weight;
    char                xfName[128];
  
    xfName[sizeof(xfName)-1] = '\0';
  
    if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      XFormStatAccum *  xfsa = NULL;
      size_t            size_inf;
      size_t            nxfsa_inf;
      size_t            nxfsa = 0;
  
      if (macro_type == mt_mean) 
      {
        xfsa   = ((Mean *)pData)->mpXFormStatAccum;
        nxfsa  = ((Mean *)pData)->mNumberOfXFormStatAccums;
        size   = ((Mean *)pData)->mVectorSize;
        vector = ((Mean *)pData)->mVector+size;
        size   = size + 1;
      } 
      else if (macro_type == mt_variance) 
      {
        xfsa   = ((Variance *)pData)->mpXFormStatAccum;
        nxfsa  = ((Variance *)pData)->mNumberOfXFormStatAccums;
        size   = ((Variance *)pData)->mVectorSize;
        vector = ((Variance *)pData)->mVector+size;
        size   = size * 2 + 1;
      }
  
      if (ud->mMmi) 
        vector += size;
  
      if (faddfloat(vector, size, ud->mWeight,     ud->fp) != size ||
        fread(&nxfsa_inf, sizeof(nxfsa_inf), 1, ud->fp) != 1) 
      {
        Error("Incompatible accumulator file: '%s'", ud->mpFileName);
      }
  
      if (!ud->mMmi) { // MMI estimation of XForm statistics has not been implemented yet
        for (i = 0; i < nxfsa_inf; i++) 
        {
          if (getc(ud->fp) != '"') 
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
          for (j=0; (c=getc(ud->fp)) != EOF && c != '"' && j < sizeof(xfName)-1; j++) 
            xfName[j] = c;
  
          xfName[j] = '\0';
          
          if (c == EOF)
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
          macro = FindMacro(&ud->mpModelSet->mXFormHash, xfName);
  
          if (fread(&size_inf, sizeof(int), 1, ud->fp) != 1) 
            Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
          if (macro != NULL) 
          {
            size = ((LinearXForm *) macro->mpData)->mInSize;
            size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
  
            if (size != size_inf)
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
  
            for (j = 0; j < nxfsa && xfsa[j].mpXForm != macro->mpData; j++)
              ;
            
            if (j < nxfsa) 
            {
              if (faddfloat(xfsa[j].mpStats, size, ud->mWeight, ud->fp) != size  ||
                  faddfloat(&xfsa[j].norm,    1, ud->mWeight, ud->fp) != 1) 
              {
                Error("Invalid accumulator file: '%s'", ud->mpFileName);
              }
            } 
            else 
            {
              macro = NULL;
            }
          }
  
          // Skip XForm accumulator
          if (macro == NULL) 
          { 
            FLOAT f;
            for (j = 0; j < size_inf+1; j++) 
              fread(&f, sizeof(f), 1, ud->fp);
          }
        }
      }
    } else if (macro_type == mt_state) {
      State *state = (State *) pData;
      if (state->mOutPdfKind == KID_DiagC) {
        FLOAT junk;
        for (i = 0; i < state->mNumberOfMixtures; i++) {
          if (ud->mMmi == 1) {
            if (faddfloat(&state->mpMixture[i].weight_accum_den, 1, ud->mWeight, ud->fp) != 1 ||
              faddfloat(&junk,                               1, ud->mWeight, ud->fp) != 1) {
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
            }
          } else if (ud->mMmi == 2) {
            if (faddfloat(&junk,                               1, ud->mWeight, ud->fp) != 1 ||
              faddfloat(&state->mpMixture[i].weight_accum_den, 1, ud->mWeight, ud->fp) != 1) {
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
            }
          } else {
            if (faddfloat(&state->mpMixture[i].weight_accum,     1, ud->mWeight, ud->fp) != 1 ||
              faddfloat(&state->mpMixture[i].weight_accum_den, 1, ud->mWeight, ud->fp) != 1) {
              Error("Incompatible accumulator file: '%s'", ud->mpFileName);
            }
          }
        }
      }
    } else if (macro_type == mt_transition) {
      FLOAT f;
      size   = SQR(((Transition *) pData)->mNStates);
      vector =     ((Transition *) pData)->matrix + size;
  
      for (i = 0; i < size; i++) {
        if (fread(&f, sizeof(FLOAT), 1, ud->fp) != 1) {
        Error("Incompatible accumulator file: '%s'", ud->mpFileName);
        }
        if (!ud->mMmi) { // MMI estimation of transition probabilities has not been implemented yet
          f += log(ud->mWeight);
          LOG_INC(vector[i], f);
        }
      }
    }
  }
  

  //**************************************************************************
  //**************************************************************************
  //   Hmm class section
  //**************************************************************************
  //**************************************************************************
  Hmm::
  Hmm(size_t nStates) : 
    mNStates(0), mpTransition(NULL), mpState(NULL)
  {
    // we allocate pointers for states. The -2 is for the non-emmiting states
    mpState = new State* [mNStates - 2];
    mNStates = nStates;
  }
  
  
  //**************************************************************************
  Hmm::
  ~Hmm()
  {
    size_t  i;
    
    // we go through each state and check whether it is pointed to by 
    // some macro
    if (mpState != NULL)
    {
      for (i = 0; i < mNStates; i++)
      {
        if ((static_cast<State *>(mpState[i]))->mpMacro == NULL)
          delete mpState[i];
      }
      
      delete [] mpState;
    }
    
    // delete transition object if not pointed to by macro
    if (mpTransition != NULL)
    {
      if (mpTransition->mpMacro == NULL)
        delete mpTransition;
    }
  }
  
  
  //**************************************************************************  
  void
  Hmm::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    size_t i;
  
    for (i = 0; i < mNStates - 2; i++) 
    {
      if (!this->state[i]->mpMacro) 
        state[i]->UpdateFromAccums(pModelSet, this);
    }
  
    if (!mpTransition->mpMacro) 
      mpTransition->UpdateFromAccums(pModelSet);
  } // UpdateFromAccums(const ModelSet * pModelSet)

  
  //**************************************************************************
  //**************************************************************************
  //   Mixture section
  //**************************************************************************
  //**************************************************************************
  void 
  Mixture::
  ComputeGConst()
  {
    FLOAT cov_det = 0;
    size_t i;
  
    for (i = 0; i < this->mpVariance->mVectorSize; i++) 
    {
      cov_det -= log(this->mpVariance->mVector[i]);
    }
    this->mGConst = cov_det + M_LOG_2PI * this->mpVariance->mVectorSize;
  }
  

  //**************************************************************************
  void 
  Mixture::
  UpdateFromAccums(const ModelSet * pModelSet) 
  {
    double  Djm;
    int     i;

    if (pModelSet->mMmiUpdate == 1 || pModelSet->mMmiUpdate == -1) 
    {
      int vec_size    = mpVariance->mVectorSize;
      FLOAT *mean_vec = mpMean->mVector;
      FLOAT *var_vec  = mpVariance->mVector;
  
      FLOAT *vac_num  = var_vec + 1 * vec_size;
      FLOAT *mac_num  = var_vec + 2 * vec_size;
      FLOAT *nrm_num  = var_vec + 3 * vec_size;
  
      FLOAT *vac_den  = vac_num + 2 * vec_size + 1;
      FLOAT *mac_den  = mac_num + 2 * vec_size + 1;
      FLOAT *nrm_den  = nrm_num + 2 * vec_size + 1;
  
      if (pModelSet->mMmiUpdate == 1) 
      {
        // I-smoothing
        for (i = 0; i < vec_size; i++) {
          mac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
          vac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
        }
        *nrm_num   += pModelSet->MMI_tauI;
      }
  
      Djm = 0.0;
      // Find minimum Djm leading to positive update of variances
      for (i = 0; i < vec_size; i++) {
        double macn_macd = mac_num[i]-mac_den[i];
        double vacn_vacd = vac_num[i]-vac_den[i];
        double nrmn_nrmd = *nrm_num - *nrm_den;
        double a  = 1/var_vec[i];
        double b  = vacn_vacd + nrmn_nrmd * (1/var_vec[i] + SQR(mean_vec[i])) -
                  2 * macn_macd * mean_vec[i];
        double c  = nrmn_nrmd * vacn_vacd - SQR(macn_macd);
        double Dd = (- b + sqrt(SQR(b) - 4 * a * c)) / (2 * a);
  
        Djm = HIGHER_OF(Djm, Dd);
      }
  
      Djm = HIGHER_OF(pModelSet->MMI_h * Djm, pModelSet->MMI_E * *nrm_den);
  
      if (pModelSet->mMmiUpdate == -1) {
        // I-smoothing
        for (i = 0; i < vec_size; i++) {
          mac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
          vac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
        }
        *nrm_num   += pModelSet->MMI_tauI;
      }
      
      for (i = 0; i < vec_size; i++) {
        double macn_macd = mac_num[i]-mac_den[i];
        double vacn_vacd = vac_num[i]-vac_den[i];
        double nrmn_nrmd = *nrm_num - *nrm_den;
  
        double new_mean = (macn_macd + Djm * mean_vec[i]) / (nrmn_nrmd + Djm);
        var_vec[i]     = 1/((vacn_vacd + Djm * (1/var_vec[i] + SQR(mean_vec[i]))) /
                        (nrmn_nrmd + Djm) - SQR(new_mean));
        mean_vec[i]    = new_mean;
  
        if (pModelSet->mpVarFloor) {
          var_vec[i] = LOWER_OF(var_vec[i], pModelSet->mpVarFloor->mVector[i]);
        }
      }
    } else if (pModelSet->mMmiUpdate == 2 || pModelSet->mMmiUpdate == -2 ) { // MFE update
      int vec_size    = mpVariance->mVectorSize;
      FLOAT *mean_vec = mpMean->mVector;
      FLOAT *var_vec  = mpVariance->mVector;
  
      FLOAT *vac_mle  = var_vec + 1 * vec_size;
      FLOAT *mac_mle  = var_vec + 2 * vec_size;
      FLOAT *nrm_mle  = var_vec + 3 * vec_size;
  
      FLOAT *vac_mfe  = vac_mle + 2 * vec_size + 1;
      FLOAT *mac_mfe  = mac_mle + 2 * vec_size + 1;
      FLOAT *nrm_mfe  = nrm_mle + 2 * vec_size + 1;
  
      if (pModelSet->mMmiUpdate == 2) {
        // I-smoothing
        for (i = 0; i < vec_size; i++) {
          mac_mfe[i] += (pModelSet->MMI_tauI / *nrm_mle * mac_mle[i]);
          vac_mfe[i] += (pModelSet->MMI_tauI / *nrm_mle * vac_mle[i]);
        }
        *nrm_mfe += pModelSet->MMI_tauI;
      }
      
      Djm = 0.0;
      // Find minimum Djm leading to positive update of variances
      for (i = 0; i < vec_size; i++) {
        double macn_macd = mac_mfe[i];
        double vacn_vacd = vac_mfe[i];
        double nrmn_nrmd = *nrm_mfe;
        double a  = 1/var_vec[i];
        double b  = vacn_vacd + nrmn_nrmd * (1/var_vec[i] + SQR(mean_vec[i])) -
                  2 * macn_macd * mean_vec[i];
        double c  = nrmn_nrmd * vacn_vacd - SQR(macn_macd);
        double Dd = (- b + sqrt(SQR(b) - 4 * a * c)) / (2 * a);
  
        Djm = HIGHER_OF(Djm, Dd);
      }
  
  //    weight_accum_den is passed using quite ugly hack that work
  //    only if mixtures are not shared by more states - MUST BE REWRITEN
      Djm = HIGHER_OF(pModelSet->MMI_h * Djm, pModelSet->MMI_E * weight_accum_den);
  
      if (pModelSet->mMmiUpdate == -2) {
        // I-smoothing
        for (i = 0; i < vec_size; i++) {
          mac_mfe[i] += (pModelSet->MMI_tauI / *nrm_mle * mac_mle[i]);
          vac_mfe[i] += (pModelSet->MMI_tauI / *nrm_mle * vac_mle[i]);
        }
        *nrm_mfe += pModelSet->MMI_tauI;
      }
  
      for (i = 0; i < vec_size; i++) {
        double macn_macd = mac_mfe[i];
        double vacn_vacd = vac_mfe[i];
        double nrmn_nrmd = *nrm_mfe;
  
        double new_mean = (macn_macd + Djm * mean_vec[i]) / (nrmn_nrmd + Djm);
        var_vec[i]     = 1/((vacn_vacd + Djm * (1/var_vec[i] + SQR(mean_vec[i]))) /
                        (nrmn_nrmd + Djm) - SQR(new_mean));
        mean_vec[i]    = new_mean;
  
        if (pModelSet->mpVarFloor) {
          var_vec[i] = LOWER_OF(var_vec[i], pModelSet->mpVarFloor->mVector[i]);
        }
      }
    } 
    else 
    {
      if (!mpVariance->mpMacro)
        mpVariance->UpdateFromAccums(pModelSet);
  
      if (!mpMean->mpMacro)
        mpMean->UpdateFromAccums(pModelSet);
    }
    ComputeGConst();
  }; // UpdateFromAccums(const ModelSet * pModelSet) 

      
  //**************************************************************************
  //**************************************************************************
  //   ModelSet section
  //**************************************************************************
  //**************************************************************************


  
  //**************************************************************************  
  //**************************************************************************  
  // Mean section
  //**************************************************************************  
  //**************************************************************************  
  void
  Mean::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    int      i;
    FLOAT *  vec  = mVector;
    FLOAT *  acc  = mVector + 1 * mVectorSize;
    FLOAT    nrm  = mVector  [2 * mVectorSize];
  
    if (pModelSet->mUpdateMask & UM_MEAN) 
    {
      if (mNumberOfXFormStatAccums == 0) 
      {
        for (i = 0; i < mVectorSize; i++) 
        {
          vec[i] = acc[i] / nrm;
        }
      } 
      else 
      { // Updating from xform statistics
        int r;
        int c;
        int pos = 0;
        
        if (!mUpdatableFromStatAccums) 
          return;
        
        // If 'mUpdatableFromStatAccums' is true then 'mNumberOfXFormStatAccums' is equal to 1
        for (i=0; i < mNumberOfXFormStatAccums; i++) 
        {
          LinearXForm * xform   = (LinearXForm *) mpXFormStatAccum[i].mpXForm;
          int           in_size = xform->mInSize;
          FLOAT *       mnv     = mpXFormStatAccum[i].mpStats;
  
          assert(xform->mXFormType == XT_LINEAR);
          assert(pos + xform->mOutSize <= mVectorSize);
          
          for (r = 0; r < xform->mOutSize; r++) 
          {
            mVector[pos + r] = 0;
            for (c = 0; c < in_size; c++) 
            {
              mVector[pos + r] += mnv[c] * xform->matrix[in_size * r + c];
            }
          }
          pos += xform->mOutSize;
        }
        
        assert(pos == mVectorSize);
      }
    }
  } // UpdateFromAccums(const ModelSet * pModelSet)

  
  //**************************************************************************  
  //**************************************************************************  
  // Variance section
  //**************************************************************************  
  //**************************************************************************  
  void
  Variance::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    size_t    i;
  
    if (pModelSet->mUpdateMask & UM_VARIANCE) 
    {
      if (mNumberOfXFormStatAccums == 0) 
      {
        FLOAT *vec  = mVector;
        FLOAT *vac  = mVector + 1 * mVectorSize;
        FLOAT *mac  = mVector + 2 * mVectorSize;
        FLOAT *nrm  = mVector + 3 * mVectorSize;
        
        for (i = 0; i < mVectorSize; i++) 
        {
          if (pModelSet->mUpdateMask & UM_OLDMEANVAR) 
          {
            vec[i] = *nrm / vac[i];
          } 
          else 
          {
            vec[i] = 1 / (vac[i] / *nrm - SQR(mac[i] / *nrm));
          }
          
          if (pModelSet->mpVarFloor && 
              pModelSet->mpVarFloor->mVectorSize == mVectorSize) 
          {
            vec[i] = LOWER_OF(vec[i], pModelSet->mpVarFloor->mVector[i]);
          }
        }
      } 
      else // Updating from xform statistics
      { 
        size_t r;
        size_t c;
        size_t t;
        size_t pos = 0;
        
        if (!mUpdatableFromStatAccums) 
          return;
  
  //    If 'mUpdatableFromStatAccums' is true then 'mNumberOfXFormStatAccums' is equal to 1
        for (i=0; i<mNumberOfXFormStatAccums; i++) 
        {
          LinearXForm * xform = (LinearXForm *) mpXFormStatAccum[i].mpXForm;
          size_t        in_size = xform->mInSize;
          FLOAT *       cov  = mpXFormStatAccum[i].mpStats + in_size;
  
          assert(xform->mXFormType == XT_LINEAR);
          assert(pos + xform->mOutSize <= mVectorSize);
          
          for (r = 0; r < xform->mOutSize; r++) 
          {
            mVector[pos + r] = 0.0;
  
            for (c = 0; c < in_size; c++) 
            {
              FLOAT aux = 0;
              for (t = 0; t <= c; t++) 
              {
                aux += cov[c * (c+1)/2 + t]    * xform->matrix[in_size * r + t];
              }
  
              for (; t < in_size; t++) 
              {
                aux += cov[t * (t+1)/2 + c]    * xform->matrix[in_size * r + t];
              }
              
              mVector[pos + r] += aux * xform->matrix[in_size * r + c];
            }
            mVector[pos + r] = 1 / mVector[pos + r];
          }
          pos += xform->mOutSize;
        }
        assert(pos == mVectorSize);
      }
    }
  }

  
  //**************************************************************************  
  //**************************************************************************  
  // Transition section
  //**************************************************************************  
  //**************************************************************************  
  void
  Transition::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    int       i;
    int       j;
    int       nstates = mNStates;
    FLOAT *   vec = matrix;
    FLOAT *   acc = matrix + 1 * SQR(mNStates);
  
    if (pModelSet->mUpdateMask & UM_TRANSITION) 
    {
      for (i=0; i < nstates; i++) 
      {
        FLOAT nrm = LOG_0;
        for (j=0; j < nstates; j++) 
        {
          LOG_INC(nrm, acc[i * nstates + j]);
        }
        
        if (nrm == LOG_0) 
          nrm = 0;
        
        for (j=0; j < nstates; j++) 
        {
          vec[i * nstates + j] = acc[i * nstates + j] - nrm; // it is in log
        }
      }
    }
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  // Transition section
  //**************************************************************************  
  //**************************************************************************  
  void 
  State::
  UpdateFromAccums(const ModelSet * pModelSet, const Hmm * pHmm) 
  {
    size_t i;
  
    if (mOutPdfKind == KID_DiagC) 
    {
      FLOAT accum_sum = 0;
  
//    if (hmm_set->mUpdateMask & UM_WEIGHT) {
      for (i = 0; i < mNumberOfMixtures; i++) 
      {
        accum_sum += mpMixture[i].weight_accum;
      }

      if (accum_sum <= 0.0) 
      {
        if (mpMacro) 
        {
          Warning("No occupation of '%s', state is not updated",
                  mpMacro->mpName);
        } 
        else 
        {
          size_t j; // find the state number
          for (j=0; j < pHmm->mNStates && pHmm->state[j] != this; j++)
            ;
          Warning("No occupation of '%s[%d]', state is not updated",
                  pHmm->mpMacro->mpName, j + 1);
        }
        return;
      }
  
        // Remove mixtures with low weight
      if (pModelSet->mUpdateMask & UM_WEIGHT) 
      {
        for (i = 0; i < mNumberOfMixtures; i++) 
        {
          if (mpMixture[i].weight_accum / accum_sum < pModelSet->mMinMixWeight) 
          {
            if (pHmm) 
            {
              size_t j;
              for (j=0; j < pHmm->mNStates && pHmm->state[j] != this; j++)
                ; // find the state number
              Warning("Discarding mixture %d of Hmm %s state %d because of too low mixture weight",
                      i, pHmm->mpMacro->mpName, j + 1);
            } 
            else 
            {
              assert(mpMacro);
              Warning("Discarding mixture %d of state %s because of too low mixture weight",
                      i, mpMacro->mpName);
            }
            accum_sum -= mpMixture[i].weight_accum;
  
            if (!mpMixture[i].estimates->mpMacro) 
            {
              mpMixture[i].estimates->Scan(MTM_ALL,NULL,ReleaseItem,NULL);
            }
  
            mpMixture[i--] = mpMixture[--mNumberOfMixtures];
            continue;
          }
        }
      }
  
      for (i = 0; i < mNumberOfMixtures; i++) 
      {
  //      printf("Weight Acc: %f\n", (float) state->mpMixture[i].weight_accum);
        if (pModelSet->mUpdateMask & UM_WEIGHT) 
        {
          mpMixture[i].weight = log(mpMixture[i].weight_accum / accum_sum);
        }
        
        if (!mpMixture[i].estimates->mpMacro) 
        {
        //!!! This is just too ugly hack
          weight_accum_den = mpMixture[i].weight_accum_den;
          mpMixture[i].estimates->UpdateFromAccums(pModelSet);
        }
      }
    }
  }
  

  //***************************************************************************
  FLOAT *
  XFormInstance::
  XFormPass(FLOAT *in_vec, int time, PropagDir dir)
  {
    if (time != UNDEF_TIME && this->time == time) return this->out_vec;
  
    this->time = time;
  
    if (this->mpInput) {
      in_vec = this->mpInput->XFormPass(in_vec, time, dir);
    }
  
    mpXForm->Evaluate(in_vec,this->out_vec,memory,dir);
  
  /*  {int i;
    for (i=0; i<xformInst->mOutSize; i++)
      printf("%.2f ", xformInst->out_vec[i]);
    printf("%s\n", dir == FORWARD ? "FORWARD" : "BACKWARD");}*/
  
    return out_vec;
  }; // XFormPass(FLOAT *in_vec, int time, PropagDir dir)
  
  
  
  
  
  
  //**************************************************************************  
  //**************************************************************************  
  // CompositeXForm section
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  CompositeXForm::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDir  direction)
  {
    size_t  i;
    size_t  j;
  
    for (i = 0; i < mNLayers; i++) 
    {
      FLOAT *in  = (i == 0)          ? pInputVector  : layer[i-1].out_vec;
      FLOAT *out = (i == mNLayers-1) ? pOutputVector : layer[i].out_vec;
      
      for (j = 0; j < layer[i].mNBlocks; j++) 
      {
        layer[i].block[j]->Evaluate(in, out, pMemory, direction);
        in  += layer[i].block[j]->mInSize;
        out += layer[i].block[j]->mOutSize;
        pMemory += layer[i].block[j]->mMemorySize;
      }
    }
    
    return pOutputVector;
  }; //Evaluate(...)

  
  
  //**************************************************************************  
  //**************************************************************************  
  // LinearXForm section
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  LinearXForm::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDir  direction)
  {
    size_t c; // column counter
    size_t r; // row counter
    
    // matrix multiplication
    for (r = 0; r < mOutSize; r++) 
    {
      pOutputVector[r] = 0.0;
      for (c = 0; c < mInSize; c++) 
      {
        pOutputVector[r] += pInputVector[c] * matrix[mInSize * r + c];
      }
    }
    return pOutputVector;  
  }; //Evaluate(...)
  

  //**************************************************************************  
  //**************************************************************************  
  // BiasXForm section
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  BiasXForm::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDir  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[i] + mVector[i];
    }
    return pOutputVector;
  }; //Evaluate(...)
  

  //**************************************************************************  
  //**************************************************************************  
  // FuncXForm section
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  FuncXForm::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDir  direction)
  {
    gFuncTable[mFuncId].funcPtr(pInputVector, pOutputVector, mOutSize);
    return pOutputVector;
  }; //Evaluate(...)
  
      
  //**************************************************************************  
  //**************************************************************************  
  // CopyXForm section
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  CopyXForm::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDir  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[indices[i]];
    }
    return pOutputVector;
  }; //Evaluate(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  // StackingXForm section
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  StackingXForm::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDir  direction)
  {
    FLOAT *stack = reinterpret_cast <FLOAT *> (pMemory);
  
    memmove(stack + (direction == BACKWARD ? mInSize : 0),
            stack + (direction ==  FORWARD ? mInSize : 0),
            (mOutSize - mInSize) * sizeof(FLOAT));
    
    memmove(stack + (direction ==  FORWARD ? mOutSize - mInSize : 0),
            pInputVector, mInSize * sizeof(FLOAT));
  
    if (!horiz_stack) 
    {
      memmove(pOutputVector, stack, mOutSize * sizeof(FLOAT));
    } 
    else 
    { // stacking == HORZ_STACK
      size_t t;
      size_t c;
      size_t stack_size = mOutSize / mInSize;
      
      for (t = 0; t < stack_size; t++) 
      {
        for (c = 0; c < mInSize; c++) 
        {
          pOutputVector[c * stack_size + t] = stack[t * mInSize + c];
        }
      }
    }
    
    return pOutputVector;
  }; //Evaluate(...)
    
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
  
      if (hmm->state[i]->mpMacro) 
      {
        fprintf(fp, "~s \"%s\"", hmm->state[i]->mpMacro->mpName);
        PutNLn(fp, binary);
      } 
      else 
      {
        WriteState(fp, binary, hmm->state[i]);
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
  void 
  ModelSet::
  WriteMmf(const std::string & rFileName, const std::string & rOutputDir,
           const std::string & rOutputExt, bool binary)
  //void WriteHMMSet(const char *mmfName, const char *out_mmf_dir,
  //               const char *out_mmf_ext, bool binary, ModelSet *hmm_set)
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
  
      if (lastFileName == NULL || (rFileName.empty() && strcmp(lastFileName, macro->mpFileName))) 
      {
        // New macro file
        lastFileName = macro->mpFileName;
        
        if (fp && fp != stdout) 
          fclose(fp);
  
        if (!strcmp(!rFileName.empty() ? rFileName.c_str() : macro->mpFileName, "-")) 
        {
          fp = stdout;
        }
        else 
        {
          MakeFileName(mmfile, !rFileName.empty() ? rFileName.c_str() : macro->mpFileName, rOutputDir.c_str(), rOutputExt.c_str());
          
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

  
  //**************************************************************************
  void
  ModelSet::
  Init(FlagType flags)
  {
    if (!my_hcreate_r(100, &mHmmHash)            ||
        !my_hcreate_r(100, &mStateHash)          ||
        !my_hcreate_r( 10, &mMixtureHash)        ||
        !my_hcreate_r( 10, &mMeanHash)           ||
        !my_hcreate_r( 10, &mVarianceHash)       ||
        !my_hcreate_r( 10, &mTransitionHash)     ||
        !my_hcreate_r( 10, &mXFormInstanceHash)  ||
        !my_hcreate_r( 10, &mXFormHash)) 
    {
      Error("Insufficient memory");
    }
  
    this->mpXFormInstances      = NULL;
    this->mpInputXForm          = NULL;
    this->mpFirstMacro          = NULL;
    this->mpLastMacro           = NULL;
    this->mInputVectorSize      = -1;
    this->mParamKind            = -1;
    this->mOutPdfKind           = KID_UNSET;
    this->mDurKind              = KID_UNSET;
    this->mNStates              = 0;
    this->mAllocAccums          = flags & MODEL_SET_WITH_ACCUM ? true : false;
    this->mTotalDelay           = 0;
    InitLogMath();
  
    //Reestimation params
    this->mMinOccurances        = 3;
    this->mMinMixWeight         = MIN_WEGIHT;
    this->mpVarFloor            = NULL;
    this->mUpdateMask           = UM_TRANSITION | UM_MEAN | UM_VARIANCE |
                                  UM_WEIGHT | UM_XFSTATS | UM_XFORM;
    this->mpXFormToUpdate       = NULL;
    this->mNumberOfXFormsToUpdate     = 0;
    this->mGaussLvl2ModelReest  = 0;
    this->mMmiUpdate            = 0;
    this->MMI_E                 = 2.0;
    this->MMI_h                 = 2.0;
    this->MMI_tauI              = 100.0;
  
    InitKwdTable();
  } // Init(...);


  //**************************************************************************  
  void
  ModelSet::
  Release()
  {
    size_t i;  
    
    Scan(MTM_REVERSE_PASS | MTM_ALL, NULL, ReleaseItem, NULL);
  
    ReleaseMacroHash(&mHmmHash);
    ReleaseMacroHash(&mStateHash);
    ReleaseMacroHash(&mMixtureHash);
    ReleaseMacroHash(&mMeanHash);
    ReleaseMacroHash(&mVarianceHash);
    ReleaseMacroHash(&mTransitionHash);
    ReleaseMacroHash(&mXFormHash);
    ReleaseMacroHash(&mXFormInstanceHash);
  
    // :KLUDGE: see the original function Release HMMSet
    for (i = 0; i < mNumberOfXFormsToUpdate; i++) 
    {
      free(mpXFormToUpdate[i].mpShellCommand);
    }
  
    free(mpXFormToUpdate);    
  } // Release();
    
  
  //**************************************************************************  
  void
  ModelSet::
  ReadAccums(const std::string & rFileName, 
             float               weight,
             long *              totFrames, 
             FLOAT *             totLogLike, 
             int                 mmiDenominatorAccums)
  {
    istkstream                in;
    FILE *                    fp;
    char                      macro_name[128];
    struct my_hsearch_data *  hash;
    unsigned int              i;
    int                       t = 0;
    int                       c;
    int                       skip_accum = 0;
    long                      occurances;
    ReadAccumUserData         ud;
    Macro *                   macro;
    int                       mtm = MTM_PRESCAN | 
                                    MTM_STATE | 
                                    MTM_MEAN | 
                                    MTM_VARIANCE | 
                                    MTM_TRANSITION;
  
    macro_name[sizeof(macro_name)-1] = '\0';
  
    // open the file
    in.open(rFileName, ios::binary);
    if (!in.good())
    {
      Error("Cannot open input accumulator file: '%s'", rFileName.c_str());
    }
    fp = in.file();
    
    if (fread(totFrames,  sizeof(long),  1, fp) != 1 ||
        fread(totLogLike, sizeof(FLOAT), 1, fp) != 1) 
    {
      Error("Invalid accumulator file: '%s'", rFileName.c_str());
    }
  
    //:KLUDGE:
    // Assignment of float to long...
    *totFrames  *= weight;
    *totLogLike *= weight;
  
    strcpy(ud.mpFileName, rFileName.c_str());
    ud.fp      = fp;
    ud.mpModelSet = this;
    ud.mWeight  = weight;
    ud.mMmi     = mmiDenominatorAccums;
  
    for (;;) 
    {
      if (skip_accum) 
      { // Skip to the begining of the next macro accumulator
        for (;;) 
        {
          while ((c = getc(fp)) != '~' && c != EOF)
            ;
          
          if (c == EOF) 
            break;
            
          if (strchr("hsmuvt", t = c = getc(fp)) &&
            (c = getc(fp)) == ' ' && (c = getc(fp)) == '"')
          {  
            break;
          }
          
          ungetc(c, fp);
        }
        
        if (c == EOF) 
          break;
      } 
      else 
      {
        if ((c = getc(fp)) == EOF) break;
        
        if (c != '~'       || !strchr("hsmuvt", t = getc(fp)) ||
          getc(fp) != ' ' || getc(fp) != '"') 
        {
          Error("Incomatible accumulator file: '%s'", rFileName.c_str());
        }
      }
  
      for (i=0; (c = getc(fp))!=EOF && c!='"' && i<sizeof(macro_name)-1; i++) 
      {
        macro_name[i] = c;
      }
      macro_name[i] = '\0';
  
      hash = t == 'h' ? &mHmmHash :
             t == 's' ? &mStateHash :
             t == 'm' ? &mMixtureHash :
             t == 'u' ? &mMeanHash :
             t == 'v' ? &mVarianceHash :
             t == 't' ? &mTransitionHash : NULL;
  
      assert(hash);
      if ((macro = FindMacro(hash, macro_name)) == NULL) 
      {
        skip_accum = 1;
        continue;
      }
  
      skip_accum = 0;
      if (fread(&occurances, sizeof(occurances), 1, fp) != 1) 
      {
        Error("Invalid accumulator file: '%s'", rFileName.c_str());
      }
  
      if (!mmiDenominatorAccums) macro->mOccurances += occurances;
  
      switch (t) 
      {
        case 'h': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
        case 's': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
        case 'm': macro->mpData->Scan(mtm, NULL, ReadAccum, &ud); break;
        case 'u': ReadAccum(mt_mean, NULL, macro->mpData, &ud);                break;
        case 'v': ReadAccum(mt_variance, NULL, macro->mpData, &ud);            break;
        case 't': ReadAccum(mt_transition, NULL, macro->mpData, &ud);          break;
        default:  assert(0);
      }
    }
    
    in.close();
    free(ud.mpFileName);    
  }; // ReadAccums(...)

  //**************************************************************************  
  void
  ModelSet::
  WriteAccums(const std::string & rFileName, 
              const std::string & rOutputDir,
              long                totFrames, 
              FLOAT               totLogLike)
  {
    FILE *                fp;
    char                  file_name[1024];
    WriteAccumUserData    ud;  
    
    MakeFileName(file_name, rFileName.c_str(), rOutputDir.c_str(), NULL);
  
    if ((fp = fopen(file_name, "wb")) == NULL) 
    {
      Error("Cannot open output file: '%s'", file_name);
    }
  
    if (fwrite(&totFrames,  sizeof(long),  1, fp) != 1 ||
      fwrite(&totLogLike, sizeof(FLOAT), 1, fp) != 1) 
    {
      Error("Cannot write accumulators to file: '%s'", file_name);
    }
  
    ud.fp  = fp;
    ud.mpFileName  = file_name;
  //  ud.mMmi = MMI_denominator_accums;
  
    Scan(MTM_PRESCAN | (MTM_ALL & ~(MTM_XFORM_INSTANCE|MTM_XFORM)),
              NULL, WriteAccum, &ud);
  
    fclose(fp);
  }; // WriteAccums(...)
  
  
  //**************************************************************************  
  void
  ModelSet::
  NormalizeAccums()
  {
    Scan(MTM_ALL & ~(MTM_XFORM_INSTANCE|MTM_XFORM), NULL,
               NormalizeAccum, NULL);
  }; // NormalizeAccums
  
    
  //**************************************************************************  
  void
  ModelSet::
  ResetAccums()
  {
    Scan(MTM_STATE | MTM_MEAN | MTM_VARIANCE | MTM_TRANSITION,
             NULL, ResetAccum, NULL);
  }; // ResetAccums()
  
  
  //**************************************************************************  
  void
  ModelSet::
  DistributeMacroOccurances()
  {
    Macro *   macro;
    size_t    i;
    size_t    j;
    size_t    k;
  
    for (k = 0; k < mHmmHash.nentries; k++) 
    {
      macro = static_cast <Macro *> (mHmmHash.entry[k]->data);
      
      if (macro->mpData->mpMacro != macro) 
      { 
        macro->mpData->mpMacro->mOccurances += macro->mOccurances;
      }
      
      Hmm *hmm = static_cast <Hmm *> (macro->mpData);
  
      for (i = 0; i < hmm->mNStates - 2; i++) 
      {
        State *state = hmm->state[i];
  
        if (state->mpMacro) 
          state->mpMacro->mOccurances += macro->mOccurances;
  
        if (state->mOutPdfKind == KID_DiagC) 
        {
          for (j = 0; j < state->mNumberOfMixtures; j++) 
          {
            Mixture *mixture = state->mpMixture[j].estimates;
  
            if (mixture->mpMacro)
              mixture->mpMacro->mOccurances += macro->mOccurances;
            
            if (mixture->mpMean->mpMacro)
              mixture->mpMean->mpMacro->mOccurances += macro->mOccurances;
            
            if (mixture->mpVariance->mpMacro) 
              mixture->mpVariance->mpMacro->mOccurances += macro->mOccurances;
          }
        }
      }
  
      if (hmm->mpTransition->mpMacro) 
        hmm->mpTransition->mpMacro->mOccurances += macro->mOccurances;
    }    
  }; // DistributeMacroOccurances()
  
  //**************************************************************************  
  void
  ModelSet::
  ComputeGlobalStats(FLOAT *observation, int time)
  {
    GlobalStatsUserData ud = {observation, time};
    Scan(MTM_STATE | MTM_MIXTURE, NULL, GlobalStats, &ud);
  }; // ComputeGlobalStats(...)
  
  
  //**************************************************************************  
  void
  ModelSet::
  UpdateFromAccums(const std::string & rOutputDir)
  {
    Macro * macro;
    size_t  i;
  
    for (i = 0; i < mHmmHash.nentries; i++) 
    {
      macro = (Macro *) mHmmHash.entry[i]->data;
  //  for (macro = hmm_set->hmm_list; macro != NULL; macro = macro->next) {
      if (macro->mpData->mpMacro != macro) 
        continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Model", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Hmm *) macro->mpData)->UpdateFromAccums(this);
      }
    }
  
    for (i = 0; i < mStateHash.nentries; i++) 
    {
      macro = (Macro *) mStateHash.entry[i]->data;
  //  for (macro = hmm_set->state_list; macro != NULL; macro = macro->next) {
      if (macro->mpData->mpMacro != macro) 
        continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("State", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((State *) macro->mpData)->UpdateFromAccums(this, NULL);
      }
    }
  
    for (i = 0; i < mMixtureHash.nentries; i++) 
    {
      macro = (Macro *) mMixtureHash.entry[i]->data;
  //  for (macro = mixture_list; macro != NULL; macro = macro->next) {
      if (macro->mpData->mpMacro != macro)
        continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Mixture", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Mixture *) macro->mpData)->UpdateFromAccums(this);
      }
    }
  
    for (i = 0; i < mMeanHash.nentries; i++) 
    {
      macro = (Macro *) mMeanHash.entry[i]->data;
  //  for (macro = mean_list; macro != NULL; macro = macro->next) {
      if (macro->mpData->mpMacro != macro) continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Mean vector", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Mean *) macro->mpData)->UpdateFromAccums(this);
      }
    }
  
    for (i = 0; i < mVarianceHash.nentries; i++) 
    {
      macro = (Macro *) mVarianceHash.entry[i]->data;
  //  for (macro = variance_list; macro != NULL; macro = macro->next) {
      if (macro->mpData->mpMacro != macro) continue;
      
      if (macro->mpName != "varFloor1") 
      {
        if (macro->mOccurances < mMinOccurances) 
        {
          WARN_FEW_EXAMPLES("Variance vector", macro->mpName, macro->mOccurances);
        } 
        else 
        {
          ((Variance *) macro->mpData)->UpdateFromAccums(this);
        }
      }
    }
  
    for (i = 0; i < mTransitionHash.nentries; i++) 
    {
      macro = (Macro *) mTransitionHash.entry[i]->data;
  //  for (macro = transition_list; macro != NULL; macro = macro->next) {
      if (macro->mpData->mpMacro != macro) continue;
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Transition matrix ", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Transition *) macro->mpData)->UpdateFromAccums(this);
      }
    }
               
  }; // UpdateHMMSetFromAccums(...)
  
  
  //**************************************************************************  
  void
  ModelSet::
  ReadHMMList(const std::string & rFileName, 
              const std::string & rInMmfDir, 
              const std::string & rInMmfExt)
  {
    struct readline_data      rld = {0};
    char *                    lhmm;
    char *                    fhmm;
    char *                    chptr;
    Macro *                   macro;
    Macro *                   macro2;
    int                       line_no = 0;
    char                      mmfile[1024];
    FILE *                    fp;
  
    if ((fp = my_fopen(rFileName.c_str(), "rt", gpHListFilter)) == NULL) 
    {
      Error("Cannot open file: '%s'", rFileName.c_str());
    }
    
    while ((lhmm = fhmm = readline(fp, &rld)) != NULL) 
    {
      line_no++;
      
      if (getHTKstr(lhmm, &chptr)) 
        Error("%s (%s:%d)", chptr, rFileName.c_str(), line_no);
      
      if (*chptr && getHTKstr(fhmm = chptr, &chptr)) 
        Error("%s (%s:%d)", chptr, rFileName.c_str(), line_no);
      
      if ((macro = FindMacro(&mHmmHash, fhmm)) == NULL) 
      {
        mmfile[0] = '\0';
        
        if (!rInMmfDir.empty()) 
          strcat(strcat(mmfile, rInMmfDir.c_str()), "/");
        
        strcat(mmfile, fhmm);
        
        if (!rInMmfExt.empty()) 
          strcat(strcat(mmfile, "."), rInMmfExt.c_str());
        
        ParseMmf(mmfile, fhmm);
  
        if ((macro = FindMacro(&mHmmHash, fhmm)) == NULL)
          Error("Definition of model '%s' not found in file '%s'", fhmm, mmfile);
      }
      
      if (lhmm != fhmm) 
      {
        current_mmf_name = NULL; // Global variable; macro will not be written to any output MMF
        macro2 = this->pAddMacro('h', lhmm);
        assert(macro2 != NULL);
        
        if (macro2->mpData != NULL) 
        {
          if (gHmmsIgnoreMacroRedefinition == 0) 
          {
            Error("Redefinition of HMM %s (%s:%d)", lhmm, rFileName.c_str(), line_no);
          } 
          else 
          {
            Warning("Redefinition of HMM %s (%s:%d) is ignored",
                    lhmm, rFileName.c_str(), line_no);
          }
        } 
        else 
        {
          macro2->mpData = macro->mpData;
        }
      }
    }
    
    if (ferror(fp) || my_fclose(fp))
      Error("Cannot read HMM list file %s", rFileName.c_str());
  }
  
  
  //**************************************************************************  
  Macro *
  ModelSet::
  pAddMacro(const char type, const std::string & rNewName)
  {
    Macro *                  macro;
    struct my_hsearch_data * hash;
    ENTRY                    e;
    ENTRY *                  ep;
  
    switch (type)
    {
      case 'h': hash = &mHmmHash; break;
      case 's': hash = &mStateHash; break;
      case 'm': hash = &mMixtureHash; break;
      case 'u': hash = &mMeanHash; break;
      case 'v': hash = &mVarianceHash; break;
      case 't': hash = &mTransitionHash; break;
      case 'j': hash = &mXFormInstanceHash; break;
      case 'x': hash = &mXFormHash; break;
      default:  hash = NULL; break;
    }
    
    if (hash == NULL) 
    {
      return NULL;
    }
  
    if ((macro = FindMacro(hash, rNewName.c_str())) != NULL) 
    {
      return macro;
    }
    
    if ((macro = (Macro *) malloc(sizeof(Macro))) == NULL ||
      (macro->mpName = strdup(rNewName.c_str())) == NULL              ||
      (macro->mpFileName = NULL, current_mmf_name
        && (macro->mpFileName = strdup(current_mmf_name)) == NULL)) 
    {
      Error("Insufficient memory");
    }
    
    e.key  = macro->mpName;
    e.data = macro;
  
    if (!my_hsearch_r(e, ENTER, &ep, hash)) 
    {
      Error("Insufficient memory");
    }
  
    macro->mpData = NULL;
    macro->mOccurances = 0;
    macro->mType = type;
  //List of all macros is made to be able to save macros in proper order
    macro->nextAll = NULL;
    macro->prevAll = mpLastMacro;
    
    if (!mpFirstMacro /* => !hmm_set->mpLastMacro  */) 
    {
      mpFirstMacro = macro;
    } 
    else 
    {
      mpLastMacro->nextAll = macro;
    }
    
    mpLastMacro = macro;
    
    return macro;
  }; // pAddMacro(...)

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
        keyword, current_mmf_name, current_mmf_line);
      }
      return (Hmm *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_BeginHMM))
      Error("Keyword <BeginHMM> expected (%s:%d)", current_mmf_name, current_mmf_line);
  
    ReadGlobalOptions(fp);
  
    if (mInputVectorSize == -1)
      Error("<VecSize> is not defined yet (%s:%d)", current_mmf_name, current_mmf_line);
  
    if (mDurKind == -1) 
      mDurKind = KID_NullD;
  
    keyword = GetString(fp, 1);
    
    if (!CheckKwd(keyword, KID_NumStates))
      Error("Keyword <NumStates> expected (%s:%d)", current_mmf_name, current_mmf_line);
  
    nstates = GetInt(fp);
  
    if ((ret = (Hmm *) malloc(sizeof(Hmm) + (nstates-3) * sizeof(State *))) == NULL) 
      Error("Insufficient memory");
  
    ret->mNStates = nstates;
  
    for (i=0; i<nstates-2; i++) 
      ret->state[i] = NULL;
  
    for (i=0; i<nstates-2; i++) 
    {
      keyword = GetString(fp, 1);
      
      if (!CheckKwd(keyword, KID_State)) 
        Error("Keyword <State> expected (%s:%d)", current_mmf_name, current_mmf_line);
  
      state_id = GetInt(fp);
  
  //    printf("%d\n", state_id);
      if (state_id < 2 || state_id >= nstates) 
        Error("State number out of the range (%s:%d)", current_mmf_name, current_mmf_line);
  
      if (ret->state[state_id-2] != NULL) 
        Error("Redefinition of state (%s:%d)", current_mmf_name, current_mmf_line);
  
      ret->state[state_id-2] = ReadState(fp, NULL);
  //    printf("\n%d: %x\n", state_id-2, ret->state[state_id-2]);
    }
  
    ret->mpTransition    = ReadTransition(fp, NULL);
    ret->mpMacro         = macro;
  
    if (ret->mpTransition->mNStates != nstates) 
      Error("Invalid transition matrix size (%s:%d)", current_mmf_name, current_mmf_line);
  
    keyword = GetString(fp, 1);
    
    if (!CheckKwd(keyword, KID_EndHMM))
      Error("Keyword <EndHMM> expected (%s:%d)", current_mmf_name, current_mmf_line);
    
    return ret;
  }

  
  //***************************************************************************
  State *
  ModelSet::
  ReadState(FILE * fp, Macro * macro)
  {
    State *ret;
    char *keyword;
    int mixture_id, i, num_mixes = 1;
    FLOAT mixture_weight;
  
  //  puts("ReadState");
    keyword = GetString(fp, 1);
    
    if (!strcmp(keyword, "~s")) 
    {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mStateHash, keyword)) == NULL) 
      {
        Error("Undefined reference to macro ~s %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
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
  
    ret = (State*) malloc(sizeof(State) + (num_mixes-1)*sizeof(ret->mpMixture[0]));
    
    if (ret == NULL) 
      Error("Insufficient memory");
  
    if (mOutPdfKind == KID_PDFObsVec) 
    {
      int range;
  
      if (!CheckKwd(keyword, KID_ObsCoef)) 
      {
        Error("Keyword <ObsCoef> expected (%s:%d)",
              current_mmf_name, current_mmf_line);
      }
      ret->PDF_obs_coef = GetInt(fp) - 1;
      range = mpInputXForm ? mpInputXForm->mOutSize
                                  : mInputVectorSize;
      if (ret->PDF_obs_coef < 0 || ret->PDF_obs_coef >= range) 
      {
        Error("Parameter <ObsCoef> is out of the range 1:%d (%s:%d)",
              range, current_mmf_name, current_mmf_line);
      }
    } 
    else 
    {
      ret->mNumberOfMixtures = num_mixes;
    //  printf("ptr: %x num_mixes: %d\n", ret, num_mixes);
  
      for (i=0; i<num_mixes; i++) 
        ret->mpMixture[i].estimates = NULL;
  
      if (CheckKwd(keyword, KID_Stream)) 
      {
        if (GetInt(fp) != 1) 
        {
          Error("Stream number out of the range (%s:%d)",
                current_mmf_name, current_mmf_line);
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
                  current_mmf_name, current_mmf_line);
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
                current_mmf_name, current_mmf_line);
        }
  
        if (ret->mpMixture[mixture_id-1].estimates != NULL) 
        {
          Error("Redefinition of mixture %d (%s:%d)",
                mixture_id, current_mmf_name, current_mmf_line);
        }
  
        ret->mpMixture[mixture_id-1].estimates = ReadMixture(fp, NULL);
        ret->mpMixture[mixture_id-1].weight = log(mixture_weight);
        ret->mpMixture[mixture_id-1].weight_accum = 0.0;
      }
    }
    
    ret->mOutPdfKind = mOutPdfKind;
    ret->mID = mNStates;
    mNStates++;
    ret->mpMacro = macro;
  //  puts("ReadState exit");
  
    return ret;
  }


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
              keyword, current_mmf_name, current_mmf_line);
      
      return (Mixture *) macro->mpData;
    }
  
    if ((ret = (Mixture *) malloc(sizeof(Mixture))) == NULL)
      Error("Insufficient memory");
  
  
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
      Error("<VecSize> is not defined yet (%s:%d)", current_mmf_name, current_mmf_line);
    } 
    else if (ret->mpMean->mVectorSize != size || ret->mpVariance->mVectorSize != size) 
    {
      Error("Invalid mean or variance vector size (%s:%d)", current_mmf_name, current_mmf_line);
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
  
    ret->mID = mNStates;
    mNStates++;
    ret->mpMacro = macro;
  
    return ret;
  }


  Mean *
  ModelSet::
  ReadMean(FILE *fp, Macro *macro)
  {
    Mean *ret;
    char *keyword;
    int vec_size, i, accum_size = 0;
  
  //  puts("ReadMean");
    keyword = GetString(fp, 1);
    if (!strcmp(keyword, "~u")) {
      keyword = GetString(fp, 1);
      if ((macro = FindMacro(&mMeanHash, keyword)) == NULL) {
        Error("Undefined reference to macro ~u %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
      }
      return  (Mean *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_Mean)) {
      Error("Keyword <Mean> expected (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    vec_size = GetInt(fp);
  
    if (mAllocAccums) accum_size = (vec_size + 1) * 2; // * 2 for MMI accums
  
    ret = (Mean *) malloc(sizeof(Mean) + (vec_size+accum_size-1) * sizeof(FLOAT));
    if (ret == NULL)  Error("Insufficient memory");
  
    ret->mVectorSize = vec_size;
  //  printf("vec_size: %d\n", vec_size);
    for (i=0; i<vec_size; i++) {
      ret->mVector[i] = GetFloat(fp);
    }
  
    ret->mpXFormStatAccum = NULL;
    ret->mNumberOfXFormStatAccums = 0;
    ret->mUpdatableFromStatAccums = true;
    ret->mpMacro = macro;
    return ret;
  }
  
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
        Error("Undefined reference to macro ~v %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
  
      }
      return (Variance *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_Variance)) {
      Error("Keyword <Variance> expected (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    vec_size = GetInt(fp);
  
    if (mAllocAccums) accum_size = (2 * vec_size + 1) * 2; // * 2 for MMI accums
  
    ret = (Variance *) malloc(sizeof(Variance) + (vec_size+accum_size-1) * sizeof(FLOAT));
    if (ret == NULL) Error("Insufficient memory");
  
    ret->mVectorSize = vec_size;
  
    for (i=0; i<vec_size; i++) {
      ret->mVector[i] = 1.0 / GetFloat(fp);
    }
  
    ret->mpXFormStatAccum = NULL;
    ret->mNumberOfXFormStatAccums = 0;
    ret->mUpdatableFromStatAccums = true;
    ret->mpMacro = macro;
    return ret;
  }
  
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
        Error("Undefined reference to macro ~t %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
      }
      return (Transition *) macro->mpData;
    }
  
    if (!CheckKwd(keyword, KID_TransP)) {
      Error("Keyword <TransP> expected (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    nstates = GetInt(fp);
  
    i = mAllocAccums ? 2 * SQR(nstates) + nstates: SQR(nstates);
  
    if ((ret = (Transition *) malloc(sizeof(Transition) + i * sizeof(ret->matrix[0]))) == NULL) {
      Error("Insufficient memory");
    }
  
    ret->mNStates = nstates;
  
    for (i=0; i < SQR(nstates); i++) {
      ret->matrix[i] = GetFloat(fp);
      ret->matrix[i] = (float) ret->matrix[i] != 0.0 ? log(ret->matrix[i]) : LOG_0;
    }
  
    ret->mpMacro = macro;
    return ret;
  }
  
  
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
        Error("Undefined reference to macro ~j %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
      
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
        Error("<MMFIdMask> different than '*' is not supported (%s:%d)", current_mmf_name, current_mmf_line);
  
      keyword = GetString(fp, 1);
    }
  
    if ((i = ReadParmKind(keyword, TRUE)) != -1) {
      if (mParamKind != -1 && mParamKind != i) 
        Error("ParamKind mismatch (%s:%d)", current_mmf_name, current_mmf_line);
  
      keyword = GetString(fp, 1);
    }
  
    if (CheckKwd(keyword, KID_LinXForm)) 
    {
      keyword = GetString(fp, 1);
  //    Error("Keyword <LinXForm> expected (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    if (!CheckKwd(keyword, KID_VecSize))
      Error("Keyword <VecSize> expected (%s:%d)", current_mmf_name, current_mmf_line);
  
    out_vec_size = GetInt(fp);
  
    ret = (XFormInstance *) malloc(sizeof(*ret)+(out_vec_size-1)*sizeof(FLOAT));
    if (ret == NULL) Error("Insufficient memory");
      ret->mpXForm = ReadXForm(fp, NULL);
  
    if (input == NULL && mInputVectorSize == -1) {
      Error("<VecSize> has not been defined yet (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    if (out_vec_size != ret->mpXForm->mOutSize /* * ret->stackSize*/) {
      Error("XFormInstance <VecSize> must equal to XForm "
            "output size (%s:%d)", current_mmf_name, current_mmf_line);
    }
  
    if (input == NULL) {
      if (ret->mpXForm->mInSize != mInputVectorSize) {
        Error("XForm input size must equal to ~o <VecSize> (%s:%d)",
              current_mmf_name, current_mmf_line);
      }
    } else {
      if (ret->mpXForm->mInSize != input->mOutSize) {
        Error("XForm input size must equal to <Input> <VecSize> (%s:%d)",
              current_mmf_name, current_mmf_line);
      }
    }
  
    ret->memory = NULL;
    if (ret->mpXForm->mMemorySize > 0 &&
      ((ret->memory = (char *) calloc(1, ret->mpXForm->mMemorySize)) == NULL)) {
      Error("Insufficient memory");
    }
  
    ret->mOutSize = out_vec_size;
    ret->next = mpXFormInstances;
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
        Error("Undefined reference to macro ~x %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
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
  
    for (i=0; i < sizeof(gFuncTable)/sizeof(*gFuncTable); i++) 
    {
      if (CheckKwd(keyword, gFuncTable[i].KID)) {
        return (XForm *) ReadFuncXForm(fp, macro, i);
      }
    }
  
    Error("Invalid XForm definition (%s:%d)", current_mmf_name, current_mmf_line);
    return NULL;
  }; //ReadXForm(FILE *fp, Macro *macro) 
  
  
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
  
    if ((ret = (CompositeXForm *) malloc(sizeof(CompositeXForm) +
                                (nlayers-1) * sizeof(ret->layer[0]))) == NULL) {
      Error("Insufficient memory");
    }
  
    ret->mMemorySize = 0;
    ret->mDelay      = 0;
    ret->mNLayers    = nlayers;
  
    for (i=0; i<nlayers; i++) ret->layer[i].block = NULL;
  
    for (i=0; i<nlayers; i++) {
      keyword = GetString(fp, 1);
      if (!CheckKwd(keyword, KID_Layer)) {
        if (nlayers > 1) {
          Error("Keyword <Layer> expected (%s:%d)", current_mmf_name, current_mmf_line);
        }
        layer_id = 1;
      } else {
        layer_id = GetInt(fp);
        keyword = GetString(fp, 1);
      }
  
      if (layer_id < 1 || layer_id > nlayers)
        Error("Layer number out of the range (%s:%d)", current_mmf_name, current_mmf_line);
      
      if (ret->layer[layer_id-1].block != NULL)
        Error("Redefinition of layer (%s:%d)", current_mmf_name, current_mmf_line);
  
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
  
      if ((block = (XForm **)  malloc(sizeof(XForm*) * nblocks)) == NULL) 
        Error("Insufficient memory");
  
      ret->layer[layer_id-1].block   = block;
      ret->layer[layer_id-1].mNBlocks = nblocks;
      
      for (j = 0; j < nblocks; j++) 
        block[j] = NULL;
  
      layer_delay = 0;
      
      for (j = 0; j < nblocks; j++) 
      {
        keyword = GetString(fp, 1);
        
        if (!CheckKwd(keyword, KID_Block)) 
        {
          if (nblocks > 1) 
            Error("Keyword <Block> expected (%s:%d)", current_mmf_name, current_mmf_line);
          
          UngetString();
          block_id = 1;
        } 
        else 
        {
          block_id = GetInt(fp);
        }
  
        if (block_id < 1 || block_id > nblocks)
          Error("Block number out of the range (%s:%d)", current_mmf_name, current_mmf_line);
  
        if (block[block_id-1] != NULL)
          Error("Redefinition of block (%s:%d)", current_mmf_name, current_mmf_line);
  
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
  
      for (j=0; j < ret->layer[i].mNBlocks; j++) 
      {
        layer_in_size  += ret->layer[i].block[j]->mInSize;
        layer_out_size += ret->layer[i].block[j]->mOutSize;
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
          Error("Output size of layer %d (%d) is smaller then input size of layer %d (%d) (%s:%d)",
              i, prev_out_size, i+1, layer_in_size, current_mmf_name, current_mmf_line);
        }
  
        if ((ret->layer[i-1].out_vec = (FLOAT *) malloc(prev_out_size * sizeof(FLOAT))) == NULL) 
        {
          Error("Insufficient memory");
        }
      }
  
      prev_out_size = layer_out_size;
    }
  
    ret->mXFormType = XT_COMPOSITE;
    ret->mpMacro = macro;
    return ret;
  }; // ReadCompositeXForm(FILE *fp, Macro *macro) 
  
  
  //***************************************************************************    
  LinearXForm *
  ModelSet::
  ReadLinearXForm(FILE *fp, Macro *macro)
  {
    LinearXForm *ret;
    int in_size, out_size, i;
  
    out_size = GetInt(fp);
    in_size = GetInt(fp);
  
    ret = (LinearXForm *) malloc(sizeof(LinearXForm)+(out_size*in_size-1)*sizeof(ret->matrix[0]));
    if (ret == NULL) Error("Insufficient memory");
  
  
    for (i=0; i < out_size * in_size; i++) {
      ret->matrix[i] = GetFloat(fp);
    }
  
    ret->mOutSize   = out_size;
    ret->mInSize    = in_size;
    ret->mMemorySize = 0;
    ret->mDelay      = 0;
    ret->mXFormType = XT_LINEAR;
    ret->mpMacro      = macro;
    return ret;
  }; // ReadLinearXForm(FILE *fp, Macro *macro)
  

  //***************************************************************************  
  BiasXForm *
  ModelSet::
  ReadBiasXForm(FILE *fp, Macro *macro)
  {
    BiasXForm *ret;
    int size, i;
  
    size = GetInt(fp);
  
    ret = (BiasXForm *) malloc(sizeof(LinearXForm)+(size-1)*sizeof(ret->mVector[0]));
    if (ret == NULL) Error("Insufficient memory");
  
  
    for (i=0; i < size; i++) {
      ret->mVector[i] = GetFloat(fp);
    }
  
    ret->mInSize = ret->mOutSize = size;
    ret->mMemorySize = 0;
    ret->mDelay      = 0;
    ret->mXFormType = XT_BIAS;
    ret->mpMacro      = macro;
    return ret;
  }; // ReadBiasXForm(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  FuncXForm *
  ModelSet::
  ReadFuncXForm(FILE *fp, Macro *macro, int funcId)
  {
    FuncXForm *ret;
    int size;
  
    size = GetInt(fp);
  
    ret = (FuncXForm *) malloc(sizeof(FuncXForm));
    if (ret == NULL) Error("Insufficient memory");
  
    ret->mFuncId  = funcId;
    ret->mInSize = ret->mOutSize = size;
    ret->mMemorySize = 0;
    ret->mDelay      = 0;
    ret->mXFormType = XT_FUNC;
    ret->mpMacro      = macro;
    return ret;
  }; // ReadFuncXForm(FILE *fp, Macro *macro, int funcId)
  
  
  //***************************************************************************
  CopyXForm *
  ModelSet::
  ReadCopyXForm(FILE *fp, Macro *macro)
  {
    CopyXForm *ret;
    int in_size, out_size, i=0, n, from, step, to;
  
    out_size = GetInt(fp);
    in_size = GetInt(fp);
  
    ret = (CopyXForm *) malloc(sizeof(CopyXForm)+(out_size-1)*sizeof(int));
    if (ret == NULL) Error("Insufficient memory");
  
    while (i < out_size) {
      RemoveSpaces(fp);
      if ((n = fscanf(fp, "%d:%d:%d", &from, &step, &to)) < 1) {
        if (ferror(fp)) {
          Error("Cannot read input file %s", current_mmf_name);
        }
        Error("Integral number expected (%s:%d)", current_mmf_name, current_mmf_line);
      }
  
      if (n == 2)      { to = step; step = 1; }
      else if (n == 1) { to = from; step = 1; }
  
      if (to < 1 || to > in_size) {
        Error("Copy index %d out of range (%s:%d)",
              to, current_mmf_name, current_mmf_line);
      }
  
      for (n = 0; n < (to-from)/step + 1; n++, i++) {
        ret->indices[i] = from + n * step - 1;
      }
    }
  
    ret->mOutSize   = out_size;
    ret->mInSize    = in_size;
    ret->mMemorySize = 0;
    ret->mDelay      = 0;
    ret->mXFormType = XT_COPY;
    ret->mpMacro      = macro;
    return ret;
  }; //ReadCopyXForm(FILE *fp, Macro *macro)
  
  
  //***************************************************************************
  StackingXForm *
  ModelSet::
  ReadStackingXForm(FILE *fp, Macro *macro)
  {
    StackingXForm *ret;
    int stack_size, in_size, out_size;
  
    stack_size = GetInt(fp);
    in_size    = GetInt(fp);
    out_size   = stack_size * in_size;
  
    ret = (StackingXForm *) malloc(sizeof(StackingXForm));
    if (ret == NULL) Error("Insufficient memory");
  
    ret->mOutSize   = out_size;
    ret->mInSize    = in_size;
    ret->mMemorySize = out_size * sizeof(FLOAT);
    ret->mDelay      = stack_size - 1;
    ret->horiz_stack= 0;
    ret->mXFormType = XT_STACKING;
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
  //      Error("Unexpected end of file %s", current_mmf_name);
      }
  
      if (CheckKwd(keyword, KID_VecSize)) 
      {
        i = GetInt(fp);
        if (mInputVectorSize != -1 && mInputVectorSize != i) 
        {
          Error("Mismatch in <VecSize> redefinition (%s:%d)",
                current_mmf_name, current_mmf_line);
        }
  
        mInputVectorSize = i;
        ret = 1;
      } 
      else if (CheckKwd(keyword, KID_StreamInfo)) 
      {
        if (GetInt(fp) != 1) 
        {
          Error("Unsupported definition of multistream (%s:%d)", current_mmf_name, current_mmf_line);
        }
  
        i = GetInt(fp);
        if (mInputVectorSize != -1 && mInputVectorSize != i) 
        {
          Error("Mismatch in <VecSize> redefinition (%s:%d)", current_mmf_name, current_mmf_line);
        }
  
        mInputVectorSize = i;
        ret = 1;
      } 
      else if ((i = ReadParmKind(keyword, TRUE)) != -1) 
      {
        if (mParamKind != -1 && mParamKind != i) 
          Error("Mismatch in paramKind redefinition (%s:%d)", current_mmf_name, current_mmf_line);
        
        mParamKind = i;
        ret = 1;
      } 
      else if ((i = ReadOutPDFKind(keyword)) != -1) 
      {
        if (mOutPdfKind != -1 && mOutPdfKind != i) {
          Error("Mismatch in outPDFKind redefinition (%s:%d)", current_mmf_name, current_mmf_line);
        }
  
        if (i != KID_PDFObsVec && i != KID_DiagC) {
          Error("Unsupported option '%s' (%s:%d)", keyword, current_mmf_name, current_mmf_line);
        }
  
        mOutPdfKind = static_cast<KeywordID> (i);
        ret = 1;
      } 
      else if ((i = ReadDurKind(keyword)) != -1) 
      {
        if (mDurKind != -1 && mDurKind != i) {
          Error("Mismatch in durKind redefinition (%s:%d)", current_mmf_name, current_mmf_line);
        }
        if (i != KID_NullD) {
          Error("Unsupported option '%s' (%s:%d)", keyword, current_mmf_name, current_mmf_line);
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
                  current_mmf_name, current_mmf_line);
        }
        
        inputXForm = ReadXFormInstance(fp, macro);
  
        if (macro->mpData != NULL) {
          Warning("Redefinition of <InputXForm> (%s:%d)",
                  current_mmf_name, current_mmf_line);
  
          // Macro is redefined. New item must be checked for compatibility with
          // the old one (vector size) All references to old
          // item must be replaced and old item must be released
  
          ReplaceItemUserData ud;
          ud.mpOldData = macro->mpData;
          ud.mpNewData = inputXForm;
          ud.mType     = 'j';
  
          this->Scan(MTM_XFORM_INSTANCE|MTM_MIXTURE,NULL,ReplaceItem, &ud);
          ud.mpOldData->Scan(MTM_REVERSE_PASS|MTM_ALL, NULL, ReleaseItem, NULL);
  
          for (unsigned int i = 0; i < mXFormInstanceHash.nentries; i++) 
          {
            if (mXFormInstanceHash.entry[i]->data == ud.mpOldData) 
            {
              mXFormInstanceHash.entry[i]->data = ud.mpNewData;
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
  
  
  void
  ModelSet::
  ParseMmf(const std::string & rFileName, char * expectHMM)
  {
    //ReadHMMSet(rName.c_str(), this, expectHMM);
  
  
    istkstream    input_stream;
    FILE *        fp;
    char *        keyword;
    Macro *       macro;
    MacroData *   data;
  
    data = NULL;
    
    current_mmf_line = 1;
    current_mmf_name = rFileName.c_str();//mmFileName;
  
    /*
    if ((fp = fopen(mmFileName, "rb")) == NULL) 
    {
      Error("Cannot open input MMF %s", mmFileName);
    }
    */
    
    // try to open the stream
    input_stream.open(rFileName, ios::in|ios::binary);
    
    if (!input_stream.good())
    {
      Error("Cannot open input MMF %s", rFileName.c_str());
    }
    
    fp = input_stream.file();
    
    
    for (;;) 
    {
      if ((keyword = GetString(fp, 0)) == NULL) 
      {
        if (ferror(fp)) 
        {
          Error("Cannot read input MMF", rFileName.c_str());
        }
        fclose(fp);
        return;
      }
  
      
      if (keyword[0] == '~' && keyword[2] == '\0' ) 
      {
        char type = keyword[1];
  
        if (type == 'o') {
          if (!ReadGlobalOptions(fp)) {
            Error("No global option defined (%s:%d)", rFileName.c_str(), current_mmf_line);
          }
        } else {
          keyword = GetString(fp, 1);
          if ((macro = pAddMacro(type, keyword)) == NULL) {
            Error("Unrecognized macro type ~%c (%s:%d)", type, rFileName.c_str(), current_mmf_line);
          }
  
          if (macro->mpData != NULL) {
            if (gHmmsIgnoreMacroRedefinition == 0) {
              Error("Redefinition of macro ~%c %s (%s:%d)", type, keyword, rFileName.c_str(), current_mmf_line);
            } else {
  //            Warning("Redefinition of macro ~%c %s (%s:%d) is ignored", type, keyword, mmFileName, current_mmf_line);
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
                    type, keyword, rFileName.c_str(), current_mmf_line);
  
            // Macro is redefined. New item must be checked for compatibility with
            // the old one (vector sizes, delays, memory sizes) All references to
            // old item must be replaced and old item must be released
            // !!! How about AllocateAccumulatorsForXFormStats() and ResetAccumsForHMMSet()
  
  
            unsigned int i;
            struct my_hsearch_data *hash = NULL;
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
              free(ud.mpOldData);
              hash = &mMeanHash;
              break;
            case 'v':
              this->Scan(MTM_MIXTURE, NULL,ReplaceItem, &ud);
              free(ud.mpOldData);
              hash = &mVarianceHash;
              break;
            case 't':
              this->Scan(MTM_HMM, NULL,ReplaceItem, &ud);
              free(ud.mpOldData);
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
            
            for (i = 0; i < hash->nentries; i++) {
              if (hash->entry[i]->data == ud.mpOldData) {
                hash->entry[i]->data = ud.mpNewData;
              }
            }
          }
        }
      } else if (CheckKwd(keyword, KID_BeginHMM)) {
        UngetString();
        if (expectHMM == NULL) {
          Error("Macro definition expected (%s:%d)",rFileName.c_str(),current_mmf_line);
        }
        //macro = AddMacroToHMMSet('h', expectHMM, hmm_set);
        macro = pAddMacro('h', expectHMM);
        
        macro->mpData = ReadHMM(fp, macro);
      } else {
        Error("Unexpected keyword %s (%s:%d)",keyword,rFileName.c_str(),current_mmf_line);
      }
    }
  }

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
          PutFlt(fp, binary, exp(state->mpMixture[i].weight));
          PutNLn(fp, binary);
        }
  
        if (state->mpMixture[i].estimates->mpMacro) 
        {
          fprintf(fp, "~m \"%s\"", state->mpMixture[i].estimates->mpMacro->mpName);
          PutNLn(fp, binary);
        } 
        else 
        {
          WriteMixture(fp, binary, state->mpMixture[i].estimates);
        }
      }
    }
  }
  
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
  void
  ModelSet::
  WriteMean(FILE *fp, bool binary, Mean *mean)
  {
    size_t    i;
  
    PutKwd(fp, binary, KID_Mean);
    PutInt(fp, binary, mean->mVectorSize);
    PutNLn(fp, binary);
  
    for (i=0; i < mean->mVectorSize; i++) {
      PutFlt(fp, binary, mean->mVector[i]);
    }
  
    PutNLn(fp, binary);
  } //WriteMean(FILE *fp, bool binary, Mean *mean)

  //***************************************************************************  
  void
  ModelSet::
  WriteVariance(FILE *fp, bool binary, Variance *variance)
  {
    size_t   i;
  
    PutKwd(fp, binary, KID_Variance);
    PutInt(fp, binary, variance->mVectorSize);
    PutNLn(fp, binary);
  
    for (i=0; i < variance->mVectorSize; i++) {
      PutFlt(fp, binary, 1/variance->mVector[i]);
    }
  
    PutNLn(fp, binary);
  }

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
        FLOAT logtp = transition->matrix[i * transition->mNStates + j];
        PutFlt(fp, binary, logtp > LOG_MIN ? exp(logtp) : 0.0);
      }
  
      PutNLn(fp, binary);
    }
  }

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
      for (i = 0; i < cxf->layer[0].mNBlocks; i++) 
      {
        if (cxf->layer[0].block[i]->mXFormType != XT_LINEAR) 
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

/*
void WriteCompositeXForm(FILE *fp, bool binary, ModelSet *hmm_set, CompositeXForm *xform);
void WriteLinearXForm(   FILE *fp, bool binary, ModelSet *hmm_set,    LinearXForm *xform);
void WriteCopyXForm(     FILE *fp, bool binary, ModelSet *hmm_set,      CopyXForm *xform);
void WriteFuncXForm(     FILE *fp, bool binary, ModelSet *hmm_set,      FuncXForm *xform);
void WriteBiasXForm(     FILE *fp, bool binary, ModelSet *hmm_set,      BiasXForm *xform);
void WriteStackingXForm( FILE *fp, bool binary, ModelSet *hmm_set,  StackingXForm *xform);
*/
  //***************************************************************************
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
      for (i = 0; i < xform->layer[0].mNBlocks; i++) 
      {
        if (xform->layer[0].block[i]->mXFormType != XT_LINEAR) 
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
        PutInt(fp, binary, xform->layer[i].mNBlocks);
        PutNLn(fp, binary);
  
        for (j = 0; j < xform->layer[i].mNBlocks; j++) 
        {
          PutInt(fp, binary, xform->layer[i].block[j]->mOutSize);
        }
  
        PutNLn(fp, binary);
      } 
      else if (xform->layer[i].mNBlocks > 1) 
      {
        PutKwd(fp, binary, KID_NumBlocks);
        PutInt(fp, binary, xform->layer[i].mNBlocks);
        PutNLn(fp, binary);
      }
  
      for (j = 0; j < xform->layer[i].mNBlocks; j++) {
        if (isHTKCompatible || xform->layer[i].mNBlocks > 1) {
          PutKwd(fp, binary, KID_Block);
          PutInt(fp, binary, j+1);
          PutNLn(fp, binary);
        }
        if (xform->layer[i].block[j]->mpMacro) {
          fprintf(fp, "~x \"%s\"", xform->layer[i].block[j]->mpMacro->mpName);
          PutNLn(fp, binary);
        } else {
          WriteXForm(fp, binary, xform->layer[i].block[j]);
        }
      }
    }
  } //WriteCompositeXForm(FILE *fp, bool binary, CompositeXForm *xform)

  
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
      PutFlt(fp, binary, xform->mVector[i]);
    }
    PutNLn(fp, binary);
  } //WriteBiasXForm(FILE *fp, bool binary, BiasXForm *xform)


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
        PutFlt(fp, binary, xform->matrix[i * xform->mInSize + j]);
      }
      PutNLn(fp, binary);
    }
  } //WriteLinearXForm(FILE *fp, bool binary, LinearXForm *xform)
  

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
  void
  ModelSet::
  WriteCopyXForm(FILE *fp, bool binary, CopyXForm *xform)
  {
    size_t      i;
    size_t      j;
    int         step = 0;
    int *       ids = xform->indices;
  
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

  void 
  ModelSet::
  AllocateAccumulatorsForXFormStats() 
  {
    mAllMixuresUpdatableFromStatAccums = true;
    Scan(MTM_XFORM_INSTANCE | MTM_MIXTURE, NULL,
         AllocateXFormStatCachesAndAccums, this);
  }


  //***************************************************************************
  void 
  ModelSet::
  WriteXFormStatsAndRunCommands(const std::string & rOutDir, bool binary)
  {
    HMMSetNodeName                nodeNameBuffer;
    WriteStatsForXFormUserData   userData;
    char                          fileName[1024];
    size_t                        i;
    size_t                        j;
    size_t                        k;
    FILE *                        fp;
  
    struct XfStatsHeader 
    {
      #define PRECISION_FLOAT    0
      #define PRECISION_DOUBLE   1
      #define PRECISION_UNKNOWN 15
      char precision;
      #define STATS_MEAN         0
      #define STATS_COV_LOW_TRI  1
      char stats_type;
      short size;
    } header = {sizeof(FLOAT) == sizeof(float)  ? PRECISION_FLOAT  :
                sizeof(FLOAT) == sizeof(double) ? PRECISION_DOUBLE :
                                                  PRECISION_UNKNOWN};
    userData.mBinary = binary;
    for (k = 0; k < mNumberOfXFormsToUpdate; k++) 
    {
      char *    ext;
      char *    shellCommand = mpXFormToUpdate[k].mpShellCommand;
      
      userData.mpXForm = (LinearXForm *) mpXFormToUpdate[k].mpXForm;
      assert(userData.mpXForm->mXFormType == XT_LINEAR);
  
      MakeFileName(fileName, userData.mpXForm->mpMacro->mpName, rOutDir.c_str(), NULL);
      ext = fileName + strlen(fileName);
      
      strcpy(ext, ".xms"); userData.mMeanFile.mpStatsN = strdup(fileName);
      strcpy(ext, ".xmo"); userData.mMeanFile.mpOccupN = strdup(fileName);
      strcpy(ext, ".xcs"); userData.mCovFile.mpStatsN = strdup(fileName);
      strcpy(ext, ".xco"); userData.mCovFile.mpOccupN = strdup(fileName);
  
      if (userData.mMeanFile.mpStatsN == NULL || userData.mMeanFile.mpOccupN == NULL||
        userData.mCovFile.mpStatsN  == NULL || userData.mCovFile.mpOccupN  == NULL) 
      {
        Error("Insufficient memory");
      }
  
      userData.mMeanFile.mpStatsP = fopen(userData.mMeanFile.mpStatsN, binary?"w":"wt");
      if (userData.mMeanFile.mpStatsP == NULL) 
      {
        Error("Cannot open output file %s",
              userData.mMeanFile.mpStatsN);
      }
  
      if (binary) 
      {
        header.stats_type = STATS_MEAN;
        header.size = userData.mpXForm->mInSize;
        
        if (!isBigEndian()) 
          swap2(header.size);
        
        if (fwrite(&header, sizeof(header), 1, userData.mMeanFile.mpStatsP) != 1) 
        {
          Error("Cannot write to file: %s", userData.mMeanFile.mpStatsN);
        }
      }
  
      userData.mMeanFile.mpOccupP = fopen(userData.mMeanFile.mpOccupN, "wt");
      if (userData.mMeanFile.mpOccupP == NULL) 
      {
        Error("Cannot open output file %s", userData.mMeanFile.mpOccupN);
      }
  
      userData.mCovFile.mpStatsP = fopen(userData.mCovFile.mpStatsN, binary?"w":"wt");
      if (userData.mCovFile.mpStatsP == NULL) 
      {
        Error("Cannot open output file %s", userData.mCovFile.mpStatsN);
      }
  
      if (binary) 
      {
        header.stats_type = STATS_COV_LOW_TRI;
        header.size = userData.mpXForm->mInSize;
        if (!isBigEndian()) 
          swap2(header.size);
        
        if (fwrite(&header, sizeof(header), 1, userData.mCovFile.mpStatsP) != 1) 
        {
          Error("Cannot write to file: %s",
                userData.mCovFile.mpStatsN);
        }
      }
  
      userData.mCovFile.mpOccupP = fopen(userData.mCovFile.mpOccupN, "wt");
      if (userData.mCovFile.mpOccupP == NULL) 
      {
        Error("Cannot open output file %s", userData.mCovFile.mpOccupN);
      }
  
      Scan(MTM_MEAN | MTM_VARIANCE, nodeNameBuffer,
           WriteStatsForXForm, &userData);
  
  
      fclose(userData.mMeanFile.mpStatsP); fclose(userData.mMeanFile.mpOccupP);
      fclose(userData.mCovFile.mpStatsP);  fclose(userData.mCovFile.mpOccupP);
      free(userData.mMeanFile.mpStatsN);   free(userData.mMeanFile.mpOccupN);
      free(userData.mCovFile.mpStatsN);    free(userData.mCovFile.mpOccupN);
  
      strcpy(ext, ".xfm");
      if ((fp = fopen(fileName, "wt")) == NULL) 
      {
        Error("Cannot open output file %s", fileName);
      }
  
      for (i=0; i < userData.mpXForm->mOutSize; i++) 
      {
        for (j=0; j < userData.mpXForm->mInSize; j++) 
        {
          fprintf(fp, " "FLOAT_FMT,
                  userData.mpXForm->matrix[i*userData.mpXForm->mInSize+j]);
        }
        fputs("\n", fp);
      }
  
      fclose(fp);
  
      if (shellCommand != NULL && *shellCommand) 
      {
        TraceLog("Executing command: %s", shellCommand);
        system(shellCommand);
      }
    }
  } // WriteXFormStatsAndRunCommands(const std::string & rOutDir, bool binary)
  //***************************************************************************

  //***************************************************************************
  void 
  ModelSet::
  ReadXFormStats(const std::string & rOutDir, bool binary) 
  {
    HMMSetNodeName                nodeNameBuffer;
    WriteStatsForXFormUserData   userData;
    char                          fileName[1024];
    size_t                        i;
    size_t                        k;
    FILE *                        fp;
  
    struct XfStatsHeader 
    {
      char    precision;
      char    stats_type;
      short   size;
    } header;
  
    userData.mBinary = binary;
    
    for (k = 0; k < mNumberOfXFormsToUpdate; k++) 
    {
      char *  ext;
      userData.mpXForm = (LinearXForm *) mpXFormToUpdate[k].mpXForm;
      assert(userData.mpXForm->mXFormType == XT_LINEAR);
  
      MakeFileName(fileName, userData.mpXForm->mpMacro->mpName, rOutDir.c_str(), NULL);
      ext = fileName + strlen(fileName);
      strcpy(ext, ".xms"); userData.mMeanFile.mpStatsN = strdup(fileName);
      strcpy(ext, ".xmo"); userData.mMeanFile.mpOccupN = strdup(fileName);
      strcpy(ext, ".xcs"); userData.mCovFile.mpStatsN = strdup(fileName);
      strcpy(ext, ".xco"); userData.mCovFile.mpOccupN = strdup(fileName);
  
      if (userData.mMeanFile.mpStatsN == NULL || userData.mMeanFile.mpOccupN == NULL||
        userData.mCovFile.mpStatsN  == NULL || userData.mCovFile.mpOccupN  == NULL) 
      {
        Error("Insufficient memory");
      }
  
      userData.mMeanFile.mpStatsP = fopen(userData.mMeanFile.mpStatsN, "r");
      
      if (userData.mMeanFile.mpStatsP == NULL) 
      {
        Error("Cannot open input file '%s'", userData.mMeanFile.mpStatsN);
      }
  
      if (binary) 
      {
        header.precision = -1;
        fread(&header, sizeof(header), 1, userData.mMeanFile.mpStatsP);
        
        if (!isBigEndian()) 
          swap2(header.size);
        
        if (ferror(userData.mMeanFile.mpStatsP))
        {
          Error("Cannot read input file '%s'", userData.mMeanFile.mpStatsN);
        }        
        else if (header.stats_type != STATS_MEAN
              || header.size       != userData.mpXForm->mInSize
              || header.precision  != (sizeof(FLOAT) == sizeof(float)
                                      ? PRECISION_FLOAT : PRECISION_DOUBLE)) 
        {
          Error("Invalid header in file '%s'", userData.mMeanFile.mpStatsN);
        }
      }
  
      userData.mMeanFile.mpOccupP = fopen(userData.mMeanFile.mpOccupN, "r");
      
      if (userData.mMeanFile.mpOccupP == NULL) 
      {
        Error("Cannot open input file '%s'", userData.mMeanFile.mpOccupN);
      }
  
      userData.mCovFile.mpStatsP = fopen(userData.mCovFile.mpStatsN, "r");
      if (userData.mCovFile.mpStatsP == NULL) 
      {
        Error("Cannot open input file '%s'", userData.mCovFile.mpStatsN);
      }
  
      if (binary) 
      {
        header.precision = -1;
        fread(&header, sizeof(header), 1, userData.mCovFile.mpStatsP);
        if (!isBigEndian()) swap2(header.size);
        if (ferror(userData.mCovFile.mpStatsP)) {
          Error("Cannot read input file '%s'", userData.mCovFile.mpStatsN);
        } else if (header.stats_type != STATS_COV_LOW_TRI
              || header.size       != userData.mpXForm->mInSize
              || header.precision  != (sizeof(FLOAT) == sizeof(float)
                                      ? PRECISION_FLOAT : PRECISION_DOUBLE)) {
          Error("Invalid header in file '%s'", userData.mCovFile.mpStatsN);
        }
      }
  
      userData.mCovFile.mpOccupP = fopen(userData.mCovFile.mpOccupN, "r");
      if (userData.mCovFile.mpOccupP == NULL) {
        Error("Cannot open output file '%s'", userData.mCovFile.mpOccupN);
      }
  
      Scan(MTM_MEAN | MTM_VARIANCE, nodeNameBuffer,
           ReadStatsForXForm, &userData);
  
  
      fclose(userData.mMeanFile.mpStatsP); fclose(userData.mMeanFile.mpOccupP);
      fclose(userData.mCovFile.mpStatsP);  fclose(userData.mCovFile.mpOccupP);
      free(userData.mMeanFile.mpStatsN);   free(userData.mMeanFile.mpOccupN);
      free(userData.mCovFile.mpStatsN);    free(userData.mCovFile.mpOccupN);
  
  
      strcpy(ext, ".xfm");
      if ((fp = fopen(fileName, "r")) == NULL) {
        Error("Cannot open input xforn file '%s'", fileName);
      }
  
      int c =0;
      for (i=0; i < userData.mpXForm->mOutSize * userData.mpXForm->mInSize; i++) {
        c |= fscanf(fp, FLOAT_FMT, &userData.mpXForm->matrix[i]) != 1;
      }
      if (ferror(fp)) {
        Error("Cannot read xform file '%s'", fileName);
      } else if (c) {
        Error("Invalid xform file '%s'", fileName);
      }
      fclose(fp);
    }
  }


  void 
  ModelSet::
  //WriteHMMStats(const char *stat_file)
  WriteHMMStats(const std::string & rFileName)
  {
    Macro *     macro;
    size_t      i = 0;
    size_t      j;
    size_t      k;
    FILE *      fp;
  
    if ((fp = fopen(rFileName.c_str(), "wt")) == NULL) 
    {
      Error("Cannot open output file: '%s'", rFileName.c_str());
    }
  
    for (macro = mpFirstMacro; macro != NULL; macro = macro->nextAll) 
    {
      if (macro->mpData->mpMacro != macro) 
        continue;
        
      if (macro->mType != mt_hmm) 
        continue;
        
      Hmm *hmm = (Hmm *) macro->mpData;
  
      fprintf(fp, "%4d%*c\"%s\" %4ld ", ++i,
                  HIGHER_OF(0,13 - strlen(macro->mpName)), ' ',
                  macro->mpName, macro->mOccurances);
  
      for (j = 0; j < hmm->mNStates-2; j++) 
      {
        State *state = hmm->state[j];
        FLOAT stOccP = 0;
        for (k = 0; k < state->mNumberOfMixtures; k++) 
        {
          stOccP += state->mpMixture[k].weight_accum;
        }
        fprintf(fp, " %10.6f", stOccP);
      }
      fputs("\n", fp);
    }
  
    fclose(fp);
  }

      
  //***************************************************************************
  void 
  ModelSet::
  //ReadXFormList(ModelSet *hmm_set, const char *xformListFileName)
  ReadXFormList(const std::string & rFileName)
  {
    char      line[1024];
    FILE *    fp;
    size_t    nlines=0;
  
    if ((fp = fopen(rFileName.c_str(), "rt")) == NULL) {
        Error("ReadXFormList: Cannot open file: '%s'", rFileName.c_str());
    }
  
  
    while (fgets(line, sizeof(line), fp)) 
    {
      char *            xformName;
      char *            makeXFormShellCommand = NULL;
      char              termChar = '\0';
      size_t            i = strlen(line);
      Macro *           macro;
      MakeXFormCommand *mxfc;
  
      nlines++;
      
      if (line[i-1] != '\n' && getc(fp) != EOF) {
        Error("ReadXFormList: Line %d is too long in file: %s",
              nlines, rFileName.c_str());
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
                termChar, nlines, rFileName.c_str());
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
              xformName, nlines, rFileName.c_str());
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
  }; //ReadXFormList(const std::string & rFileName);
  
  
  void 
  ModelSet::
  ResetXFormInstances()
  {
    XFormInstance *inst;
    for (inst = mpXFormInstances; inst != NULL; inst = inst->next) 
    {
      inst->statCacheTime = UNDEF_TIME;
      inst->time          = UNDEF_TIME;
    }
  }

    
  void 
  ModelSet::
  UpdateStacks(FLOAT *obs, int time,  PropagDir dir) 
  {
    XFormInstance *inst;
    
    for (inst = mpXFormInstances; inst != NULL; inst = inst->next) 
    {
      if (inst->mpXForm->mDelay > 0) {
        XFormPass(inst, obs, time, dir);
      }
    }
  }
      
    
  struct my_hsearch_data 
  ModelSet::
  MakeCIPhoneHash()
  {
    unsigned int              i;
    int                       nCIphns = 0;
    struct my_hsearch_data    tmpHash;
    struct my_hsearch_data    retHash;
    ENTRY                     e;
    ENTRY *                   ep;
  
    if (!my_hcreate_r(100, &tmpHash)) 
      Error("Insufficient memory");
  
    // Find CI HMMs and put them into hash
    for (i = 0; i < mHmmHash.nentries; i++) 
    {
      Macro *macro = (Macro *) mHmmHash.entry[i]->data;
      if (strpbrk(macro->mpName, "+-")) continue;
  
      e.key  = macro->mpName;
      e.data = macro->mpData;
      
      if (!my_hsearch_r(e, ENTER, &ep, &tmpHash)) 
        Error("Insufficient memory");
      
      nCIphns++;
    }
  
    // Find CD HMMs and mark corresponding CI HMMs in the hash
    for (i = 0; i < mHmmHash.nentries; i++) 
    {
      char *    ciname;
      char      chr;
      int       cinlen;
      Macro *   macro = (Macro *) mHmmHash.entry[i]->data;
      
      if (!strpbrk(macro->mpName, "+-")) 
        continue;
  
      ciname = strrchr(macro->mpName, '-');
      
      if (ciname == NULL) ciname = macro->mpName;
      else ciname++;
      
      cinlen = strcspn(ciname, "+");
      chr = ciname[cinlen];
      ciname[cinlen] = '\0';
      e.key  = ciname;
      my_hsearch_r(e, FIND, &ep, &tmpHash);
      ciname[cinlen] = chr;
  
      if (ep != NULL && ep->data != 0) 
      {
        ep->data = NULL;
        nCIphns--;
      }
    }
    
    if (!my_hcreate_r(nCIphns, &retHash)) 
      Error("Insufficient memory");
  
    // To the new hash, collect only those CI HMMs from tmpHash not having CD versions
    for (i = 0; i < tmpHash.nentries; i++) 
    {
      if (tmpHash.entry[i]->data != NULL) 
      {
        Hmm *hmm  = (Hmm *) tmpHash.entry[i]->data;
        int isTee = hmm->mpTransition->matrix[hmm->mNStates - 1] > LOG_MIN;
        e.key  = tmpHash.entry[i]->key;
        e.data = (void *) isTee;
        
        if (!my_hsearch_r(e, ENTER, &ep, &retHash)) 
          Error("Insufficient memory");
      }
    }
    
    my_hdestroy_r(&tmpHash, 0);
    return retHash;
  }
  
  
}; //namespace STK  
  
