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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "hmms.h"
#include "common.h"

#define SQR(x) ((x) * (x))

static int hmm_read_binary;
static int current_mmf_line = 1;
static int string_unget = 0;
static const char *current_mmf_name;
int hmms_ignore_macro_redefinition=1;

HMM *ReadHMM(FILE *fp, HMMSet *hmm_set, Macro *macro);
State *ReadState(FILE *fp, HMMSet *hmm_set, Macro *macro);
Mixture *ReadMixture(FILE *fp, HMMSet *hmm_set, Macro *macro);
Mean *ReadMean(FILE *fp, HMMSet *hmm_set, Macro *macro);
Variance *ReadVariance(FILE *fp, HMMSet *hmm_set, Macro *macro);
Transition *ReadTransition(FILE *fp, HMMSet *hmm_set, Macro *macro);
XformInstance *ReadXformInstance(FILE *fp, HMMSet *hmm_set, Macro *macro);
Xform *ReadXform(FILE *fp, HMMSet *hmm_set, Macro *macro);
int ReadGlobalOptions(FILE *fp, HMMSet *hmm_set);
Macro *AddMacroToHMMSet(char type, const char *name, HMMSet *hmm_set);
Macro *FindMacro(struct my_hsearch_data *macro_hash, const char *name);
void ComputeGConst(Mixture *mix);
char *GetString(FILE *fp, int eofNotExpected);
void UngetString(void);
int GetInt(FILE *fp);
FLOAT GetFloat(FILE *fp);
void RemoveSpaces(FILE *fp);
void ReleaseItem(int macro_type, HMMSetNodeName, void *data, void *userData);
void ReplaceItem(int macro_type, HMMSetNodeName, void *data, void *userData);

static int checkKwd(const char *str, KeywordID kwdID);

char *Kwds[KID_MaxKwdID] = {0};

void initKwdTab()
{
  Kwds[KID_BeginHMM   ] = "BeginHMM";    Kwds[KID_Use        ] = "Use";
  Kwds[KID_EndHMM     ] = "EndHMM";      Kwds[KID_NumMixes   ] = "NumMixes";
  Kwds[KID_NumStates  ] = "NumStates";   Kwds[KID_StreamInfo ] = "StreamInfo";
  Kwds[KID_VecSize    ] = "VecSize";     Kwds[KID_NullD      ] = "NullD";
  Kwds[KID_PoissonD   ] = "PoissonD";    Kwds[KID_GammaD     ] = "GammaD";
  Kwds[KID_RelD       ] = "RelD";        Kwds[KID_GenD       ] = "GenD";
  Kwds[KID_DiagC      ] = "DiagC";       Kwds[KID_FullC      ] = "FullC";
  Kwds[KID_XformC     ] = "XformC";      Kwds[KID_State      ] = "State";
  Kwds[KID_TMix       ] = "TMix";        Kwds[KID_Mixture    ] = "Mixture";
  Kwds[KID_Stream     ] = "Stream";      Kwds[KID_SWeights   ] = "SWeights";
  Kwds[KID_Mean       ] = "Mean";        Kwds[KID_Variance   ] = "Variance";
  Kwds[KID_InvCovar   ] = "InvCovar";    Kwds[KID_Xform      ] = "Xform";
  Kwds[KID_GConst     ] = "GConst";      Kwds[KID_Duration   ] = "Duration";
  Kwds[KID_InvDiagC   ] = "InvDiagC";    Kwds[KID_TransP     ] = "TransP";
  Kwds[KID_DProb      ] = "DProb";       Kwds[KID_LLTC       ] = "LLTC";
  Kwds[KID_LLTCovar   ] = "LLTCovar";    Kwds[KID_XformKind  ] = "XformKind";
  Kwds[KID_ParentXform] = "ParentXform"; Kwds[KID_NumXforms  ] = "NumXforms";
  Kwds[KID_XformSet   ] = "XformSet";    Kwds[KID_LinXform   ] = "LinXform";
  Kwds[KID_Offset     ] = "Offset";      Kwds[KID_Bias       ] = "Bias";
  Kwds[KID_BlockInfo  ] = "BlockInfo";   Kwds[KID_Block      ] = "Block";
  Kwds[KID_BaseClass  ] = "BaseClass";   Kwds[KID_Class      ] = "Class";
  Kwds[KID_XformWgtSet] = "XformWgtSet"; Kwds[KID_ClassXform ] = "ClassXform";
  Kwds[KID_MMFIDMask  ] = "MMFIDMask";   Kwds[KID_Parameters ] = "Parameters";
  Kwds[KID_NumClasses ] = "NumClasses";  Kwds[KID_AdaptKind  ] = "AdaptKind";
  Kwds[KID_Prequal    ] = "Prequal";     Kwds[KID_InputXform ] = "InputXform";
  Kwds[KID_RClass     ] = "RClass";      Kwds[KID_RegTree    ] = "RegTree";
  Kwds[KID_Node       ] = "Node";        Kwds[KID_TNode      ] = "TNode";
  Kwds[KID_HMMSetID   ] = "HMMSetID";    Kwds[KID_ParmKind   ] = "ParmKind";

  /* Non-HTK keywords */
  Kwds[KID_FrmExt     ] = "FrmExt";      Kwds[KID_PDFObsVec  ] = "PDFObsVec";
  Kwds[KID_ObsCoef    ] = "ObsCoef";     Kwds[KID_Input      ] = "Input";
  Kwds[KID_NumLayers  ] = "NumLayers";   Kwds[KID_NumBlocks  ] = "NumBlocks";
  Kwds[KID_Layer      ] = "Layer";       Kwds[KID_Copy       ] = "Copy";
  Kwds[KID_Stacking   ] = "Stacking";
  /* Numeric functions - FuncXform*/
  Kwds[KID_Sigmoid    ] = "Sigmoid";     Kwds[KID_Log        ] = "Log";
  Kwds[KID_Exp        ] = "Exp";         Kwds[KID_Sqrt       ] = "Sqrt";
  Kwds[KID_SoftMax    ] = "SoftMax";

}

struct {
  void (*funcPtr)(FLOAT *, FLOAT*, int);
  KeywordID KID;
} FuncTable[] = {
  {sigmoid_vec, KID_Sigmoid}, {log_vec,     KID_Log},
  {exp_vec,     KID_Exp},     {sqrt_vec,    KID_Sqrt},
  {softmax_vec, KID_SoftMax},
};

void InitHMMSet(HMMSet *hmm_set, int reest)
{
  if(!my_hcreate_r(100, &hmm_set->hmm_hash)            ||
     !my_hcreate_r(100, &hmm_set->state_hash)          ||
     !my_hcreate_r( 10, &hmm_set->mixture_hash)        ||
     !my_hcreate_r( 10, &hmm_set->mean_hash)           ||
     !my_hcreate_r( 10, &hmm_set->variance_hash)       ||
     !my_hcreate_r( 10, &hmm_set->transition_hash)     ||
     !my_hcreate_r( 10, &hmm_set->Xform_instance_hash) ||
     !my_hcreate_r( 10, &hmm_set->Xform_hash)) {
     Error("Insufficient memory");
  }

  hmm_set->Xform_instances     = NULL;
  hmm_set->inputXform          = NULL;
  hmm_set->first_macro         = NULL;
  hmm_set->last_macro          = NULL;
  hmm_set->in_vec_size         = -1;
  hmm_set->param_kind          = -1;
  hmm_set->outPDF_kind         = KID_UNSET;
  hmm_set->dur_kind            = KID_UNSET;
  hmm_set->nstates             = 0;
  hmm_set->nmixtures           = 0;
  hmm_set->alloc_accums        = reest;
  hmm_set->totalDelay          = 0;
  InitLogMath();

  //Reestimation params
  hmm_set->minOccurances       = 3;
  hmm_set->minMixWeight        = MIN_WEGIHT;
  hmm_set->varFloor            = NULL;
  hmm_set->updateMask          = UM_TRANSITION | UM_MEAN | UM_VARIANCE |
                                 UM_WEIGHT | UM_XFSTATS | UM_XFORM;
  hmm_set->xformToUpdate       = NULL;
  hmm_set->nxformsToUpdate     = 0;
  hmm_set->gaussLvl2ModelReest = 0;
  hmm_set->MMIUpdate           = 0;
  hmm_set->MMI_E               = 2.0;
  hmm_set->MMI_h               = 2.0;
  hmm_set->MMI_tauI            = 100.0;


  initKwdTab();
}


typedef struct {
  void *old_data;
  void *new_data;
  int  type;
} ReplaceItem_UserData;

void ReadHMMSet(const char *mmFileName, HMMSet *hmm_set, char *expectHMM)
{
  FILE *fp;
  char *keyword;
  Macro *macro;
  void *data;


  current_mmf_line = 1;
  current_mmf_name   = mmFileName;

//  puts("ReadHMMSet");
  if((fp = fopen(mmFileName, "rb")) == NULL) {
    Error("Cannot open input MMF %s", mmFileName);
  }

  for(;;) {
    if((keyword = GetString(fp, 0)) == NULL) {
      if(ferror(fp)) {
        Error("Cannot read input MMF", mmFileName);
      }
      fclose(fp);
      return;
    }

    if(keyword[0] == '~' && keyword[2] == '\0' ) {
      char type = keyword[1];

      if(type == 'o') {
        if(!ReadGlobalOptions(fp, hmm_set)) {
          Error("No global option defined (%s:%d)", mmFileName, current_mmf_line);
        }
      } else {
        keyword = GetString(fp, 1);
        if((macro = AddMacroToHMMSet(type, keyword, hmm_set)) == NULL) {
          Error("Unrecognized macro type ~%c (%s:%d)", type, mmFileName, current_mmf_line);
        }

        if(macro->data != NULL) {
          if(hmms_ignore_macro_redefinition == 0) {
            Error("Redefinition of macro ~%c %s (%s:%d)", type, keyword, mmFileName, current_mmf_line);
          } else {
//            Warning("Redefinition of macro ~%c %s (%s:%d) is ignored", type, keyword, mmFileName, current_mmf_line);
          }
        }

        data = type == 'h' ? (void *) ReadHMM          (fp, hmm_set, macro) :
               type == 's' ? (void *) ReadState        (fp, hmm_set, macro) :
               type == 'm' ? (void *) ReadMixture      (fp, hmm_set, macro) :
               type == 'u' ? (void *) ReadMean         (fp, hmm_set, macro) :
               type == 'v' ? (void *) ReadVariance     (fp, hmm_set, macro) :
               type == 't' ? (void *) ReadTransition   (fp, hmm_set, macro) :
               type == 'j' ? (void *) ReadXformInstance(fp, hmm_set, macro) :
               type == 'x' ? (void *) ReadXform        (fp, hmm_set, macro) :
                                      NULL;

        assert(data != NULL);

        if(macro->data == NULL) {
          macro->data = data;
        } else {
          Warning("Redefinition of macro ~%c %s (%s:%d)",
                  type, keyword, mmFileName, current_mmf_line);

          // Macro is redefined. New item must be checked for compatibility with
          // the old one (vector sizes, delays, memory sizes) All references to
          // old item must be replaced and old item must be released
          // !!! How about AllocateAccumulatorsForXformStats() and ResetAccumsForHMMSet()


          int i;
          struct my_hsearch_data *hash = NULL;
          ReplaceItem_UserData ud;
          ud.old_data = macro->data;
          ud.new_data = data;
          ud.type     = type;

          switch(type) {
           case 'h':
             ScanHMM((HMM*) ud.old_data, mtm_revpass | mtm_all, NULL,ReleaseItem,NULL);
             hash = &hmm_set->hmm_hash;
             break;
           case 's':
             ScanHMMSet(hmm_set, mtm_hmm, NULL,ReplaceItem, &ud);
             ScanState((State*) ud.old_data, mtm_revpass | mtm_all, NULL,ReleaseItem,NULL);
             hash = &hmm_set->state_hash;
             break;
           case 'm':
             ScanHMMSet(hmm_set, mtm_state, NULL,ReplaceItem, &ud);
             ScanMixture((Mixture*) ud.old_data, mtm_revpass | mtm_all, NULL,ReleaseItem,NULL);
             hash = &hmm_set->mixture_hash;
             break;
           case 'u':
             ScanHMMSet(hmm_set, mtm_mixture, NULL,ReplaceItem, &ud);
             free(ud.old_data);
             hash = &hmm_set->mean_hash;
             break;
           case 'v':
             ScanHMMSet(hmm_set, mtm_mixture, NULL,ReplaceItem, &ud);
             free(ud.old_data);
             hash = &hmm_set->variance_hash;
             break;
           case 't':
             ScanHMMSet(hmm_set, mtm_hmm, NULL,ReplaceItem, &ud);
             free(ud.old_data);
             hash = &hmm_set->transition_hash;
             break;
           case 'j':
             ScanHMMSet(hmm_set, mtm_XformInstance | mtm_mixture,NULL,ReplaceItem, &ud);
             ScanXformInstance((XformInstance*) ud.old_data,mtm_revpass|mtm_all,NULL,ReleaseItem,NULL);
             hash = &hmm_set->Xform_instance_hash;
             break;
           case 'x':
             ScanHMMSet(hmm_set, mtm_Xform | mtm_XformInstance,NULL,ReplaceItem, &ud);
             ScanXform((Xform*) ud.old_data, mtm_revpass | mtm_all, NULL,ReleaseItem,NULL);
             hash = &hmm_set->Xform_hash;
             break;
          }
          for(i = 0; i < hash->nentries; i++) {
            if(hash->entry[i]->data == ud.old_data) {
              hash->entry[i]->data = ud.new_data;
            }
          }
        }
      }
    } else if(checkKwd(keyword, KID_BeginHMM)) {
      UngetString();
      if(expectHMM == NULL) {
        Error("Macro definition expected (%s:%d)",mmFileName,current_mmf_line);
      }
      macro = AddMacroToHMMSet('h', expectHMM, hmm_set);
      macro->data = ReadHMM(fp, hmm_set, macro);
    } else {
      Error("Unexpected keyword %s (%s:%d)",keyword,mmFileName,current_mmf_line);
    }
  }
}

void ReplaceItem(int macro_type, HMMSetNodeName nodeName, void *data, void *userData)
{
  ReplaceItem_UserData *ud = (ReplaceItem_UserData *) userData;
  int i, j;

  if(macro_type == 'h') {
    HMM *hmm = (HMM *) data;
    if(ud->type == 's') {
      for(i = 0; i < hmm->nstates-2; i++) {
        if(hmm->state[i] == ud->old_data) {
          hmm->state[i] = (State*) ud->new_data;
        }
      }
    } else if(hmm->transition == ud->old_data) {
      hmm->transition = (Transition*) ud->new_data;
    }
  } else if(macro_type == 's') {
    State *state = (State *) data;
    if(state->outPDF_kind != KID_PDFObsVec) {
      for(i = 0; i < state->num_mixtures; i++) {
        if(state->mixture[i].estimates == ud->old_data) {
          state->mixture[i].estimates = (Mixture*) ud->new_data;
        }
      }
    }
  } else if(macro_type == 'm') {
    Mixture *mixture = (Mixture *) data;
    if(mixture->mean      == ud->old_data) mixture->mean       = (Mean*)          ud->new_data;
    if(mixture->variance  == ud->old_data) mixture->variance   = (Variance*)      ud->new_data;
    if(mixture->inputXform== ud->old_data) mixture->inputXform = (XformInstance*) ud->new_data;
  } else if(macro_type == 'x') {
    CompositeXform *cxf = (CompositeXform *) data;
    if(cxf->xform_type == XT_COMPOSITE) {
      for(i = 0; i < cxf->nlayers; i++) {
        for(j = 0; j < cxf->layer[i].nblocks; j++) {
          if(cxf->layer[i].block[j] == ud->old_data) {
            cxf->layer[i].block[j] = (Xform*) ud->new_data;
          }
        }
      }
    }
  } else if(macro_type == 'j') {
    XformInstance *xformInstance = (XformInstance *) data;
    if(xformInstance->input == ud->old_data) xformInstance->input = (XformInstance*) ud->new_data;
    if(xformInstance->xform == ud->old_data) xformInstance->xform = (Xform*)         ud->new_data;
  }
}


HMM *ReadHMM(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  HMM *ret;
  char *keyword;
  int nstates, i, state_id;

//  puts("ReadHMM");

  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~h")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->hmm_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~h %s (%s:%d)",
      keyword, current_mmf_name, current_mmf_line);
    }
    return (HMM *) macro->data;
  }

  if(!checkKwd(keyword, KID_BeginHMM)) {
     Error("Keyword <BeginHMM> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  ReadGlobalOptions(fp, hmm_set);

  if(hmm_set->in_vec_size == -1) {
    Error("<VecSize> is not defined yet (%s:%d)", current_mmf_name, current_mmf_line);
  }

  if(hmm_set->dur_kind == -1) hmm_set->dur_kind = KID_NullD;

  keyword = GetString(fp, 1);
  if(!checkKwd(keyword, KID_NumStates)) {
    Error("Keyword <NumStates> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  nstates = GetInt(fp);

  if((ret = (HMM *) malloc(sizeof(HMM) + (nstates-3) * sizeof(State *))) == NULL) {
    Error("Insufficient memory");
  }

  ret->nstates = nstates;

  for(i=0; i<nstates-2; i++) ret->state[i] = NULL;

  for(i=0; i<nstates-2; i++) {
    keyword = GetString(fp, 1);
    if(!checkKwd(keyword, KID_State)) {
      Error("Keyword <State> expected (%s:%d)", current_mmf_name, current_mmf_line);
    }

    state_id = GetInt(fp);

//    printf("%d\n", state_id);
    if(state_id < 2 || state_id >= nstates) {
      Error("State number out of the range (%s:%d)", current_mmf_name, current_mmf_line);
    }

    if(ret->state[state_id-2] != NULL) {
      Error("Redefinition of state (%s:%d)", current_mmf_name, current_mmf_line);
    }

    ret->state[state_id-2] = ReadState(fp, hmm_set, NULL);
//    printf("\n%d: %x\n", state_id-2, ret->state[state_id-2]);
  }

  ret->transition    = ReadTransition(fp, hmm_set, NULL);
  ret->macro         = macro;

  if(ret->transition->nstates != nstates) {
    Error("Invalid transition matrix size (%s:%d)", current_mmf_name, current_mmf_line);
  }

  keyword = GetString(fp, 1);
  if(!checkKwd(keyword, KID_EndHMM)) {
    Error("Keyword <EndHMM> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }
  return ret;
}

State *ReadState(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  State *ret;
  char *keyword;
  int mixture_id, i, num_mixes = 1;
  FLOAT mixture_weight;

//  puts("ReadState");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~s")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->state_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~s %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
    }
    return (State *) macro->data;
  }

  if(hmm_set->outPDF_kind == -1) {
    hmm_set->outPDF_kind = checkKwd(keyword, KID_ObsCoef)
                           ? KID_PDFObsVec : KID_DiagC;
  }

  if(hmm_set->outPDF_kind == KID_PDFObsVec) {
    num_mixes = 0;
  } else if(checkKwd(keyword, KID_NumMixes)) {
    num_mixes = GetInt(fp);

    keyword = GetString(fp, 1);
  }

  ret = (State*) malloc(sizeof(State) + (num_mixes-1)*sizeof(ret->mixture[0]));
  if(ret == NULL) Error("Insufficient memory");

  if(hmm_set->outPDF_kind == KID_PDFObsVec) {
    int range;

    if(!checkKwd(keyword, KID_ObsCoef)) {
      Error("Keyword <ObsCoef> expected (%s:%d)",
            current_mmf_name, current_mmf_line);
    }
    ret->PDF_obs_coef = GetInt(fp) - 1;
    range = hmm_set->inputXform ? hmm_set->inputXform->out_size
                                : hmm_set->in_vec_size;
    if(ret->PDF_obs_coef < 0 || ret->PDF_obs_coef >= range) {
      Error("Parameter <ObsCoef> is out of the range 1:%d (%s:%d)",
            range, current_mmf_name, current_mmf_line);
    }
  } else {
    ret->num_mixtures = num_mixes;
  //  printf("ptr: %x num_mixes: %d\n", ret, num_mixes);

    for(i=0; i<num_mixes; i++) ret->mixture[i].estimates = NULL;

    if(checkKwd(keyword, KID_Stream)) {
      if(GetInt(fp) != 1) {
        Error("Stream number out of the range (%s:%d)",
              current_mmf_name, current_mmf_line);
      }
    } else {
      UngetString();
    }

    for(i=0; i<num_mixes; i++) {
      keyword = GetString(fp, 1);
      if(!checkKwd(keyword, KID_Mixture)) {
        if(num_mixes > 1) {
          Error("Keyword <Mixture> expected (%s:%d)",
                current_mmf_name, current_mmf_line);
        }
        UngetString();
        mixture_id = 1;
        mixture_weight = 1.0;
      } else {
        mixture_id = GetInt(fp);
        mixture_weight = GetFloat(fp);
      }

      if(mixture_id < 1 || mixture_id > num_mixes) {
        Error("Mixture number out of the range (%s:%d)",
              current_mmf_name, current_mmf_line);
      }

      if(ret->mixture[mixture_id-1].estimates != NULL) {
        Error("Redefinition of mixture %d (%s:%d)",
              mixture_id, current_mmf_name, current_mmf_line);
      }

      ret->mixture[mixture_id-1].estimates = ReadMixture(fp, hmm_set, NULL);
      ret->mixture[mixture_id-1].weight = log(mixture_weight);
      ret->mixture[mixture_id-1].weight_accum = 0.0;
    }
  }
  ret->outPDF_kind = hmm_set->outPDF_kind;
  ret->state_id = hmm_set->nstates;
  hmm_set->nstates++;
  ret->macro = macro;
//  puts("ReadState exit");
  return ret;
}


Mixture *ReadMixture(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  Mixture *ret;
  char *keyword;
  int size;

//  puts("ReadMixture");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~m")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->mixture_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~m %s (%s:%d)",
            keyword, current_mmf_name, current_mmf_line);
    }
    return (Mixture *) macro->data;
  }

  if((ret = (Mixture *) malloc(sizeof(Mixture))) == NULL) {
    Error("Insufficient memory");
  }


  if(checkKwd(keyword, KID_InputXform)) {
    ret->inputXform = ReadXformInstance(fp, hmm_set, NULL);
  } else {
    ret->inputXform = hmm_set->inputXform;
    UngetString();
  }

  ret->mean = ReadMean(fp, hmm_set, NULL);
  ret->variance = ReadVariance(fp, hmm_set, NULL);

  size = ret->inputXform ? ret->inputXform->out_size : hmm_set->in_vec_size;
  if(size == -1) {
    Error("<VecSize> is not defined yet (%s:%d)", current_mmf_name, current_mmf_line);
  } else if(ret->mean->vec_size != size || ret->variance->vec_size != size) {
    Error("Invalid mean or variance vector size (%s:%d)", current_mmf_name, current_mmf_line);
  }

  if((keyword = GetString(fp, 0)) != NULL && checkKwd(keyword, KID_GConst)) {
    ret->g_const = GetFloat(fp);
  } else {
    ComputeGConst(ret);
    if(keyword != NULL) UngetString();
  }

  ret->mixture_id = hmm_set->nmixtures;
  hmm_set->nmixtures++;
  ret->macro = macro;

  return ret;
}

Mean *ReadMean(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  Mean *ret;
  char *keyword;
  int vec_size, i, accum_size = 0;

//  puts("ReadMean");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~u")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->mean_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~u %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
    }
    return  (Mean *) macro->data;
  }

  if(!checkKwd(keyword, KID_Mean)) {
     Error("Keyword <Mean> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  vec_size = GetInt(fp);

  if(hmm_set->alloc_accums) accum_size = (vec_size + 1) * 2; // * 2 for MMI accums

  ret = (Mean *) malloc(sizeof(Mean) + (vec_size+accum_size-1) * sizeof(FLOAT));
  if(ret == NULL)  Error("Insufficient memory");

  ret->vec_size = vec_size;
//  printf("vec_size: %d\n", vec_size);
  for(i=0; i<vec_size; i++) {
    ret->vector[i] = GetFloat(fp);
  }

  ret->xformStatAccum = NULL;
  ret->nxformStatAccums = 0;
  ret->updatableFromStatAccums = 1;
  ret->macro = macro;
  return ret;
}

Variance *ReadVariance(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  Variance *ret;
  char *keyword;
  int vec_size, i, accum_size = 0;

//  puts("ReadVariance");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~v")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->variance_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~v %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);

    }
    return (Variance *) macro->data;
  }

  if(!checkKwd(keyword, KID_Variance)) {
     Error("Keyword <Variance> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  vec_size = GetInt(fp);

  if(hmm_set->alloc_accums) accum_size = (2 * vec_size + 1) * 2; // * 2 for MMI accums

  ret = (Variance *) malloc(sizeof(Variance) + (vec_size+accum_size-1) * sizeof(FLOAT));
  if(ret == NULL) Error("Insufficient memory");

  ret->vec_size = vec_size;

  for(i=0; i<vec_size; i++) {
    ret->vector[i] = 1.0 / GetFloat(fp);
  }

  ret->xformStatAccum = NULL;
  ret->nxformStatAccums = 0;
  ret->updatableFromStatAccums = 1;
  ret->macro = macro;
  return ret;
}

Transition *ReadTransition(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  Transition *ret;
  char *keyword;
  int nstates, i;

//  puts("ReadTransition");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~t")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->transition_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~t %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
    }
    return (Transition *) macro->data;
  }

  if(!checkKwd(keyword, KID_TransP)) {
     Error("Keyword <TransP> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  nstates = GetInt(fp);

  i = hmm_set->alloc_accums ? 2 * SQR(nstates) + nstates: SQR(nstates);

  if((ret = (Transition *) malloc(sizeof(Transition) + i * sizeof(ret->matrix[0]))) == NULL) {
    Error("Insufficient memory");
  }

  ret->nstates = nstates;

  for(i=0; i < SQR(nstates); i++) {
    ret->matrix[i] = GetFloat(fp);
    ret->matrix[i] = (float) ret->matrix[i] != 0.0 ? log(ret->matrix[i]) : LOG_0;
  }

  ret->macro = macro;
  return ret;
}

CompositeXform *ReadCompositeXform(FILE *fp, HMMSet *hmm_set, Macro *macro);
LinearXform       *ReadLinearXform(FILE *fp, HMMSet *hmm_set, Macro *macro);
CopyXform           *ReadCopyXform(FILE *fp, HMMSet *hmm_set, Macro *macro);
BiasXform           *ReadBiasXform(FILE *fp, HMMSet *hmm_set, Macro *macro);
FuncXform           *ReadFuncXform(FILE *fp, HMMSet *hmm_set, Macro *macro, int funcId);
StackingXform   *ReadStackingXform(FILE *fp, HMMSet *hmm_set, Macro *macro);

Xform *ReadXform(FILE *fp, HMMSet *hmm_set, Macro *macro) {
  char *keyword;
  int i;

//  puts("ReadXform");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~x")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->Xform_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~x %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
    }
    return (Xform *) macro->data;
  }

  if(checkKwd(keyword, KID_Xform)) {
    return (Xform *) ReadLinearXform(fp, hmm_set, macro);
  }

  if(checkKwd(keyword, KID_Bias)) {
    return (Xform *) ReadBiasXform(fp, hmm_set, macro);
  }

  if(checkKwd(keyword, KID_Copy)) {
    return (Xform *) ReadCopyXform(fp, hmm_set, macro);
  }

  if(checkKwd(keyword, KID_Stacking)) {
    return (Xform *) ReadStackingXform(fp, hmm_set, macro);
  }

  if(checkKwd(keyword, KID_NumLayers) ||
     checkKwd(keyword, KID_NumBlocks) ||
     checkKwd(keyword, KID_BlockInfo))
  {
    UngetString();
    return (Xform *) ReadCompositeXform(fp, hmm_set, macro);
  }

  for(i=0; i < sizeof(FuncTable)/sizeof(*FuncTable); i++) {
    if(checkKwd(keyword, FuncTable[i].KID)) {
      return (Xform *) ReadFuncXform(fp, hmm_set, macro, i);
    }
  }

  Error("Invalid Xform definition (%s:%d)", current_mmf_name, current_mmf_line);
  return NULL;
}

CompositeXform *ReadCompositeXform(FILE *fp, HMMSet *hmm_set, Macro *macro) {
  CompositeXform *ret;
  Xform **block;
  char *keyword;
  int i, j, layer_delay, layer_id, nlayers,
      block_id, nblocks, prev_out_size = 0;

  keyword = GetString(fp, 1);
  if(checkKwd(keyword, KID_NumLayers)) {
    nlayers = GetInt(fp);
  } else {
    nlayers = 1;
    UngetString();
  }

  if((ret = (CompositeXform *) malloc(sizeof(CompositeXform) +
                               (nlayers-1) * sizeof(ret->layer[0]))) == NULL) {
    Error("Insufficient memory");
  }

  ret->memorySize = 0;
  ret->delay      = 0;
  ret->nlayers    = nlayers;

  for(i=0; i<nlayers; i++) ret->layer[i].block = NULL;

  for(i=0; i<nlayers; i++) {
    keyword = GetString(fp, 1);
    if(!checkKwd(keyword, KID_Layer)) {
      if(nlayers > 1) {
        Error("Keyword <Layer> expected (%s:%d)", current_mmf_name, current_mmf_line);
      }
      layer_id = 1;
    } else {
      layer_id = GetInt(fp);
      keyword = GetString(fp, 1);
    }

    if(layer_id < 1 || layer_id > nlayers) {
      Error("Layer number out of the range (%s:%d)", current_mmf_name, current_mmf_line);
    }

    if(ret->layer[layer_id-1].block != NULL) {
      Error("Redefinition of layer (%s:%d)", current_mmf_name, current_mmf_line);
    }

    if(checkKwd(keyword, KID_NumBlocks)) {
      nblocks = GetInt(fp);
    }  else if(checkKwd(keyword, KID_BlockInfo)) {
      nblocks = GetInt(fp);
      for(j = 0; j < nblocks; j++) GetInt(fp); //Blocks' output sizes are not needed
    }  else {
      nblocks = 1;
      UngetString();
    }

    if((block = (Xform **)  malloc(sizeof(Xform*) * nblocks)) == NULL) {
      Error("Insufficient memory");
    }

    ret->layer[layer_id-1].block   = block;
    ret->layer[layer_id-1].nblocks = nblocks;
    for(j = 0; j < nblocks; j++) block[j] = NULL;

    layer_delay = 0;
    for(j = 0; j < nblocks; j++) {
      keyword = GetString(fp, 1);
      if(!checkKwd(keyword, KID_Block)) {
        if(nblocks > 1) {
          Error("Keyword <Block> expected (%s:%d)", current_mmf_name, current_mmf_line);
        }
        UngetString();
        block_id = 1;
      } else {
        block_id = GetInt(fp);
      }

      if(block_id < 1 || block_id > nblocks) {
        Error("Block number out of the range (%s:%d)", current_mmf_name, current_mmf_line);
      }

      if(block[block_id-1] != NULL) {
        Error("Redefinition of block (%s:%d)", current_mmf_name, current_mmf_line);
       }

      block[block_id-1] = ReadXform(fp, hmm_set, NULL);
      ret->memorySize += block[block_id-1]->memorySize;
      layer_delay = HIGHER_OF(layer_delay ,block[block_id-1]->delay);
    }
    ret->delay += layer_delay;
  }

  for(i=0; i<nlayers; i++) {
    int layer_in_size  = 0, layer_out_size = 0;

    for(j=0; j < ret->layer[i].nblocks; j++) {
      layer_in_size  += ret->layer[i].block[j]->in_size;
      layer_out_size += ret->layer[i].block[j]->out_size;
    }

    if(i == nlayers-1) ret->out_size = layer_out_size;
    if(i == 0)         ret->in_size  = layer_in_size;
    else {
      if(prev_out_size < layer_in_size) {
        Error("Output size of layer %d (%d) is smaller then input size of layer %d (%d) (%s:%d)",
             i, prev_out_size, i+1, layer_in_size, current_mmf_name, current_mmf_line);
      }

      if((ret->layer[i-1].out_vec = (FLOAT *) malloc(prev_out_size * sizeof(FLOAT))) == NULL) {
        Error("Insufficient memory");
      }
    }

    prev_out_size = layer_out_size;
  }

  ret->xform_type = XT_COMPOSITE;
  ret->macro = macro;
  return ret;
}


LinearXform *ReadLinearXform(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  LinearXform *ret;
  int in_size, out_size, i;

  out_size = GetInt(fp);
  in_size = GetInt(fp);

  ret = (LinearXform *) malloc(sizeof(LinearXform)+(out_size*in_size-1)*sizeof(ret->matrix[0]));
  if(ret == NULL) Error("Insufficient memory");


  for(i=0; i < out_size * in_size; i++) {
    ret->matrix[i] = GetFloat(fp);
  }

  ret->out_size   = out_size;
  ret->in_size    = in_size;
  ret->memorySize = 0;
  ret->delay      = 0;
  ret->xform_type = XT_LINEAR;
  ret->macro      = macro;
  return ret;
}

BiasXform *ReadBiasXform(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  BiasXform *ret;
  int size, i;

  size = GetInt(fp);

  ret = (BiasXform *) malloc(sizeof(LinearXform)+(size-1)*sizeof(ret->vector[0]));
  if(ret == NULL) Error("Insufficient memory");


  for(i=0; i < size; i++) {
    ret->vector[i] = GetFloat(fp);
  }

  ret->in_size = ret->out_size = size;
  ret->memorySize = 0;
  ret->delay      = 0;
  ret->xform_type = XT_BIAS;
  ret->macro      = macro;
  return ret;
}

FuncXform *ReadFuncXform(FILE *fp, HMMSet *hmm_set, Macro *macro, int funcId)
{
  FuncXform *ret;
  int size;

  size = GetInt(fp);

  ret = (FuncXform *) malloc(sizeof(FuncXform));
  if(ret == NULL) Error("Insufficient memory");

  ret->funcId  = funcId;
  ret->in_size = ret->out_size = size;
  ret->memorySize = 0;
  ret->delay      = 0;
  ret->xform_type = XT_FUNC;
  ret->macro      = macro;
  return ret;
}



CopyXform *ReadCopyXform(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  CopyXform *ret;
  int in_size, out_size, i=0, n, from, step, to;

  out_size = GetInt(fp);
  in_size = GetInt(fp);

  ret = (CopyXform *) malloc(sizeof(CopyXform)+(out_size-1)*sizeof(int));
  if(ret == NULL) Error("Insufficient memory");

  while(i < out_size) {
    RemoveSpaces(fp);
    if((n = fscanf(fp, "%d:%d:%d", &from, &step, &to)) < 1) {
      if(ferror(fp)) {
        Error("Cannot read input file %s", current_mmf_name);
      }
      Error("Integral number expected (%s:%d)", current_mmf_name, current_mmf_line);
    }

    if(n == 2)      { to = step; step = 1; }
    else if(n == 1) { to = from; step = 1; }

    if(to < 1 || to > in_size) {
      Error("Copy index %d out of range (%s:%d)",
            to, current_mmf_name, current_mmf_line);
    }

    for(n = 0; n < (to-from)/step + 1; n++, i++) {
      ret->indices[i] = from + n * step - 1;
    }
  }

  ret->out_size   = out_size;
  ret->in_size    = in_size;
  ret->memorySize = 0;
  ret->delay      = 0;
  ret->xform_type = XT_COPY;
  ret->macro      = macro;
  return ret;
}

StackingXform *ReadStackingXform(FILE *fp, HMMSet *hmm_set, Macro *macro)
{
  StackingXform *ret;
  int stack_size, in_size, out_size;

  stack_size = GetInt(fp);
  in_size    = GetInt(fp);
  out_size   = stack_size * in_size;

  ret = (StackingXform *) malloc(sizeof(StackingXform));
  if(ret == NULL) Error("Insufficient memory");

  ret->out_size   = out_size;
  ret->in_size    = in_size;
  ret->memorySize = out_size * sizeof(FLOAT);
  ret->delay      = stack_size - 1;
  ret->horiz_stack= 0;
  ret->xform_type = XT_STACKING;
  ret->macro      = macro;
  return ret;
}

int IsXformIn1stLayer(Xform *xform, Xform *topXform)
{
  int i;
  if(topXform == NULL)  return 0;
  if(topXform == xform) return 1;

  if(topXform->xform_type == XT_COMPOSITE) {
    CompositeXform *cxf = (CompositeXform *) topXform;
    for(i=0; i<cxf->layer[0].nblocks; i++) {
      if(IsXformIn1stLayer(xform, cxf->layer[0].block[i])) return 1;
    }
  }
  return 0;
}

int Is1Layer1BlockLinearXform(Xform *xform)
{
  CompositeXform *cxf = (CompositeXform *) xform;
  if(cxf == NULL)                                   return 0;
  if(cxf->xform_type == XT_LINEAR)                  return 1;
  if(cxf->xform_type != XT_COMPOSITE)               return 0;
  if(cxf->nlayers > 1 || cxf->layer[0].nblocks > 1) return 0;
  return Is1Layer1BlockLinearXform(cxf->layer[0].block[0]);
}

XformInstance *ReadXformInstance(FILE *fp, HMMSet *hmm_set, Macro *macro) {
  XformInstance *ret, *input = NULL;
  char *keyword;
  int out_vec_size = -1, i;

//  puts("ReadXformInstance");
  keyword = GetString(fp, 1);
  if(!strcmp(keyword, "~j")) {
    keyword = GetString(fp, 1);
    if((macro = FindMacro(&hmm_set->Xform_instance_hash, keyword)) == NULL) {
      Error("Undefined reference to macro ~j %s (%s:%d)", keyword, current_mmf_name, current_mmf_line);
    }
    return (XformInstance *) macro->data;
  }

  if(checkKwd(keyword, KID_Input)) {
    input = ReadXformInstance(fp, hmm_set, NULL);
    keyword = GetString(fp, 1);
  }

  if(checkKwd(keyword, KID_MMFIDMask)) {
    keyword = GetString(fp, 1);
    if(strcmp(keyword, "*")) {
      Error("<MMFIdMask> different than '*' is not supported (%s:%d)", current_mmf_name, current_mmf_line);
    }

    keyword = GetString(fp, 1);
  }

  if((i = ReadParmKind(keyword, TRUE)) != -1) {
    if(hmm_set->param_kind != -1 && hmm_set->param_kind != i) {
      Error("ParamKind mismatch (%s:%d)", current_mmf_name, current_mmf_line);
    }

    keyword = GetString(fp, 1);
  }

  if(checkKwd(keyword, KID_LinXform)) {
    keyword = GetString(fp, 1);
//    Error("Keyword <LinXform> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  if(!checkKwd(keyword, KID_VecSize)) {
      Error("Keyword <VecSize> expected (%s:%d)", current_mmf_name, current_mmf_line);
  }

  out_vec_size = GetInt(fp);

  ret = (XformInstance *) malloc(sizeof(*ret)+(out_vec_size-1)*sizeof(FLOAT));
  if(ret == NULL) Error("Insufficient memory");
    ret->xform = ReadXform(fp, hmm_set, NULL);

  if(input == NULL && hmm_set->in_vec_size == -1) {
    Error("<VecSize> has not been defined yet (%s:%d)", current_mmf_name, current_mmf_line);
  }

  if(out_vec_size != ret->xform->out_size /* * ret->stackSize*/) {
    Error("XformInstance <VecSize> must equal to Xform "
          "output size (%s:%d)", current_mmf_name, current_mmf_line);
  }

  if(input == NULL) {
    if(ret->xform->in_size != hmm_set->in_vec_size) {
      Error("Xform input size must equal to ~o <VecSize> (%s:%d)",
            current_mmf_name, current_mmf_line);
    }
  } else {
    if(ret->xform->in_size != input->out_size) {
      Error("Xform input size must equal to <Input> <VecSize> (%s:%d)",
            current_mmf_name, current_mmf_line);
    }
  }

  ret->memory = NULL;
  if(ret->xform->memorySize > 0 &&
    ((ret->memory = (char *) calloc(1, ret->xform->memorySize)) == NULL)) {
    Error("Insufficient memory");
  }

  ret->out_size = out_vec_size;
  ret->next = hmm_set->Xform_instances;
  hmm_set->Xform_instances = ret;
  ret->input = input;
  ret->macro = macro;
//  puts("ReadXformInstance exit");

  ret->nxformStatCaches = 0;
  ret->xformStatCache   = NULL;
  ret->totalDelay = ret->xform->delay + (input ? input->totalDelay : 0);

  hmm_set->totalDelay = HIGHER_OF(hmm_set->totalDelay, ret->totalDelay);
  return ret;
}


KeywordID ReadOutPDFKind(char *str)
{
  if(     checkKwd(str, KID_DiagC))    return KID_DiagC;
  else if(checkKwd(str, KID_InvDiagC)) return KID_InvDiagC;
  else if(checkKwd(str, KID_FullC))    return KID_FullC;
  else if(checkKwd(str, KID_XformC))   return KID_XformC;
  else if(checkKwd(str, KID_LLTC))     return KID_LLTC;
  else if(checkKwd(str, KID_PDFObsVec))return KID_PDFObsVec;
  else return KID_UNSET;
}

KeywordID ReadDurKind(char *str)
{
  if(     checkKwd(str, KID_NullD))    return KID_NullD;
  else if(checkKwd(str, KID_PoissonD)) return KID_PoissonD;
  else if(checkKwd(str, KID_GammaD))   return KID_GammaD;
  else if(checkKwd(str, KID_GenD))     return KID_GenD;
  else return KID_UNSET;
}


int ReadGlobalOptions(FILE *fp, HMMSet *hmm_set)
{
  int i, ret = 0;
  char *keyword;

//puts("ReadGlobalOptions");
  for(;;) {
    if((keyword = GetString(fp, 0)) == NULL) {
      return ret;
//      Error("Unexpected end of file %s", current_mmf_name);
    }

    if(checkKwd(keyword, KID_VecSize)) {
      i = GetInt(fp);
      if(hmm_set->in_vec_size != -1 && hmm_set->in_vec_size != i) {
        Error("Mismatch in <VecSize> redefinition (%s:%d)",
              current_mmf_name, current_mmf_line);
      }

      hmm_set->in_vec_size = i;
      ret = 1;
    } else if(checkKwd(keyword, KID_StreamInfo)) {
      if(GetInt(fp) != 1) {
        Error("Unsupported definition of multistream (%s:%d)", current_mmf_name, current_mmf_line);
      }

      i = GetInt(fp);
      if(hmm_set->in_vec_size != -1 && hmm_set->in_vec_size != i) {
        Error("Mismatch in <VecSize> redefinition (%s:%d)", current_mmf_name, current_mmf_line);
      }

      hmm_set->in_vec_size = i;
      ret = 1;
    } else if((i = ReadParmKind(keyword, TRUE)) != -1) {
      if(hmm_set->param_kind != -1 && hmm_set->param_kind != i) {
        Error("Mismatch in paramKind redefinition (%s:%d)", current_mmf_name, current_mmf_line);
      }
      hmm_set->param_kind = i;
      ret = 1;
    } else if((i = ReadOutPDFKind(keyword)) != -1) {
      if(hmm_set->outPDF_kind != -1 && hmm_set->outPDF_kind != i) {
        Error("Mismatch in outPDFKind redefinition (%s:%d)", current_mmf_name, current_mmf_line);
      }

      if(i != KID_PDFObsVec && i != KID_DiagC) {
        Error("Unsupported option '%s' (%s:%d)", keyword, current_mmf_name, current_mmf_line);
      }

      hmm_set->outPDF_kind = static_cast<KeywordID> (i);
      ret = 1;
    } else if((i = ReadDurKind(keyword)) != -1) {
      if(hmm_set->dur_kind != -1 && hmm_set->dur_kind != i) {
        Error("Mismatch in durKind redefinition (%s:%d)", current_mmf_name, current_mmf_line);
      }
      if(i != KID_NullD) {
        Error("Unsupported option '%s' (%s:%d)", keyword, current_mmf_name, current_mmf_line);
      }
      hmm_set->dur_kind = static_cast<KeywordID> (i);
      ret = 1;
    } else if(checkKwd(keyword, KID_HMMSetID)) {
      ret = 1;
    } else if(checkKwd(keyword, KID_InputXform)) {
      XformInstance *inputXform;
      Macro *macro = AddMacroToHMMSet(mt_XformInstance, DEFAULT_XFORM_NAME, hmm_set);
      if(macro->data != NULL) {
        if(hmms_ignore_macro_redefinition == 0) {
          Error("Redefinition of <InputXform> (%s:%d)",
                current_mmf_name, current_mmf_line);
        }
      }
      inputXform = ReadXformInstance(fp, hmm_set, macro);

      if(macro->data != NULL) {
        Warning("Redefinition of <InputXform> (%s:%d)",
                current_mmf_name, current_mmf_line);

        // Macro is redefined. New item must be checked for compatibility with
        // the old one (vector size) All references to old
        // item must be replaced and old item must be released

        ReplaceItem_UserData ud;
        ud.old_data = macro->data;
        ud.new_data = inputXform;
        ud.type     = 'j';

        ScanHMMSet(hmm_set, mtm_XformInstance|mtm_mixture,NULL,ReplaceItem, &ud);
        ScanXformInstance(static_cast<XformInstance*>(ud.old_data),
	                  mtm_revpass|mtm_all,
			  NULL,
			  ReleaseItem,
			  NULL);

        for(i = 0; i < hmm_set->Xform_instance_hash.nentries; i++) {
          if(hmm_set->Xform_instance_hash.entry[i]->data == ud.old_data) {
            hmm_set->Xform_instance_hash.entry[i]->data = ud.new_data;
          }
        }
      } else {
        hmm_set->inputXform = inputXform;
        macro->data = inputXform;
      }
      ret = 1;
//    } else if(checkKwd(keyword, KID_LinXform)) {
//      hmm_set->linXform = ReadXformInstance(fp, hmm_set, NULL);
      ret = 1;
    } else {
      UngetString();
      return ret;
    }
  }
}

Macro *AddMacroToHMMSet(char type, const char *name, HMMSet *hmm_set)
{
  Macro *macro;
  struct my_hsearch_data *hash;
  ENTRY e, *ep;

//  puts("AddMacroToHMMSet");
  hash = type == 'h' ? &hmm_set->hmm_hash :
         type == 's' ? &hmm_set->state_hash :
         type == 'm' ? &hmm_set->mixture_hash :
         type == 'u' ? &hmm_set->mean_hash :
         type == 'v' ? &hmm_set->variance_hash :
//          type == 'i' ? &hmm_set->variance_hash :
         type == 't' ? &hmm_set->transition_hash :
         type == 'j' ? &hmm_set->Xform_instance_hash :
         type == 'x' ? &hmm_set->Xform_hash : NULL;


  if(hash == NULL) {
    return NULL;
  }

  if((macro = FindMacro(hash, name)) != NULL) {
    return macro;
  }
  if((macro = (Macro *) malloc(sizeof(Macro))) == NULL ||
     (macro->name = strdup(name)) == NULL              ||
     (macro->file = NULL, current_mmf_name
      && (macro->file = strdup(current_mmf_name)) == NULL)) {
    Error("Insufficient memory");
  }
  e.key  = macro->name;
  e.data = macro;

  if(!my_hsearch_r(e, ENTER, &ep, hash)) {
    Error("Insufficient memory");
  }

  macro->data = NULL;
//  macro->next = *first;
//  *first = macro;
  macro->occurances = 0;
  macro->type = type;
//List of all macros is made to be able to save macros in proper order
  macro->nextAll = NULL;
  macro->prevAll = hmm_set->last_macro;
  if(!hmm_set->first_macro /* => !hmm_set->last_macro  */) {
    hmm_set->first_macro = macro;
  } else {
    hmm_set->last_macro->nextAll = macro;
  }
  hmm_set->last_macro = macro;
//  puts("AddMacroToHMMSet exit");
  return macro;
}

Macro *FindMacro(struct my_hsearch_data *macro_hash, const char *name) {
  ENTRY e, *ep;
  e.key = (char *) name;
  my_hsearch_r(e, FIND, &ep, macro_hash);
  return (Macro *) (ep ? ep->data : NULL);
}

void ReleaseMacroHash(struct my_hsearch_data *macro_hash) {
  int i;
  for(i = 0; i < macro_hash->nentries; i++) {
    Macro *macro = (Macro *) macro_hash->entry[i]->data;
    free(macro->name);
    free(macro->file);
    free(macro);
    macro_hash->entry[i]->data = NULL;
  }
  my_hdestroy_r(macro_hash, 0);
}

void ReleaseItem(int macro_type,HMMSetNodeName nodeName,void *data,void *userData)
{
    int i;
  CompositeXform *cxf = (CompositeXform *) data;
  if(macro_type == 'x' && cxf->xform_type == XT_COMPOSITE) {
    for(i = 0; i < cxf->nlayers-1; i++) free(cxf->layer[i].out_vec);
    for(i = 0; i < cxf->nlayers;   i++) free(cxf->layer[i].block);
  }
  if(macro_type == 'j') free(((XformInstance *) data)->memory);

  free(data);
}

void ReleaseHMMSet(HMMSet *hmm_set)
{
  int i;

  ScanHMMSet(hmm_set, mtm_revpass | mtm_all, NULL, ReleaseItem, NULL);

  ReleaseMacroHash(&hmm_set->hmm_hash);
  ReleaseMacroHash(&hmm_set->state_hash);
  ReleaseMacroHash(&hmm_set->mixture_hash);
  ReleaseMacroHash(&hmm_set->mean_hash);
  ReleaseMacroHash(&hmm_set->variance_hash);
  ReleaseMacroHash(&hmm_set->transition_hash);
  ReleaseMacroHash(&hmm_set->Xform_hash);
  ReleaseMacroHash(&hmm_set->Xform_instance_hash);

  for(i = 0; i < hmm_set->nxformsToUpdate; i++) {
    free(hmm_set->xformToUpdate->shellCommand);
  }

  free(hmm_set->xformToUpdate);
}

void ComputeGConst(Mixture *mix)
{
  FLOAT cov_det = 0;
  int i;

  for(i = 0; i < mix->variance->vec_size; i++) {
    cov_det -= log(mix->variance->vector[i]);
  }
  mix->g_const = cov_det + M_LOG_2PI * mix->variance->vec_size;
}


int checkKwd(const char *str, KeywordID kwdID)
{
  const char *chptr;
  if(str[0] == ':') {
    hmm_read_binary = 1;
    return str[1] == kwdID;
  }

  if(str[0] != '<') return 0;
  for(chptr = Kwds[kwdID], str++; *chptr; chptr++, str++) {
    if(toupper(*chptr) != toupper(*str)) return 0;
  }

  if(str[0] != '>') return 0;

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
  if(string_unget) {
    string_unget = 0;

//    puts(buffer);
    return buffer;
  }

  RemoveSpaces(fp);

  ch = getc(fp);
  if(ch == '\"' || ch == '\'' ) {
    char termChar = ch;

    while((ch = getc(fp)) != EOF && ch != termChar && chptr-buffer < sizeof(buffer)-1) {
      if(ch == '\n') {
        ++lines;
      }
      *chptr++ = ch;
    }

    if(ch == EOF && ferror(fp)) {
      Error("Cannot read input file %s", current_mmf_name);
    }

    if(ch != termChar) {
      Error("Unterminated string constant (%s:%d)", current_mmf_name, current_mmf_line);
    }
    current_mmf_line += lines;
  } else if(ch == '<') {
    *chptr++ = '<';
    while((ch = getc(fp)) != EOF && !isspace(ch) && ch != '>' && chptr-buffer < sizeof(buffer)-1) {
      *chptr++ = ch;
    }

    if(ch == EOF && ferror(fp)) {
      Error("Cannot read input file %s", current_mmf_name);
    }

    if(ch != '>') {
      Error("Unterminated keyword %s (%s:%d)", buffer, current_mmf_name, current_mmf_line);
    }

    *chptr++ = '>';
  } else if(ch == ':') {
    *chptr++ = ':';
    *chptr++ = ch = getc(fp);

    if(ch == EOF){
     if(ferror(fp)) Error("Cannot read input file %s", current_mmf_name);
     else           Error("Unexpected end of file %s", current_mmf_name);
    }
  } else {
    while(ch != EOF && !isspace(ch) && chptr-buffer < sizeof(buffer)-1) {
      *chptr++ = ch;
      ch = getc(fp);
    }

    if(ch != EOF) {
      ungetc(ch, fp);
    } else if(ferror(fp)) {
      Error("Cannot read input file %s", current_mmf_name);
    }

    if(chptr == buffer) {
      if(eofNotExpected) {
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
  
  if(hmm_read_binary) {
    cc = fread(&ret, sizeof(short), 1, fp);
    if(!isBigEndian()) swap2(ret);
  } else {
    RemoveSpaces(fp);
    cc = fscanf(fp, "%hd", &ret);
  }
  
  if(cc != 1) {
    if(ferror(fp)) {
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
  
  if(hmm_read_binary) {
    cc = fread(&ret, sizeof(float), 1, fp);
    if(!isBigEndian()) swap4(ret);  
  } else {
    RemoveSpaces(fp);
    cc = fscanf(fp, "%f", &ret);
  }
   
  if(cc != 1) {
    if(ferror(fp)) {
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
  while(isspace(ch = getc(fp))) {
    if(ch == '\n') {
      ++current_mmf_line;
    }
  }
  if(ch != EOF) {
    ungetc(ch, fp);
  }
}

int qsmacrocmp(const void *a, const void *b) {
  return strcmp(((Macro *)a)->name, ((Macro *)b)->name);
}

/*int MakeHMMLookupTable(HMMSet *hmm_set,
                       char ***sorted_HMM_names, HMM ***lookup_table)
{
  Macro *macros, *macro;
  int i, nHMMs = 0;

  for(macro = hmm_set->hmm_list; macro != NULL; macro = macro->next) nHMMs++;

  macros = (Macro *) malloc(nHMMs * sizeof(Macro));
  *sorted_HMM_names = (char **) malloc(nHMMs * sizeof(char *));
  *lookup_table     =  (HMM **) malloc(nHMMs * sizeof(HMM *));

  if(!macros || !*sorted_HMM_names || !*lookup_table) {
    Error("Insufficient memory");
  }

  for(i=0, macro = hmm_set->hmm_list; macro != NULL; i++, macro = macro->next) {
    macros[i] = *macro;
  }

  qsort(macros, nHMMs, sizeof(Macro), qsmacrocmp);

  for(i=0; i < nHMMs; i++) {
    (*sorted_HMM_names)[i] =     macros[i].name;
    (*lookup_table)[i] = (HMM *) macros[i].data;
  }

  free(macros);

  return nHMMs;
}*/

void WriteHMM(FILE *fp, int binary, HMMSet *hmm_set, HMM *hmm);
void WriteState(FILE *fp, int binary, HMMSet *hmm_set, State *state);
void WriteMixture(FILE *fp, int binary, HMMSet *hmm_set, Mixture *mixture);
void WriteMean(FILE *fp, int binary, HMMSet *hmm_set, Mean *mean);
void WriteVariance(FILE *fp, int binary, HMMSet *hmm_set, Variance *variance);
void WriteTransition(FILE *fp, int binary, HMMSet *hmm_set, Transition *transition);
void WriteXform(FILE *fp, int binary, HMMSet *hmm_set, Xform *xform);
void WriteXformInstance(FILE *fp, int binary, HMMSet *hmm_set, XformInstance *xformInstance);
void WriteGlobalOptions(FILE *fp, int binary, HMMSet *hmm_set);

void PutKwd(FILE *fp, int binary, KeywordID kwdID)
{
  if(binary) {
    putc(':', fp);
    putc(kwdID, fp);
  } else {
    putc('<', fp);
    fputs(Kwds[kwdID], fp);
    fputs("> ", fp);
  }
}

void PutInt(FILE *fp, int binary, int i)
{
  if(binary) {
    short b = i;
    if(!isBigEndian()) swap2(b);
    fwrite(&b, sizeof(short), 1, fp);
  } else {
    fprintf(fp, "%d ", i);
  }
}

void PutFlt(FILE *fp, int binary, FLOAT f)
{
  if(binary) {
    float b = f;
    if(!isBigEndian()) swap4(b);
    fwrite(&b, sizeof(float), 1, fp);
  } else {
    fprintf(fp, FLOAT_FMT" ", f);
  }
}

void PutNLn(FILE *fp, int binary)
{
  if(!binary) putc('\n', fp);
}


void WriteHMMSet(const char *mmfName, const char *out_mmf_dir,
                 const char *out_mmf_ext, int binary, HMMSet *hmm_set)
{
  FILE *fp = NULL;
  Macro *macro;
  char mmfile[1024];
  char *lastFileName = NULL;
  int waitingForNonXform = 1;

  for(macro = hmm_set->first_macro; macro != NULL; macro = macro->nextAll) {
    if(macro->file == NULL) continue; // Artificial macro not read from file

    if(lastFileName == NULL || (!mmfName && strcmp(lastFileName, macro->file))) {
      // New macro file
      lastFileName = macro->file;
      if(fp && fp != stdout) fclose(fp);

      if(!strcmp(mmfName ? mmfName : macro->file, "-")) fp = stdout;
      else {
        MakeFileName(mmfile, mmfName ? mmfName : macro->file, out_mmf_dir, out_mmf_ext);
        if((fp  = fopen(mmfile, "wb")) == NULL) {
          Error("Cannot open output MMF %s", mmfile);
        }
      }
      waitingForNonXform = 1;
      WriteGlobalOptions(fp, binary, hmm_set);
    }

    if(macro->data == hmm_set->inputXform &&
       !strcmp(macro->name, DEFAULT_XFORM_NAME)) {
      fputs("~o ", fp);
      PutKwd(fp, binary, KID_InputXform);
      PutNLn(fp, binary);
    } else {
      fprintf(fp, "~%c \"%s\"", macro->type, macro->name);
      PutNLn(fp, binary);
    }

    if(*(Macro **)macro->data != macro) {
      fprintf(fp, " ~%c \"%s\"", macro->type, (*(Macro **)macro->data)->name);
      PutNLn(fp, binary);
    } else {
      switch(macro->type) {
        case 'x': WriteXform        (fp, binary, hmm_set, static_cast <Xform*>         (macro->data)); break;
        case 'j': WriteXformInstance(fp, binary, hmm_set, static_cast <XformInstance*> (macro->data)); break;
        case 'u': WriteMean         (fp, binary, hmm_set, static_cast <Mean*>          (macro->data)); break;
        case 'v': WriteVariance     (fp, binary, hmm_set, static_cast <Variance*>      (macro->data)); break;
        case 't': WriteTransition   (fp, binary, hmm_set, static_cast <Transition*>    (macro->data)); break;
        case 'm': WriteMixture      (fp, binary, hmm_set, static_cast <Mixture*>       (macro->data)); break;
        case 's': WriteState        (fp, binary, hmm_set, static_cast <State*>         (macro->data)); break;
        case 'h': WriteHMM          (fp, binary, hmm_set, static_cast <HMM*>           (macro->data)); break;
      }
    }
  }
  if(fp && fp != stdout) fclose(fp);
}

void WriteGlobalOptions(FILE *fp, int binary, HMMSet *hmm_set)
{
  char parmkindstr[64];

  fputs("~o ", fp);
  PutKwd(fp, binary, KID_VecSize);
  PutInt(fp, binary, hmm_set->in_vec_size);

  if(ParmKind2Str(hmm_set->param_kind, parmkindstr)) {
    fprintf(fp, "<%s> ", parmkindstr);
  }
  if(hmm_set->outPDF_kind != -1) {
    PutKwd(fp, binary, hmm_set->outPDF_kind);
  }
  if(hmm_set->outPDF_kind != -1) {
    PutKwd(fp, binary, hmm_set->dur_kind);
  }
  PutNLn(fp, binary);
}


void WriteHMM(FILE *fp, int binary, HMMSet *hmm_set, HMM *hmm)
{
  int i;

  PutKwd(fp, binary, KID_BeginHMM);
  PutNLn(fp, binary);
  PutKwd(fp, binary, KID_NumStates);
  PutInt(fp, binary, hmm->nstates);
  PutNLn(fp, binary);

  for(i=0; i < hmm->nstates-2; i++) {
    PutKwd(fp, binary, KID_State);
    PutInt(fp, binary, i+2);
    PutNLn(fp, binary);

    if(hmm->state[i]->macro) {
      fprintf(fp, "~s \"%s\"", hmm->state[i]->macro->name);
      PutNLn(fp, binary);
    } else {
      WriteState(fp, binary, hmm_set, hmm->state[i]);
    }
  }
  if(hmm->transition->macro) {
    fprintf(fp, "~t \"%s\"", hmm->transition->macro->name);
    PutNLn(fp, binary);
  } else {
    WriteTransition(fp, binary, hmm_set, hmm->transition);
  }
  PutKwd(fp, binary, KID_EndHMM);
  PutNLn(fp, binary);
}

void WriteState(FILE *fp, int binary, HMMSet *hmm_set, State *state)
{
  int i;

  if(hmm_set->outPDF_kind == KID_PDFObsVec) {
    PutKwd(fp, binary, KID_ObsCoef);
    PutInt(fp, binary, state->PDF_obs_coef);
    PutNLn(fp, binary);
  } else {
    if(state->num_mixtures > 1) {
      PutKwd(fp, binary, KID_NumMixes);
      PutInt(fp, binary, state->num_mixtures);
      PutNLn(fp, binary);
    }

    for(i=0; i < state->num_mixtures; i++) {
      if(state->num_mixtures > 1) {
        PutKwd(fp, binary, KID_Mixture);
        PutInt(fp, binary, i+1);
        PutFlt(fp, binary, exp(state->mixture[i].weight));
        PutNLn(fp, binary);
      }

      if(state->mixture[i].estimates->macro) {
        fprintf(fp, "~m \"%s\"", state->mixture[i].estimates->macro->name);
        PutNLn(fp, binary);
      } else {
        WriteMixture(fp, binary, hmm_set, state->mixture[i].estimates);
      }
    }
  }
}

void WriteMixture(FILE *fp, int binary, HMMSet *hmm_set, Mixture *mixture)
{
  if(mixture->inputXform != hmm_set->inputXform) {
    PutKwd(fp, binary, KID_InputXform);
    if(mixture->inputXform->macro) {
      fprintf(fp, "~j \"%s\"", mixture->inputXform->macro->name);
      PutNLn(fp, binary);
    } else {
      WriteXformInstance(fp, binary, hmm_set, mixture->inputXform);
    }
  }
  if(mixture->mean->macro) {
    fprintf(fp, "~u \"%s\"", mixture->mean->macro->name);
    PutNLn(fp, binary);
  } else {
    WriteMean(fp, binary, hmm_set, mixture->mean);
  }
  if(mixture->variance->macro) {
    fprintf(fp, "~v \"%s\"", mixture->variance->macro->name);
    PutNLn(fp, binary);
  } else {
    WriteVariance(fp, binary, hmm_set, mixture->variance);
  }
  PutKwd(fp, binary, KID_GConst);
  PutFlt(fp, binary, mixture->g_const);
  PutNLn(fp, binary);
}

void WriteMean(FILE *fp, int binary, HMMSet *hmm_set, Mean *mean)
{
  int i;

  PutKwd(fp, binary, KID_Mean);
  PutInt(fp, binary, mean->vec_size);
  PutNLn(fp, binary);

  for(i=0; i < mean->vec_size; i++) {
    PutFlt(fp, binary, mean->vector[i]);
  }

  PutNLn(fp, binary);
}

void WriteVariance(FILE *fp, int binary, HMMSet *hmm_set, Variance *variance)
{
  int i;

  PutKwd(fp, binary, KID_Variance);
  PutInt(fp, binary, variance->vec_size);
  PutNLn(fp, binary);

  for(i=0; i < variance->vec_size; i++) {
    PutFlt(fp, binary, 1/variance->vector[i]);
  }

  PutNLn(fp, binary);
}

void WriteTransition(FILE *fp, int binary, HMMSet *hmm_set, Transition *transition)
{
  int i, j;

  PutKwd(fp, binary, KID_TransP);
  PutInt(fp, binary, transition->nstates);
  PutNLn(fp, binary);

  for(i=0; i < transition->nstates; i++) {
    for(j=0; j < transition->nstates; j++) {
      FLOAT logtp = transition->matrix[i * transition->nstates + j];
      PutFlt(fp, binary, logtp > LOG_MIN ? exp(logtp) : 0.0);
    }

    PutNLn(fp, binary);
  }
}

void WriteXformInstance(FILE *fp, int binary, HMMSet *hmm_set, XformInstance *xformInstance)
{
  int i, isHTKCompatible = 1;
  char parmkindstr[64];

  CompositeXform *cxf = (CompositeXform *) xformInstance->xform;
  if(xformInstance->input != NULL || cxf == NULL || cxf->macro ||
     cxf->xform_type != XT_COMPOSITE || cxf->nlayers != 1) {
    isHTKCompatible = 0;
  } else for(i = 0; i < cxf->layer[0].nblocks; i++) {
    if(cxf->layer[0].block[i]->xform_type != XT_LINEAR) {
      isHTKCompatible = 0;
      break;
    }
  }

  if(xformInstance->input != NULL) {
    PutKwd(fp, binary, KID_Input);
    PutNLn(fp, binary);

    if(xformInstance->input->macro) {
      fprintf(fp, "~j \"%s\"", xformInstance->input->macro->name);
      PutNLn(fp, binary);
    } else {
      WriteXformInstance(fp, binary, hmm_set, xformInstance->input);
    }
  }

  if(isHTKCompatible) {
    PutKwd(fp, binary, KID_MMFIDMask);
    fputs("* ", fp);

    if(ParmKind2Str(hmm_set->param_kind, parmkindstr)) {
      fprintf(fp, "<%s> ", parmkindstr);
    }

    PutKwd(fp, binary, KID_LinXform);
  }

  PutKwd(fp, binary, KID_VecSize);
  PutInt(fp, binary, xformInstance->out_size);
  PutNLn(fp, binary);

  if(xformInstance->xform->macro) {
    fprintf(fp, "~x \"%s\"", xformInstance->xform->macro->name);
    PutNLn(fp, binary);
  } else {
    WriteXform(fp, binary, hmm_set, xformInstance->xform);
  }
}


void WriteCompositeXform(FILE *fp, int binary, HMMSet *hmm_set, CompositeXform *xform);
void WriteLinearXform(   FILE *fp, int binary, HMMSet *hmm_set,    LinearXform *xform);
void WriteCopyXform(     FILE *fp, int binary, HMMSet *hmm_set,      CopyXform *xform);
void WriteFuncXform(     FILE *fp, int binary, HMMSet *hmm_set,      FuncXform *xform);
void WriteBiasXform(     FILE *fp, int binary, HMMSet *hmm_set,      BiasXform *xform);
void WriteStackingXform( FILE *fp, int binary, HMMSet *hmm_set,  StackingXform *xform);

void WriteXform(FILE *fp, int binary, HMMSet *hmm_set, Xform *xform)
{
  typedef FLOAT * (*pWrtFunc)(FILE *, int, HMMSet *, Xform *);
  XformType type = xform->xform_type;
  pWrtFunc wrtFunc = type == XT_LINEAR    ? (pWrtFunc) WriteLinearXform    :
                     type == XT_COPY      ? (pWrtFunc) WriteCopyXform      :
                     type == XT_FUNC      ? (pWrtFunc) WriteFuncXform      :
                     type == XT_BIAS      ? (pWrtFunc) WriteBiasXform      :
                     type == XT_STACKING  ? (pWrtFunc) WriteStackingXform  :
                     type == XT_COMPOSITE ? (pWrtFunc) WriteCompositeXform :
                     NULL;

  assert(wrtFunc);
  wrtFunc(fp, binary, hmm_set, xform);
}


void WriteCompositeXform(FILE *fp, int binary, HMMSet *hmm_set, CompositeXform *xform)
{
  int i, j, isHTKCompatible = 1;

  if(xform->macro || xform->nlayers != 1) {
    isHTKCompatible = 0;
  } else for(i = 0; i < xform->layer[0].nblocks; i++) {
    if(xform->layer[0].block[i]->xform_type != XT_LINEAR) {
      isHTKCompatible = 0;
      break;
    }
  }
  if(xform->nlayers > 1) {
    PutKwd(fp, binary, KID_NumLayers);
    PutInt(fp, binary, xform->nlayers);
    PutNLn(fp, binary);
  }

  for(i=0; i < xform->nlayers; i++) {
    if(xform->nlayers > 1) {
      PutKwd(fp, binary, KID_Layer);
      PutInt(fp, binary, i+1);
      PutNLn(fp, binary);
    }

    if(isHTKCompatible) {
      PutKwd(fp, binary, KID_BlockInfo);
      PutInt(fp, binary, xform->layer[i].nblocks);
      PutNLn(fp, binary);

      for(j = 0; j < xform->layer[i].nblocks; j++) {
        PutInt(fp, binary, xform->layer[i].block[j]->out_size);
      }

      PutNLn(fp, binary);
    } else if(xform->layer[i].nblocks > 1) {
      PutKwd(fp, binary, KID_NumBlocks);
      PutInt(fp, binary, xform->layer[i].nblocks);
      PutNLn(fp, binary);
    }

    for(j = 0; j < xform->layer[i].nblocks; j++) {
      if(isHTKCompatible || xform->layer[i].nblocks > 1) {
        PutKwd(fp, binary, KID_Block);
        PutInt(fp, binary, j+1);
        PutNLn(fp, binary);
      }
      if(xform->layer[i].block[j]->macro) {
        fprintf(fp, "~x \"%s\"", xform->layer[i].block[j]->macro->name);
        PutNLn(fp, binary);
      } else {
        WriteXform(fp, binary, hmm_set, xform->layer[i].block[j]);
      }
    }
  }
}

void WriteFuncXform(FILE *fp, int binary, HMMSet *hmm_set, FuncXform *xform)
{
  PutKwd(fp, binary, FuncTable[xform->funcId].KID);
  PutInt(fp, binary, xform->out_size);
  PutNLn(fp, binary);
}


void WriteBiasXform(FILE *fp, int binary, HMMSet *hmm_set, BiasXform *xform)
{
  int i;

  PutKwd(fp, binary, KID_Bias);
  PutInt(fp, binary, xform->out_size);
  PutNLn(fp, binary);
  for(i=0; i < xform->out_size; i++) {
    PutFlt(fp, binary, xform->vector[i]);
  }
  PutNLn(fp, binary);
}


void WriteLinearXform(FILE *fp, int binary, HMMSet *hmm_set, LinearXform *xform)
{
  int i, j;

  PutKwd(fp, binary, KID_Xform);
  PutInt(fp, binary, xform->out_size);
  PutInt(fp, binary, xform->in_size);
  PutNLn(fp, binary);
  for(i=0; i < xform->out_size; i++) {
    for(j=0; j < xform->in_size; j++) {
      PutFlt(fp, binary, xform->matrix[i * xform->in_size + j]);
    }
    PutNLn(fp, binary);
  }
}

void WriteStackingXform(FILE *fp, int binary, HMMSet *hmm_set, StackingXform *xform)
{
  PutKwd(fp, binary, KID_Stacking);
  PutInt(fp, binary, xform->out_size / xform->in_size);
  PutInt(fp, binary, xform->in_size);
  PutNLn(fp, binary);
}

void WriteCopyXform(FILE *fp, int binary, HMMSet *hmm_set, CopyXform *xform)
{
  int i, j, step = 0;
  int *ids = xform->indices;

  fprintf(fp, "<Copy> %d %d\n", xform->out_size, xform->in_size);

  for(i=0; i < xform->out_size; i++) {
    if(i + 1 < xform->out_size) {
      step = ids[i+1] - ids[i];
    }

    for(j = i + 2; j < xform->out_size && ids[j] - ids[j-1] == step; j++);

    if(step == 1 && j > i + 2) {
      fprintf(fp, " %d:%d",    ids[i]+1, ids[j-1]+1);
      i = j-1;
    } else if(j > i + 3) {
      fprintf(fp, " %d:%d:%d", ids[i]+1, step, ids[j-1]+1);
      i = j-1;
    } else {
      fprintf(fp, " %d", ids[i]+1);
    }
  }
  fputs("\n", fp);
}

void ResetAccum(int macro_type, HMMSetNodeName nodeName,
                void *data, void *userData) {
  int i, size;
  FLOAT *vector = NULL;

  if(macro_type == mt_mean || macro_type == mt_variance) {
    if(macro_type == mt_mean) {
      size   = ((Mean *)data)->vec_size;
      vector = ((Mean *)data)->vector + size;
      size   = (size + 1) * 2;
    } else if(macro_type == mt_variance) {
      size   = ((Variance *)data)->vec_size;
      vector = ((Variance *)data)->vector + size;
      size   = (size * 2 + 1) * 2;
    }

    for(i = 0; i < size; i++) vector[i] = 0;

  } else if(macro_type == mt_state) {
    State *state = (State *) data;
    if(state->outPDF_kind == KID_DiagC) {
      for(i = 0; i < state->num_mixtures; i++) {
        state->mixture[i].weight_accum     = 0;
        state->mixture[i].weight_accum_den = 0;
      }
    }
  } else if(macro_type == mt_transition) {
    size   = SQR(((Transition *) data)->nstates);
    vector = ((Transition *) data)->matrix + size;

    for(i = 0; i < size; i++) vector[i] = LOG_0;
  }
}

void ResetAccumsForHMMSet(HMMSet *hmm_set)
{
  ScanHMMSet(hmm_set, mtm_state | mtm_mean | mtm_variance | mtm_transition,
             NULL, ResetAccum, NULL);
}


void UpdateMeanFromAccums(HMMSet *hmm_set, Mean *mean)
{
  int i;
  FLOAT *vec  = mean->vector;
  FLOAT *acc  = mean->vector + 1 * mean->vec_size;
  FLOAT  nrm  = mean->vector  [2 * mean->vec_size];

  if(hmm_set->updateMask & UM_MEAN) {
    if(mean->nxformStatAccums == 0) {
      for(i = 0; i < mean->vec_size; i++) {
        vec[i] = acc[i] / nrm;
      }
    } else { // Updating from xform statistics
      int r, c, pos = 0;
      if(!mean->updatableFromStatAccums) {
        //Warning("Mean 'xxx' is not updated");
        return;
      }
//    If 'updatableFromStatAccums' is true then 'nxformStatAccums' is equal to 1
      for(i=0; i<mean->nxformStatAccums; i++) {
        LinearXform *xform = (LinearXform *) mean->xformStatAccum[i].xform;
        int in_size = xform->in_size;
        FLOAT *mnv  = mean->xformStatAccum[i].stats;

        assert(xform->xform_type == XT_LINEAR);
        assert(pos + xform->out_size <= mean->vec_size);
        for(r = 0; r < xform->out_size; r++) {
          mean->vector[pos + r] = 0;
          for(c = 0; c < in_size; c++) {
            mean->vector[pos + r] += mnv[c] * xform->matrix[in_size * r + c];
          }
        }
        pos += xform->out_size;
      }
      assert(pos == mean->vec_size);
    }
  }
}

void UpdateVarianceFromAccums(HMMSet *hmm_set, Variance *variance)
{
  int i;

  if(hmm_set->updateMask & UM_VARIANCE) {
    if(variance->nxformStatAccums == 0) {
      FLOAT *vec  = variance->vector;
      FLOAT *vac  = variance->vector + 1 * variance->vec_size;
      FLOAT *mac  = variance->vector + 2 * variance->vec_size;
      FLOAT *nrm  = variance->vector + 3 * variance->vec_size;
      for(i = 0; i < variance->vec_size; i++) {
        if(hmm_set->updateMask & UM_OLDMEANVAR) {
          vec[i] = *nrm / vac[i];
        } else {
          vec[i] = 1 / (vac[i] / *nrm - SQR(mac[i] / *nrm));
        }
        if(hmm_set->varFloor && 
           hmm_set->varFloor->vec_size == variance->vec_size) {
          vec[i] = LOWER_OF(vec[i], hmm_set->
          varFloor->vector[i]);
        }
      }
    } else { // Updating from xform statistics
      int r, c, t, pos = 0;
      if(!variance->updatableFromStatAccums) {
        //Warning("Variance 'xxx' is not updated");
        return;
      }
//    If 'updatableFromStatAccums' is true then 'nxformStatAccums' is equal to 1
      for(i=0; i<variance->nxformStatAccums; i++) {
        LinearXform *xform = (LinearXform *) variance->xformStatAccum[i].xform;
        int in_size = xform->in_size;
        FLOAT *cov  = variance->xformStatAccum[i].stats + in_size;

        assert(xform->xform_type == XT_LINEAR);
        assert(pos + xform->out_size <= variance->vec_size);
        for(r = 0; r < xform->out_size; r++) {
          variance->vector[pos + r] = 0.0;

          for(c = 0; c < in_size; c++) {
            FLOAT aux = 0;
            for(t = 0; t <= c; t++) {
              aux += cov[c * (c+1)/2 + t]    * xform->matrix[in_size * r + t];
            }

            for(; t < in_size; t++) {
              aux += cov[t * (t+1)/2 + c]    * xform->matrix[in_size * r + t];
            }
            variance->vector[pos + r] += aux * xform->matrix[in_size * r + c];
          }
          variance->vector[pos + r] = 1 / variance->vector[pos + r];
        }
        pos += xform->out_size;
      }
      assert(pos == variance->vec_size);
    }
  }
}

void UpdateTransitionFromAccums(HMMSet *hmm_set, Transition *transition)
{  int i, j, nstates = transition->nstates;
  FLOAT *vec  = transition->matrix;
  FLOAT *acc  = transition->matrix + 1 * SQR(transition->nstates);

  if(hmm_set->updateMask & UM_TRANSITION) {
    for(i=0; i < nstates; i++) {
      FLOAT nrm = LOG_0;
      for(j=0; j < nstates; j++) {
        LOG_INC(nrm, acc[i * nstates + j]);
      }
      if(nrm == LOG_0) nrm = 0;
      for(j=0; j < nstates; j++) {
        vec[i * nstates + j] = acc[i * nstates + j] - nrm; // it is in log
      }
    }
  }
}


FLOAT weight_accum_den;

void UpdateMixtureFromAccums(HMMSet *hmm_set, Mixture *mixture)
{
    double Djm;
    int i;

  if(hmm_set->MMIUpdate == 1 || hmm_set->MMIUpdate == -1) {
    int vec_size    = mixture->variance->vec_size;
    FLOAT *mean_vec = mixture->mean->vector;
    FLOAT *var_vec  = mixture->variance->vector;

    FLOAT *vac_num  = var_vec + 1 * vec_size;
    FLOAT *mac_num  = var_vec + 2 * vec_size;
    FLOAT *nrm_num  = var_vec + 3 * vec_size;

    FLOAT *vac_den  = vac_num + 2 * vec_size + 1;
    FLOAT *mac_den  = mac_num + 2 * vec_size + 1;
    FLOAT *nrm_den  = nrm_num + 2 * vec_size + 1;

    if(hmm_set->MMIUpdate == 1) {
      // I-smoothing
      for(i = 0; i < vec_size; i++) {
        mac_num[i] *= (*nrm_num + hmm_set->MMI_tauI) / *nrm_num;
        vac_num[i] *= (*nrm_num + hmm_set->MMI_tauI) / *nrm_num;
      }
      *nrm_num   += hmm_set->MMI_tauI;
    }

    Djm = 0.0;
    // Find minimum Djm leading to positive update of variances
    for(i = 0; i < vec_size; i++) {
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

    Djm = HIGHER_OF(hmm_set->MMI_h * Djm, hmm_set->MMI_E * *nrm_den);

    if(hmm_set->MMIUpdate == -1) {
      // I-smoothing
      for(i = 0; i < vec_size; i++) {
        mac_num[i] *= (*nrm_num + hmm_set->MMI_tauI) / *nrm_num;
        vac_num[i] *= (*nrm_num + hmm_set->MMI_tauI) / *nrm_num;
      }
      *nrm_num   += hmm_set->MMI_tauI;
    }
    
    for(i = 0; i < vec_size; i++) {
      double macn_macd = mac_num[i]-mac_den[i];
      double vacn_vacd = vac_num[i]-vac_den[i];
      double nrmn_nrmd = *nrm_num - *nrm_den;

      double new_mean = (macn_macd + Djm * mean_vec[i]) / (nrmn_nrmd + Djm);
      var_vec[i]     = 1/((vacn_vacd + Djm * (1/var_vec[i] + SQR(mean_vec[i]))) /
                       (nrmn_nrmd + Djm) - SQR(new_mean));
      mean_vec[i]    = new_mean;

      if(hmm_set->varFloor) {
        var_vec[i] = LOWER_OF(var_vec[i], hmm_set->varFloor->vector[i]);
      }
    }
  } else if(hmm_set->MMIUpdate == 2 || hmm_set->MMIUpdate == -2 ) { // MFE update
    int vec_size    = mixture->variance->vec_size;
    FLOAT *mean_vec = mixture->mean->vector;
    FLOAT *var_vec  = mixture->variance->vector;

    FLOAT *vac_mle  = var_vec + 1 * vec_size;
    FLOAT *mac_mle  = var_vec + 2 * vec_size;
    FLOAT *nrm_mle  = var_vec + 3 * vec_size;

    FLOAT *vac_mfe  = vac_mle + 2 * vec_size + 1;
    FLOAT *mac_mfe  = mac_mle + 2 * vec_size + 1;
    FLOAT *nrm_mfe  = nrm_mle + 2 * vec_size + 1;

    if(hmm_set->MMIUpdate == 2) {
      // I-smoothing
      for(i = 0; i < vec_size; i++) {
        mac_mfe[i] += (hmm_set->MMI_tauI / *nrm_mle * mac_mle[i]);
        vac_mfe[i] += (hmm_set->MMI_tauI / *nrm_mle * vac_mle[i]);
      }
      *nrm_mfe += hmm_set->MMI_tauI;
    }
    
    Djm = 0.0;
    // Find minimum Djm leading to positive update of variances
    for(i = 0; i < vec_size; i++) {
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
    Djm = HIGHER_OF(hmm_set->MMI_h * Djm, hmm_set->MMI_E * weight_accum_den);

    if(hmm_set->MMIUpdate == -2) {
      // I-smoothing
      for(i = 0; i < vec_size; i++) {
        mac_mfe[i] += (hmm_set->MMI_tauI / *nrm_mle * mac_mle[i]);
        vac_mfe[i] += (hmm_set->MMI_tauI / *nrm_mle * vac_mle[i]);
      }
      *nrm_mfe += hmm_set->MMI_tauI;
    }

    for(i = 0; i < vec_size; i++) {
      double macn_macd = mac_mfe[i];
      double vacn_vacd = vac_mfe[i];
      double nrmn_nrmd = *nrm_mfe;

      double new_mean = (macn_macd + Djm * mean_vec[i]) / (nrmn_nrmd + Djm);
      var_vec[i]     = 1/((vacn_vacd + Djm * (1/var_vec[i] + SQR(mean_vec[i]))) /
                       (nrmn_nrmd + Djm) - SQR(new_mean));
      mean_vec[i]    = new_mean;

      if(hmm_set->varFloor) {
        var_vec[i] = LOWER_OF(var_vec[i], hmm_set->varFloor->vector[i]);
      }
    }
  } else {
    if(!mixture->variance->macro) {
      UpdateVarianceFromAccums(hmm_set, mixture->variance);
    }

    if(!mixture->mean->macro) {
      UpdateMeanFromAccums(hmm_set, mixture->mean);
    }
  }
  ComputeGConst(mixture);
}

void UpdateStateFromAccums(HMMSet *hmm_set, State *state, HMM *hmm) {
  int i;

  if(state->outPDF_kind == KID_DiagC) {
    FLOAT accum_sum = 0;

//    if(hmm_set->updateMask & UM_WEIGHT) {
      for(i = 0; i < state->num_mixtures; i++) {
        accum_sum += state->mixture[i].weight_accum;
      }

      if(accum_sum <= 0.0) {
        if(state->macro) {
          Warning("No occupation of '%s', state is not updated",
                  state->macro->name);
        } else {
          int j; // find the state number
          for(j=0; j < hmm->nstates && hmm->state[j] != state; j++);
          Warning("No occupation of '%s[%d]', state is not updated",
                  hmm->macro->name, j + 1);
        }
        return;
      }

      // Remove mixtures with low weight
    if(hmm_set->updateMask & UM_WEIGHT) {
      for(i = 0; i < state->num_mixtures; i++) {
        if(state->mixture[i].weight_accum / accum_sum < hmm_set->minMixWeight) {
          if(hmm) {
            int j; for(j=0; j < hmm->nstates && hmm->state[j] != state; j++); // find the state number
            Warning("Discarding mixture %d of HMM %s state %d because of too low mixture weight",
                    i, hmm->macro->name, j + 1);
          } else {
            assert(state->macro);
            Warning("Discarding mixture %d of state %s because of too low mixture weight",
                    i, state->macro->name);
          }
          accum_sum -= state->mixture[i].weight_accum;

          if(!state->mixture[i].estimates->macro) {
            ScanMixture(state->mixture[i].estimates,mtm_all,NULL,ReleaseItem,NULL);
          }

          state->mixture[i--] = state->mixture[--state->num_mixtures];
          continue;
        }
      }
    }

    for(i = 0; i < state->num_mixtures; i++) {
//      printf("Weight Acc: %f\n", (float) state->mixture[i].weight_accum);
      if(hmm_set->updateMask & UM_WEIGHT) {
        state->mixture[i].weight = log(state->mixture[i].weight_accum / accum_sum);
      }
      if(!state->mixture[i].estimates->macro) {
//!!! This is just too ugly hack
weight_accum_den = state->mixture[i].weight_accum_den;
        UpdateMixtureFromAccums(hmm_set, state->mixture[i].estimates);
      }
    }
  }
}

void UpdateHMMFromAccums(HMMSet *hmm_set, HMM *hmm) {
  int i;

  for(i = 0; i < hmm->nstates - 2; i++) {
    if(!hmm->state[i]->macro) {
      UpdateStateFromAccums(hmm_set, hmm->state[i], hmm);
    }
  }

  if(!hmm->transition->macro) {
    UpdateTransitionFromAccums(hmm_set, hmm->transition);
  }
}

#define WARN_FEW_EXAMPLES(type, name, exs) \
  Warning(type" %s is not updated (%s%ld example%s)", \
          name, exs == 0 ? "" : "only ", exs, exs == 1 ? "" : "s")

void UpdateHMMSetFromAccums(const char *out_dir, HMMSet *hmm_set) {
  Macro *macro;
  int i;
/*char fileName[1024];
  FILE *fp;


  if(hmm_set->updateMask & UM_XFORM) {

    for(j = 0; j < hmm_set->nxformsToUpdate; j++) {
      char *ext;
      LinearXform *xform = (LinearXform *) hmm_set->xformToUpdate[j].xform;
      assert(xform->xform_type == XT_LINEAR);
      MakeFileName(fileName, xform->macro->name, out_dir, NULL);
      ext = fileName + strlen(fileName);
      strcpy(ext, ".xfm");

      if((fp = fopen(fileName, "rt")) == NULL) {
        Error("Cannot open input xforn file %s", fileName);
      }

      for(i=0; i < xform->out_size * xform->in_size; i++) {
        if(fscanf(fp, FLOAT_FMT, &xform->matrix[i]) != 1) {
          Error("Cannot read xform file %s", fileName);
        }
      }

      fclose(fp);
    }

  } */

  for(i = 0; i < hmm_set->hmm_hash.nentries; i++) {
    macro = (Macro *) hmm_set->hmm_hash.entry[i]->data;
//  for(macro = hmm_set->hmm_list; macro != NULL; macro = macro->next) {
    if(*(Macro **)macro->data != macro) continue;
    if(macro->occurances < hmm_set->minOccurances) {
      WARN_FEW_EXAMPLES("Model", macro->name, macro->occurances);
    } else {
      UpdateHMMFromAccums(hmm_set, (HMM *) macro->data);
    }
  }

  for(i = 0; i < hmm_set->state_hash.nentries; i++) {
    macro = (Macro *) hmm_set->state_hash.entry[i]->data;
//  for(macro = hmm_set->state_list; macro != NULL; macro = macro->next) {
    if(*(Macro **)macro->data != macro) continue;
    if(macro->occurances < hmm_set->minOccurances) {
      WARN_FEW_EXAMPLES("State", macro->name, macro->occurances);
    } else {
      UpdateStateFromAccums(hmm_set, (State *) macro->data, NULL);
    }
  }

  for(i = 0; i < hmm_set->mixture_hash.nentries; i++) {
    macro = (Macro *) hmm_set->mixture_hash.entry[i]->data;
//  for(macro = hmm_set->mixture_list; macro != NULL; macro = macro->next) {
    if(*(Macro **)macro->data != macro) continue;
    if(macro->occurances < hmm_set->minOccurances) {
      WARN_FEW_EXAMPLES("Mixture", macro->name, macro->occurances);
    } else {
      UpdateMixtureFromAccums(hmm_set, (Mixture *) macro->data);
    }
  }

  for(i = 0; i < hmm_set->mean_hash.nentries; i++) {
    macro = (Macro *) hmm_set->mean_hash.entry[i]->data;
//  for(macro = hmm_set->mean_list; macro != NULL; macro = macro->next) {
    if(*(Macro **)macro->data != macro) continue;
    if(macro->occurances < hmm_set->minOccurances) {
      WARN_FEW_EXAMPLES("Mean vector", macro->name, macro->occurances);
    } else {
      UpdateMeanFromAccums(hmm_set, (Mean *) macro->data);
    }
  }

  for(i = 0; i < hmm_set->variance_hash.nentries; i++) {
    macro = (Macro *) hmm_set->variance_hash.entry[i]->data;
//  for(macro = hmm_set->variance_list; macro != NULL; macro = macro->next) {
    if(*(Macro **)macro->data != macro) continue;
    if(strcmp(macro->name, "varFloor1")) {
      if(macro->occurances < hmm_set->minOccurances) {
        WARN_FEW_EXAMPLES("Variance vector", macro->name, macro->occurances);
      } else {
        UpdateVarianceFromAccums(hmm_set, (Variance *) macro->data);
      }
    }
  }

  for(i = 0; i < hmm_set->transition_hash.nentries; i++) {
    macro = (Macro *) hmm_set->transition_hash.entry[i]->data;
//  for(macro = hmm_set->transition_list; macro != NULL; macro = macro->next) {
    if(*(Macro **)macro->data != macro) continue;
    if(macro->occurances < hmm_set->minOccurances) {
      WARN_FEW_EXAMPLES("Transition matrix ", macro->name, macro->occurances);
    } else {
      UpdateTransitionFromAccums(hmm_set, (Transition *) macro->data);
    }
  }
}

void DistributeMacroOccurances(HMMSet *hmm_set) {
  Macro *macro;
  int i, j, k;

  for(k = 0; k < hmm_set->hmm_hash.nentries; k++) {
    macro = (Macro *) hmm_set->hmm_hash.entry[k]->data;

    if(*(Macro **)macro->data != macro) { // Not a physical model
      (*(Macro **)macro->data)->occurances += macro->occurances;
    }

    HMM *hmm = (HMM *) macro->data;

    for(i = 0; i < hmm->nstates - 2; i++) {
      State *state = hmm->state[i];

      if(state->macro) {
        state->macro->occurances += macro->occurances;
      }

      if(state->outPDF_kind == KID_DiagC) {
        for(j = 0; j < state->num_mixtures; j++) {
          Mixture *mixture = state->mixture[j].estimates;

          if(mixture->macro) {
            mixture->macro->occurances += macro->occurances;
          }
          if(mixture->mean->macro) {
            mixture->mean->macro->occurances += macro->occurances;
          }
          if(mixture->variance->macro) {
            mixture->variance->macro->occurances += macro->occurances;
          }
        }
      }
    }

    if(hmm->transition->macro) {
      hmm->transition->macro->occurances += macro->occurances;
    }
  }
}

typedef struct {
  FLOAT *observation;
  int time;
} GlobalStats_UserData;

void GlobalStats(int macro_type, HMMSetNodeName nn, void *data, void *userData)
{
  int i;
  GlobalStats_UserData *ud = (GlobalStats_UserData *) userData;

  if(macro_type == mt_state) {
    State * state = (State *) data;
    if(state->outPDF_kind == KID_DiagC) {
      for(i = 0; i < state->num_mixtures; i++) {
        state->mixture[i].weight_accum += 1;
      }
    }
  } else { // macro_type == mt_mixture
    Mixture *mixture = (Mixture *) data;
    int vec_size = mixture->mean->vec_size;
    FLOAT *obs=XformPass(mixture->inputXform,ud->observation,ud->time,FORWARD);
    for(i = 0; i < vec_size; i++) {
      mixture->mean->vector[vec_size + i] += obs[i];
    }
    mixture->mean->vector[2 * vec_size] += 1;

    for(i = 0; i < vec_size; i++) {
      mixture->variance->vector[vec_size  +i] += SQR(obs[i]);
      mixture->variance->vector[2*vec_size+i] += obs[i];
    }
    mixture->variance->vector[3 * vec_size] += 1;
  }
}

void GlobalStatsForHMMSet(HMMSet *hmm_set, FLOAT *observation, int time) {
  GlobalStats_UserData ud = {observation, time};
  ScanHMMSet(hmm_set, mtm_state | mtm_mixture, NULL, GlobalStats, &ud);
}


FLOAT    *LinearXformEval(   LinearXform *, FLOAT *, FLOAT *, char *, PropagDir);
FLOAT      *CopyXformEval(     CopyXform *, FLOAT *, FLOAT *, char *, PropagDir);
FLOAT      *FuncXformEval(     FuncXform *, FLOAT *, FLOAT *, char *, PropagDir);
FLOAT      *BiasXformEval(     BiasXform *, FLOAT *, FLOAT *, char *, PropagDir);
FLOAT  *StackingXformEval( StackingXform *, FLOAT *, FLOAT *, char *, PropagDir);
FLOAT *CompositeXformEval(CompositeXform *, FLOAT *, FLOAT *, char *, PropagDir);

void ResetXformInstances(HMMSet *hmm_set)
{
  XformInstance *inst;
  for(inst = hmm_set->Xform_instances; inst != NULL; inst = inst->next) {
    inst->statCacheTime = UNDEF_TIME;
    inst->time          = UNDEF_TIME;
  }
}


FLOAT *XformEval(Xform *xform, FLOAT *in_vec, FLOAT *out_vec,
                 char *memory, PropagDir dir)
{
  typedef FLOAT * (*pEvalFunc)(Xform *, FLOAT *, FLOAT *, char *, PropagDir);
  XformType type = xform->xform_type;
  pEvalFunc evalFunc = type == XT_LINEAR    ? (pEvalFunc) LinearXformEval    :
                       type == XT_COPY      ? (pEvalFunc) CopyXformEval      :
                       type == XT_FUNC      ? (pEvalFunc) FuncXformEval      :
                       type == XT_BIAS      ? (pEvalFunc) BiasXformEval      :
                       type == XT_STACKING  ? (pEvalFunc) StackingXformEval  :
                       type == XT_COMPOSITE ? (pEvalFunc) CompositeXformEval :
                       NULL;
  assert(evalFunc);
  return evalFunc(xform, in_vec, out_vec, memory, dir);
}

FLOAT *LinearXformEval(LinearXform *xform, FLOAT *in_vec, FLOAT *out_vec,
                       char *memory, PropagDir dir)
{
  int c, r;
  for(r = 0; r < xform->out_size; r++) {
    out_vec[r] = 0.0;
    for(c = 0; c < xform->in_size; c++) {
      out_vec[r] += in_vec[c] * xform->matrix[xform->in_size * r + c];
    }
  }
  return out_vec;
}

FLOAT *CopyXformEval(CopyXform *xform, FLOAT *in_vec, FLOAT *out_vec,
                     char *memory, PropagDir dir)
{
  int i;
  for(i = 0; i < xform->out_size; i++) {
    out_vec[i] = in_vec[xform->indices[i]];
  }
  return out_vec;
}

FLOAT *FuncXformEval(FuncXform *xform, FLOAT *in_vec, FLOAT *out_vec,
                     char *memory, PropagDir dir)
{
  FuncTable[xform->funcId].funcPtr(in_vec, out_vec, xform->out_size);
  return out_vec;
}

FLOAT *BiasXformEval(BiasXform *xform, FLOAT *in_vec, FLOAT *out_vec,
                     char *memory, PropagDir dir)
{
  int i;
  for(i = 0; i < xform->out_size; i++) {
    out_vec[i] = in_vec[i] + xform->vector[i];
  }
  return out_vec;
}



FLOAT *StackingXformEval(StackingXform *xform, FLOAT *in_vec,
                       FLOAT *out_vec, char *memory, PropagDir dir)
{
  FLOAT *stack = (FLOAT *) memory;

  memmove(stack + (dir == BACKWARD ? xform->in_size : 0),
          stack + (dir ==  FORWARD ? xform->in_size : 0),
          (xform->out_size - xform->in_size) * sizeof(FLOAT));
//        out_vec += (stack_size - 1) * out_size; FORWARD
  memmove(stack + (dir ==  FORWARD ? xform->out_size - xform->in_size : 0),
          in_vec, xform->in_size * sizeof(FLOAT));

  if(!xform->horiz_stack) {
    memmove(out_vec, stack, xform->out_size * sizeof(FLOAT));
  } else { // stacking == HORZ_STACK
    int t, c, stack_size = xform->out_size / xform->in_size;
    for(t = 0; t < stack_size; t++) {
      for(c = 0; c < xform->in_size; c++) {
        out_vec[c * stack_size + t] = stack[t * xform->in_size + c];
      }
    }
  }
  return out_vec;
}


FLOAT *CompositeXformEval(CompositeXform *xform, FLOAT *in_vec, FLOAT *out_vec,
                          char *memory, PropagDir dir)
{
  int i, j;

  for(i = 0; i < xform->nlayers; i++) {
    FLOAT *in =  i == 0                ? in_vec  : xform->layer[i-1].out_vec;
    FLOAT *out = i == xform->nlayers-1 ? out_vec : xform->layer[i]  .out_vec;
    
    for(j = 0; j < xform->layer[i].nblocks; j++) {
      XformEval(xform->layer[i].block[j], in, out, memory, dir);
      in  += xform->layer[i].block[j]->in_size;
      out += xform->layer[i].block[j]->out_size;
      memory += xform->layer[i].block[j]->memorySize;
    }
  }
  return out_vec;
}

FLOAT *XformPass(XformInstance *xformInst, FLOAT *in_vec, int time, PropagDir dir)
{
  if(xformInst == NULL) return in_vec;

  if(time != UNDEF_TIME && xformInst->time == time) return xformInst->out_vec;

  xformInst->time = time;

  if(xformInst->input) {
    in_vec = XformPass(xformInst->input, in_vec, time, dir);
  }

  XformEval(xformInst->xform,in_vec,xformInst->out_vec,xformInst->memory,dir);

/*  {int i;
  for(i=0; i<xformInst->out_size; i++)
    printf("%.2f ", xformInst->out_vec[i]);
  printf("%s\n", dir == FORWARD ? "FORWARD" : "BACKWARD");}*/

  return xformInst->out_vec;
}

void UpdateStacks(HMMSet *hmm_set, FLOAT *obs, int time,  PropagDir dir) {
  XformInstance *inst;
  for(inst = hmm_set->Xform_instances; inst != NULL; inst = inst->next) {
    if(inst->xform->delay > 0) {
      XformPass(inst, obs, time, dir);
    }
  }
}

void ScanHMMSet(HMMSet *hmm_set, int mask, HMMSetNodeName nodeName,
                ScanAction action, void *userData)
{
  Macro *macro;
  if(nodeName != NULL) strcpy(nodeName+sizeof(HMMSetNodeName)-4, "...");

  for(macro = mask & mtm_revpass ? hmm_set->last_macro : hmm_set->first_macro;
      macro != NULL;
      macro = mask & mtm_revpass ? macro->prevAll      : macro->nextAll) {
    if(*(Macro **) macro->data != macro) continue;
    if(nodeName != NULL) {
      strncpy(nodeName, macro->name, sizeof(HMMSetNodeName)-4);
    }

    switch(macro->type) {
      case mt_Xform:
        if(!(mask & mtm_Xform)) break;
        ScanXform(static_cast <Xform*> (macro->data), mask, nodeName, action, userData);
        break;
      case mt_XformInstance:
        if(!(mask & (mtm_Xform | mtm_XformInstance))) break;
        ScanXformInstance(static_cast <XformInstance*> (macro->data), mask, nodeName, action, userData);
        break;
      case mt_mean:
        if(!(mask & mtm_mean)) break;
        action(mt_mean, nodeName, macro->data, userData);
        break;
      case mt_variance:
        if(!(mask & mtm_variance)) break;
        action(mt_variance, nodeName, macro->data, userData);
        break;
      case mt_transition:
        if(!(mask & mtm_transition)) break;
        action(mt_transition, nodeName, macro->data, userData);
        break;
      case mt_mixture:
        if(!(mask & (mtm_all & ~(mtm_state | mtm_hmm | mtm_transition)))) break;
        ScanMixture(static_cast <Mixture*> (macro->data), mask, nodeName, action, userData);
        break;
      case mt_state:
        if(!(mask & (mtm_all & ~(mtm_hmm | mtm_transition)))) break;
        ScanState(static_cast <State*> (macro->data), mask, nodeName, action, userData);
        break;
      case mt_hmm:
        if(!(mask & mtm_all)) break;
        ScanHMM(static_cast <HMM*> (macro->data), mask, nodeName, action, userData);
        break;
      default: assert(0);
    }
  }
}

void ScanHMM(HMM *hmm, int mask, HMMSetNodeName nodeName,
             ScanAction action, void *userData)
{
  int i, n = 0;
  char *chptr = NULL;

  if(nodeName != NULL) {
   n = strlen(nodeName);
   chptr = nodeName + n;
   n = sizeof(HMMSetNodeName) - 4 - n;
  }

  if(mask & mtm_hmm && mask & mtm_prescan) {
    action(mt_hmm, nodeName, hmm, userData);
  }

  if(mask & (mtm_all & ~(mtm_hmm | mtm_transition))) {
    for(i=0; i < hmm->nstates-2; i++) {
      if(!hmm->state[i]->macro) {
        if(n > 0 ) snprintf(chptr, n, ".state[%d]", i+2);
        ScanState(hmm->state[i], mask, nodeName, action, userData);
      }
    }
  }

  if(mask & mtm_transition && !hmm->transition->macro) {
    if(n > 0) strncpy(chptr, ".transP", n);
    action(mt_transition, nodeName, hmm->transition, userData);
  }

  if(mask & mtm_hmm && !(mask & mtm_prescan)) {
    if(n > 0) chptr = '\0';
    action(mt_hmm, nodeName, hmm, userData);
  }
}

void ScanState(State *state, int mask, HMMSetNodeName nodeName,
               ScanAction action, void *userData)
{
  int i, n = 0;
  char *chptr = NULL;

  if(nodeName != NULL) {
   n = strlen(nodeName);
   chptr = nodeName + n;
   n = sizeof(HMMSetNodeName) - 4 - n;
  }

  if(mask & mtm_state && mask & mtm_prescan) {
    action(mt_state, nodeName, state, userData);
  }

  if(state->outPDF_kind != KID_PDFObsVec &&
     mask & (mtm_all & ~(mtm_state | mtm_hmm | mtm_transition))) {
    for(i=0; i < state->num_mixtures; i++) {
      if(!state->mixture[i].estimates->macro) {
        if(n > 0 ) snprintf(chptr, n, ".mix[%d]", i+1);
        ScanMixture(state->mixture[i].estimates, mask, nodeName,
                    action, userData);
      }
    }
  }
  if(mask & mtm_state && !(mask & mtm_prescan)) {
    if(n > 0) chptr = '\0';
    action(mt_state, nodeName, state, userData);
  }
}

void ScanMixture(Mixture *mixture, int mask,
                 HMMSetNodeName nodeName, ScanAction action, void *userData)
{
  int n = 0;
  char *chptr = NULL;
  
  if(nodeName != NULL) {
   n = strlen(nodeName);
   chptr = nodeName + n;
   n = sizeof(HMMSetNodeName) - 4 - n;
  }

  if(mask & mtm_mixture && mask & mtm_prescan) {
    action(mt_mixture, nodeName, mixture, userData);
  }

  if(mask & mtm_mean && !mixture->mean->macro) {
    if(n > 0) strncpy(chptr, ".mean", n);
    action(mt_mean, nodeName, mixture->mean, userData);
  }

  if(mask & mtm_variance && !mixture->variance->macro) {
    if(n > 0) strncpy(chptr, ".cov", n);
    action(mt_variance, nodeName, mixture->variance, userData);
  }

  if(mask & mtm_XformInstance && mixture->inputXform &&
     !mixture->inputXform->macro) {
    if(n > 0) strncpy(chptr, ".input", n);
    ScanXformInstance(mixture->inputXform, mask, nodeName, action, userData);
  }

  if(mask & mtm_mixture && !(mask & mtm_prescan)) {
    if(n > 0) chptr = '\0';
    action(mt_mixture, nodeName, mixture, userData);
  }
}

void ScanXformInstance(XformInstance *xformInstance, int mask,
                       HMMSetNodeName nodeName, ScanAction action,
                       void *userData)
{
  int n = 0;
  char *chptr = NULL;
  
  if(nodeName != NULL) {
   n = strlen(nodeName);
   chptr = nodeName + n;
   n = sizeof(HMMSetNodeName) - 4 - n;
  }

  if(mask & mtm_XformInstance && mask & mtm_prescan) {
    action(mt_XformInstance, nodeName, xformInstance, userData);
  }
  
  if(xformInstance->input != NULL) {
    if(!xformInstance->input->macro) {
    if(n > 0) strncpy(chptr, ".input", n);
      ScanXformInstance(xformInstance->input, mask, nodeName, action, userData);
    }
  }

  if(mask & mtm_Xform && xformInstance->xform != NULL && 
     !xformInstance->xform->macro) {
    if(n > 0) strncpy(chptr, ".xform", n);
    ScanXform(xformInstance->xform, mask, nodeName, action, userData);
  }
  
  if(mask & mtm_XformInstance && !(mask & mtm_prescan)) {
    if(n > 0) chptr = '\0';
    action(mt_XformInstance, nodeName, xformInstance, userData);
  }
}

void ScanXform(Xform *xform, int mask, HMMSetNodeName nodeName,
               ScanAction action, void *userData)
{
  int n = 0;
  char *chptr = NULL;

  if(nodeName != NULL) {
   n = strlen(nodeName);
   chptr = nodeName + n;
   n = sizeof(HMMSetNodeName) - 4 - n;
  }

  if(mask & mtm_prescan) {
    action(mt_Xform, nodeName, xform, userData);
  }
  
  if(xform->xform_type == XT_COMPOSITE) {
    CompositeXform *cxf = (CompositeXform *) xform;
    int i, j;

    for(i=0; i < cxf->nlayers; i++) {
      for(j = 0; j < cxf->layer[i].nblocks; j++) {
        if(!cxf->layer[i].block[j]->macro) {
          if(n > 0) snprintf(chptr, n, ".part[%d,%d]", i+1, j+1);
          ScanXform(cxf->layer[i].block[j], mask, nodeName, action, userData);
        }
      }
    }
  }

  if(!(mask & mtm_prescan)) {
    if(n > 0) chptr = '\0';
    action(mt_Xform, nodeName, xform, userData);
  }
}


enum StatType {MEAN_STATS, COV_STATS};
void AllocXformStatAccums(XformStatAccum **xformStatAccum,
                          int *nxformStatAccums,
                          XformInstance *xformInstance,
                          enum StatType stat_type) {
  int i, j;
  if(xformInstance == NULL) return;

  for(i = 0; i < xformInstance->nxformStatCaches; i++) {
    XformStatCache *xfsc = &xformInstance->xformStatCache[i];
    XformStatAccum *xfsa = *xformStatAccum;

    for(j = 0; j < *nxformStatAccums; j++, xfsa++) {
      if(xfsa->xform == xfsc->xform) break;
    }

    if(j == *nxformStatAccums) {
      int size = xfsc->xform->in_size; //mean : mean+covariance
      size = (stat_type == MEAN_STATS) ? size : size+size*(size+1)/2;

      *xformStatAccum =
        (XformStatAccum *) realloc(*xformStatAccum,
                                   sizeof(XformStatAccum) * ++*nxformStatAccums);

      if(*xformStatAccum == NULL) {
        Error("Insufficient memory");
      }

      xfsa = *xformStatAccum + *nxformStatAccums - 1;

      if((xfsa->stats = (FLOAT *) malloc(sizeof(FLOAT) * size)) == NULL) {
        Error("Insufficient memory");
      }

      xfsa->xform    = xfsc->xform;
      xfsa->norm     = 0.0;
      for(j = 0; j < size; j++) xfsa->stats[j] = 0.0;
    }
  }
}

void AllocateXformStatCachesAndAccums(int macro_type, HMMSetNodeName nodeName,
                                      void *data, void *hmm_set)
{
  if(macro_type == mt_XformInstance) {
    //Allocate Xform stat caches for XformInstance

    XformInstance *xfi =(XformInstance *) data;
    int i, j;

    for(i=0; i < ((HMMSet *) hmm_set)->nxformsToUpdate; i++) {
      Xform *xform = ((HMMSet *) hmm_set)->xformToUpdate[i].xform;
      int instanceContainXfrom = IsXformIn1stLayer(xform, xfi->xform);

      //Does instance one level up contain cache for this xform
      XformStatCache *upperLevelStats = NULL;
      if(xfi->input != NULL) {
        for(j=0; j < xfi->input->nxformStatCaches; j++) {
          if(xfi->input->xformStatCache[j].xform == xform) {
            upperLevelStats = &xfi->input->xformStatCache[j];
            break;
          }
        }
      }

      if(instanceContainXfrom || upperLevelStats != NULL) {
        XformStatCache *xfsc;

        xfi->xformStatCache = (XformStatCache *)
          realloc(xfi->xformStatCache,
                  sizeof(XformStatCache) * ++xfi->nxformStatCaches);

        if(xfi->xformStatCache == NULL) {
          Error("Insufficient memory");
        }

        xfsc = &xfi->xformStatCache[xfi->nxformStatCaches-1];

        if(instanceContainXfrom) {
          int size = xform->in_size;
          size = size+size*(size+1)/2;

          if((xfsc->stats = (FLOAT *) malloc(sizeof(FLOAT) * size))==NULL) {
            Error("Insufficient memory");
          }
        } else {
          xfsc->stats = upperLevelStats->stats;
        }

        xfsc->norm = 0;
        xfsc->xform = xform;
        xfsc->upperLevelStats = upperLevelStats;
      }
    }
  } else if(macro_type == mt_mixture) {
    //Allocate Xform stat accumulators for mean and covariance

    Mixture *mix = (Mixture *) data;
    AllocXformStatAccums(&mix->mean->xformStatAccum,
                         &mix->mean->nxformStatAccums,
                         mix->inputXform, MEAN_STATS);

    AllocXformStatAccums(&mix->variance->xformStatAccum,
                         &mix->variance->nxformStatAccums,
                         mix->inputXform, COV_STATS);

    if(mix->inputXform == NULL || mix->inputXform->nxformStatCaches == 0)
      return;
      
    if(mix->inputXform->nxformStatCaches != 1 ||
       !Is1Layer1BlockLinearXform(mix->inputXform->xform) ||
       mix->inputXform->xformStatCache[0].upperLevelStats != NULL) {
      mix->variance->updatableFromStatAccums = 0;
      mix->mean    ->updatableFromStatAccums = 0;
      ((HMMSet *) hmm_set)->allMixuresUpdatableFromStatAccums = 0;  
    } else if(mix->mean->nxformStatAccums != 1) {
      assert(mix->mean->nxformStatAccums > 1);
      mix->mean->updatableFromStatAccums = 0;
      ((HMMSet *) hmm_set)->allMixuresUpdatableFromStatAccums = 0;
    } else if(mix->variance->nxformStatAccums != 1) {
      assert(mix->variance->nxformStatAccums > 1);
      mix->variance->updatableFromStatAccums = 0;
      ((HMMSet *) hmm_set)->allMixuresUpdatableFromStatAccums = 0;
    }
  }
}

void AllocateAccumulatorsForXformStats(HMMSet *hmm_set) {
  hmm_set->allMixuresUpdatableFromStatAccums = 1;
  ScanHMMSet(hmm_set, mtm_XformInstance | mtm_mixture, NULL,
             AllocateXformStatCachesAndAccums, hmm_set);
}

void NormalizeStatsForXform(int macro_type, HMMSetNodeName nodeName,
                            void *data, void *userData) {
  XformStatAccum *xfsa = NULL;
  int i, j, k, nxfsa = 0, size;
  FLOAT *mean, *cov, inorm;

  if(macro_type == mt_mean) {
    xfsa  = ((Mean *)data)->xformStatAccum;
    nxfsa = ((Mean *)data)->nxformStatAccums;
  } else if(macro_type == mt_variance) {
    xfsa  = ((Variance *)data)->xformStatAccum;
    nxfsa = ((Variance *)data)->nxformStatAccums;
  }

  for(i = 0; i < nxfsa; i++) {
    size = xfsa[i].xform->in_size;
    mean = xfsa[i].stats;
    cov  = xfsa[i].stats + size;
    inorm = 1.0 / xfsa[i].norm;

    for(j = 0; j < size; j++) mean[j] *= inorm; //normalize means

    if(macro_type == mt_variance) {
      for(k=0; k < size; k++) {
        for(j=0; j <= k; j++) {                 //normalize covariances
          cov[k*(k+1)/2+j] = cov[k*(k+1)/2+j] * inorm - mean[k] * mean[j];
        }
      }
    }
  }
}

typedef struct {
  FILE *statsP;
  char *statsN;
  FILE *occupP;
  char *occupN;
} WriteStatsForXform_Files;

typedef struct {
  LinearXform *xform;
  WriteStatsForXform_Files meanFile;
  WriteStatsForXform_Files covFile;
  int binary;
} WriteStatsForXform_UserData;

void WriteStatsForXform(int macro_type, HMMSetNodeName nodeName,
                        void *data, void *userData) {
  WriteStatsForXform_Files *file = NULL;
  XformStatAccum *xfsa = NULL;
  int i, j, k, nxfsa = 0, cc = 0, size;
  FLOAT *mean, *cov;
  WriteStatsForXform_UserData *ud = (WriteStatsForXform_UserData *) userData;

  if(macro_type == mt_mean) {
    file  = &ud->meanFile;
    xfsa  = ((Mean *)data)->xformStatAccum;
    nxfsa = ((Mean *)data)->nxformStatAccums;
  } else if(macro_type == mt_variance) {
    file  = &ud->covFile;
    xfsa  = ((Variance *)data)->xformStatAccum;
    nxfsa = ((Variance *)data)->nxformStatAccums;
  }

  for(i = 0; i < nxfsa && xfsa[i].xform != (Xform *) ud->xform; i++);
  if(i == nxfsa) return;

  if(fprintf(file->occupP, "%s "FLOAT_FMT"\n", nodeName, xfsa[i].norm) < 0) {
    Error("Cannot write to file: %s", file->occupN);
  }

  size = xfsa[i].xform->in_size;
  mean = xfsa[i].stats;
  cov  = xfsa[i].stats + size;

  if(macro_type == mt_mean) {
    if(ud->binary) {
      if(!isBigEndian()) for(i = 0; i < size; i++) swapFLOAT(mean[i]);
      cc |= fwrite(mean, sizeof(FLOAT), size, file->statsP) != size;
      if(!isBigEndian()) for(i = 0; i < size; i++) swapFLOAT(mean[i]);
    } else {
      for(j=0;j<size;j++) {
        cc |= fprintf(file->statsP, FLOAT_FMT" ", mean[j]) < 0;
      }
      cc |= fputs("\n", file->statsP) == EOF;
    }
  } else {
    if(ud->binary) {
      size = size*(size+1)/2;
      if(!isBigEndian()) for(i = 0; i < size; i++) swapFLOAT(cov[i]);
      cc |= fwrite(cov, sizeof(FLOAT), size, file->statsP) != size;
      if(!isBigEndian()) for(i = 0; i < size; i++) swapFLOAT(cov[i]);
    } else{
      for(k=0; k < size; k++) {
        for(j=0;j<=k;j++) {
          cc |= fprintf(file->statsP, FLOAT_FMT" ", cov[k*(k+1)/2+j]) < 0;
        }

        for(;j<size; j++) {
          cc |= fprintf(file->statsP, FLOAT_FMT" ", cov[j*(j+1)/2+k]) < 0;
        }

        cc |= fputs("\n", file->statsP) == EOF;
      }
      cc |= fputs("\n", file->statsP) == EOF;
    }
  }

  if(cc) {
    Error("Cannot write to file %s", file->statsN);
  }
}

void WriteXformStatsAndRunCommands(const char *out_dir, int binary,
                                   HMMSet *hmm_set) {
  HMMSetNodeName nodeNameBuffer;
  WriteStatsForXform_UserData userData;
  char fileName[1024];
  int i, j, k;
  FILE *fp;

  struct XfStatsHeader {
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
  userData.binary = binary;
  for(k = 0; k < hmm_set->nxformsToUpdate; k++) {
    char *ext, *shellCommand = hmm_set->xformToUpdate[k].shellCommand;
    userData.xform = (LinearXform *) hmm_set->xformToUpdate[k].xform;
    assert(userData.xform->xform_type == XT_LINEAR);

    MakeFileName(fileName, userData.xform->macro->name, out_dir, NULL);
    ext = fileName + strlen(fileName);
    strcpy(ext, ".xms"); userData.meanFile.statsN = strdup(fileName);
    strcpy(ext, ".xmo"); userData.meanFile.occupN = strdup(fileName);
    strcpy(ext, ".xcs"); userData.covFile.statsN = strdup(fileName);
    strcpy(ext, ".xco"); userData.covFile.occupN = strdup(fileName);

    if(userData.meanFile.statsN == NULL || userData.meanFile.occupN == NULL||
       userData.covFile.statsN  == NULL || userData.covFile.occupN  == NULL) {
      Error("Insufficient memory");
    }

    userData.meanFile.statsP = fopen(userData.meanFile.statsN, binary?"w":"wt");
    if(userData.meanFile.statsP == NULL) {
      Error("Cannot open output file %s",
            userData.meanFile.statsN);
    }

    if(binary) {
      header.stats_type = STATS_MEAN;
      header.size = userData.xform->in_size;
      if(!isBigEndian()) swap2(header.size);
      if(fwrite(&header, sizeof(header), 1, userData.meanFile.statsP) != 1) {
        Error("Cannot write to file: %s", userData.meanFile.statsN);
      }
    }

    userData.meanFile.occupP = fopen(userData.meanFile.occupN, "wt");
    if(userData.meanFile.occupP == NULL) {
      Error("Cannot open output file %s", userData.meanFile.occupN);
    }

    userData.covFile.statsP = fopen(userData.covFile.statsN, binary?"w":"wt");
    if(userData.covFile.statsP == NULL) {
      Error("Cannot open output file %s", userData.covFile.statsN);
    }

    if(binary) {
      header.stats_type = STATS_COV_LOW_TRI;
      header.size = userData.xform->in_size;
      if(!isBigEndian()) swap2(header.size);
      if(fwrite(&header, sizeof(header), 1, userData.covFile.statsP) != 1) {
        Error("Cannot write to file: %s",
              userData.covFile.statsN);
      }
    }

    userData.covFile.occupP = fopen(userData.covFile.occupN, "wt");
    if(userData.covFile.occupP == NULL) {
      Error("Cannot open output file %s", userData.covFile.occupN);
    }

    ScanHMMSet(hmm_set, mtm_mean | mtm_variance, nodeNameBuffer,
               WriteStatsForXform, &userData);


    fclose(userData.meanFile.statsP); fclose(userData.meanFile.occupP);
    fclose(userData.covFile.statsP);  fclose(userData.covFile.occupP);
    free(userData.meanFile.statsN);   free(userData.meanFile.occupN);
    free(userData.covFile.statsN);    free(userData.covFile.occupN);

    strcpy(ext, ".xfm");
    if((fp = fopen(fileName, "wt")) == NULL) {
      Error("Cannot open output file %s", fileName);
    }

    for(i=0; i < userData.xform->out_size; i++) {
      for(j=0; j < userData.xform->in_size; j++) {
        fprintf(fp, " "FLOAT_FMT,
                userData.xform->matrix[i*userData.xform->in_size+j]);
      }
      fputs("\n", fp);
    }

    fclose(fp);

    if(shellCommand != NULL && *shellCommand) {
      TraceLog("Executing command: %s", shellCommand);
      system(shellCommand);
    }
  }
}

void ReadStatsForXform(int macro_type, HMMSetNodeName nodeName,
                        void *data, void *userData) {
  char buff[128];
  WriteStatsForXform_Files *file = NULL;
  XformStatAccum *xfsa = NULL;
  int i, j, k, nxfsa = 0, cc = 0, size;
  FLOAT *mean, *cov, f;
  WriteStatsForXform_UserData *ud = (WriteStatsForXform_UserData *) userData;

  if(macro_type == mt_mean) {
    file  = &ud->meanFile;
    xfsa  = ((Mean *)data)->xformStatAccum;
    nxfsa = ((Mean *)data)->nxformStatAccums;
  } else if(macro_type == mt_variance) {
    file  = &ud->covFile;
    xfsa  = ((Variance *)data)->xformStatAccum;
    nxfsa = ((Variance *)data)->nxformStatAccums;
  }

  for(i = 0; i < nxfsa && xfsa[i].xform != (Xform *) ud->xform; i++);
  if(i == nxfsa) return;

  j = fscanf(file->occupP, "%128s "FLOAT_FMT"\n", buff, &xfsa[i].norm);
  if(j < 1) {
    Error("Unexpected end of file: %s", file->occupN);
  } else if(strcmp(buff, nodeName)) {
    Error("'%s' expected but '%s' found in file: %s",nodeName,buff,file->occupN);
  } else if(j < 2) {
    Error("Decimal number expected after '%s'in file: %s", buff, file->occupN);
  }
  size = xfsa[i].xform->in_size;
  mean = xfsa[i].stats;
  cov  = xfsa[i].stats + size;

  if(macro_type == mt_mean) {
    if(ud->binary) {
      j = fread(mean, sizeof(FLOAT), size, file->statsP);
      cc |= j != size;
      if(!isBigEndian()) for(i = 0; i < size; i++) swapFLOAT(mean[i]);
    } else {
      for(j=0;j<size;j++) {
        cc |= fscanf(file->statsP, FLOAT_FMT" ", &mean[j]) != 1;
      }
    }
  } else {
    if(ud->binary) {
      size = size*(size+1)/2;
      cc |= fread(cov, sizeof(FLOAT), size, file->statsP) != size;
      if(!isBigEndian()) for(i = 0; i < size; i++) swapFLOAT(cov[i]);
    } else{
      for(k=0; k < size; k++) {
        for(j=0;j<k;j++) {
          cc |= fscanf(file->statsP, FLOAT_FMT" ", &f) != 1;
          if(f != cov[k*(k+1)/2+j]) {
            Error("Covariance matrix '%s' in file '%s' must be symetric",
                  nodeName, file->statsP);
          }
        }

        for(;j<size; j++) {
          cc |= fscanf(file->statsP, FLOAT_FMT" ", &cov[j*(j+1)/2+k]) != 1;
        }
      }
    }
  }

  if(ferror(file->statsP)) {
    Error("Cannot read file '%s'", file->statsN);
  } else if(cc) {
    Error("Invalid file with Xform statistics '%s'", file->statsN);
  }
}

void ReadXformStats(const char *out_dir, int binary, HMMSet *hmm_set) {
  HMMSetNodeName nodeNameBuffer;
  WriteStatsForXform_UserData userData;
  char fileName[1024];
  int i, k;
  FILE *fp;

  struct XfStatsHeader {
    char precision;
    char stats_type;
    short size;
  } header;

  userData.binary = binary;
  for(k = 0; k < hmm_set->nxformsToUpdate; k++) {
    char *ext;
    userData.xform = (LinearXform *) hmm_set->xformToUpdate[k].xform;
    assert(userData.xform->xform_type == XT_LINEAR);

    MakeFileName(fileName, userData.xform->macro->name, out_dir, NULL);
    ext = fileName + strlen(fileName);
    strcpy(ext, ".xms"); userData.meanFile.statsN = strdup(fileName);
    strcpy(ext, ".xmo"); userData.meanFile.occupN = strdup(fileName);
    strcpy(ext, ".xcs"); userData.covFile.statsN = strdup(fileName);
    strcpy(ext, ".xco"); userData.covFile.occupN = strdup(fileName);

    if(userData.meanFile.statsN == NULL || userData.meanFile.occupN == NULL||
       userData.covFile.statsN  == NULL || userData.covFile.occupN  == NULL) {
      Error("Insufficient memory");
    }

    userData.meanFile.statsP = fopen(userData.meanFile.statsN, "r");
    if(userData.meanFile.statsP == NULL) {
      Error("Cannot open input file '%s'", userData.meanFile.statsN);
    }

    if(binary) {
      header.precision = -1;
      fread(&header, sizeof(header), 1, userData.meanFile.statsP);
      if(!isBigEndian()) swap2(header.size);
      if(ferror(userData.meanFile.statsP)) {
        Error("Cannot read input file '%s'", userData.meanFile.statsN);
      } else if(header.stats_type != STATS_MEAN
             || header.size       != userData.xform->in_size
             || header.precision  != (sizeof(FLOAT) == sizeof(float)
                                     ? PRECISION_FLOAT : PRECISION_DOUBLE)) {
        Error("Invalid header in file '%s'", userData.meanFile.statsN);
      }
    }

    userData.meanFile.occupP = fopen(userData.meanFile.occupN, "r");
    if(userData.meanFile.occupP == NULL) {
      Error("Cannot open input file '%s'", userData.meanFile.occupN);
    }

    userData.covFile.statsP = fopen(userData.covFile.statsN, "r");
    if(userData.covFile.statsP == NULL) {
      Error("Cannot open input file '%s'", userData.covFile.statsN);
    }

    if(binary) {
      header.precision = -1;
      fread(&header, sizeof(header), 1, userData.covFile.statsP);
      if(!isBigEndian()) swap2(header.size);
      if(ferror(userData.covFile.statsP)) {
        Error("Cannot read input file '%s'", userData.covFile.statsN);
      } else if(header.stats_type != STATS_COV_LOW_TRI
             || header.size       != userData.xform->in_size
             || header.precision  != (sizeof(FLOAT) == sizeof(float)
                                     ? PRECISION_FLOAT : PRECISION_DOUBLE)) {
        Error("Invalid header in file '%s'", userData.covFile.statsN);
      }
    }

    userData.covFile.occupP = fopen(userData.covFile.occupN, "r");
    if(userData.covFile.occupP == NULL) {
      Error("Cannot open output file '%s'", userData.covFile.occupN);
    }

    ScanHMMSet(hmm_set, mtm_mean | mtm_variance, nodeNameBuffer,
               ReadStatsForXform, &userData);


    fclose(userData.meanFile.statsP); fclose(userData.meanFile.occupP);
    fclose(userData.covFile.statsP);  fclose(userData.covFile.occupP);
    free(userData.meanFile.statsN);   free(userData.meanFile.occupN);
    free(userData.covFile.statsN);    free(userData.covFile.occupN);


    strcpy(ext, ".xfm");
    if((fp = fopen(fileName, "r")) == NULL) {
      Error("Cannot open input xforn file '%s'", fileName);
    }

    int c =0;
    for(i=0; i < userData.xform->out_size * userData.xform->in_size; i++) {
      c |= fscanf(fp, FLOAT_FMT, &userData.xform->matrix[i]) != 1;
    }
    if(ferror(fp)) {
      Error("Cannot read xform file '%s'", fileName);
    } else if(c) {
      Error("Invalid xform file '%s'", fileName);
    }
    fclose(fp);
  }
}


void WriteHMMStats(const char *stat_file, HMMSet *hmm_set)
{
  Macro *macro;
  int i = 0, j, k;
  FILE *fp;

  if((fp = fopen(stat_file, "wt")) == NULL) {
    Error("Cannot open output file: '%s'", stat_file);
  }

  for(macro = hmm_set->first_macro; macro != NULL; macro = macro->nextAll) {
    if(*(Macro **)macro->data != macro) continue;
    if(macro->type != mt_hmm) continue;
    HMM *hmm = (HMM *) macro->data;

    fprintf(fp, "%4d%*c\"%s\" %4ld ", ++i,
                HIGHER_OF(0,13-strlen(macro->name)), ' ',
                macro->name, macro->occurances);

    for(j = 0; j < hmm->nstates-2; j++) {
      State *state = hmm->state[j];
      FLOAT stOccP = 0;
      for(k = 0; k < state->num_mixtures; k++) {
         stOccP += state->mixture[k].weight_accum;
      }
      fprintf(fp, " %10.6f", stOccP);
    }
    fputs("\n", fp);
  }

  fclose(fp);
}

typedef struct {
  FILE *fp;
  char *fn;
  int  mmi;
} WriteAccumUserData;

void WriteAccum(int macro_type, HMMSetNodeName nodeName,
                void *data, void *userData) {
  int i, size;
  FLOAT *vector = NULL;
//  FILE *fp = (FILE *) userData;
  WriteAccumUserData *ud = (WriteAccumUserData * ) userData;
  Macro *macro = *(Macro **) data;

  if(macro &&
    (fprintf(ud->fp, "~%c \"%s\"", macro->type, macro->name) < 0 ||
     fwrite(&macro->occurances, sizeof(macro->occurances), 1, ud->fp) != 1)) {
    Error("Cannot write accumulators to file: '%s'", ud->fn);
  }

  if(macro_type == mt_mean || macro_type == mt_variance) {
    XformStatAccum *xfsa = NULL;
    int nxfsa = 0;

    if(macro_type == mt_mean) {
      xfsa   = ((Mean *)data)->xformStatAccum;
      nxfsa  = ((Mean *)data)->nxformStatAccums;
      size   = ((Mean *)data)->vec_size;
      vector = ((Mean *)data)->vector+size;
      size   = size + 1;
    } else if(macro_type == mt_variance) {
      xfsa   = ((Variance *)data)->xformStatAccum;
      nxfsa  = ((Variance *)data)->nxformStatAccums;
      size   = ((Variance *)data)->vec_size;
      vector = ((Variance *)data)->vector+size;
      size   = size * 2 + 1;
    }

//    if(ud->mmi) vector += size; // Move to MMI accums, which follows ML accums

    if(fwrite(vector, sizeof(FLOAT), size, ud->fp) != size ||
       fwrite(&nxfsa, sizeof(nxfsa),    1, ud->fp) != 1) {
      Error("Cannot write accumulators to file: '%s'", ud->fn);
    }

//    if(!ud->mmi) { // MMI estimation of Xform statistics has not been implemented yet
    for(i = 0; i < nxfsa; i++) {
      size = xfsa[i].xform->in_size;
      size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
      assert(xfsa[i].xform->macro != NULL);
      if(fprintf(ud->fp, "\"%s\"", xfsa[i].xform->macro->name) < 0 ||
        fwrite(&size,         sizeof(int),      1, ud->fp) != 1    ||
        fwrite(xfsa[i].stats, sizeof(FLOAT), size, ud->fp) != size ||
        fwrite(&xfsa[i].norm, sizeof(FLOAT),    1, ud->fp) != 1) {
        Error("Cannot write accumulators to file: '%s'", ud->fn);
      }
    }
//    }
  } else if(macro_type == mt_state) {
    State *state = (State *) data;
    if(state->outPDF_kind == KID_DiagC) {
      for(i = 0; i < state->num_mixtures; i++) {
        if(fwrite(&state->mixture[i].weight_accum,
                  sizeof(FLOAT), 1, ud->fp) != 1 ||
           fwrite(&state->mixture[i].weight_accum_den,
                  sizeof(FLOAT), 1, ud->fp) != 1) {
          Error("Cannot write accumulators to file: '%s'", ud->fn);
        }
      }
    }
  } else if(macro_type == mt_transition) {
    size   = SQR(((Transition *) data)->nstates);
    vector = ((Transition *) data)->matrix + size;

    if(fwrite(vector, sizeof(FLOAT), size, ud->fp) != size) {
      Error("Cannot write accumulators to file: '%s'", ud->fn);
    }
  }
}


void WriteAccums(const char *accfName, const char *out_dir, HMMSet *hmm_set,
                 long totFrames, FLOAT totLogLike) {
  FILE *fp;
  char fileName[1024];
  MakeFileName(fileName, accfName, out_dir, NULL);
  WriteAccumUserData ud;

  if((fp = fopen(fileName, "wb")) == NULL) {
    Error("Cannot open output file: '%s'", fileName);
  }

  if(fwrite(&totFrames,  sizeof(long),  1, fp) != 1 ||
     fwrite(&totLogLike, sizeof(FLOAT), 1, fp) != 1) {
    Error("Cannot write accumulators to file: '%s'", fileName);
  }

  ud.fp  = fp;
  ud.fn  = fileName;
//  ud.mmi = MMI_denominator_accums;

  ScanHMMSet(hmm_set, mtm_prescan | (mtm_all & ~(mtm_XformInstance|mtm_Xform)),
             NULL, WriteAccum, &ud);

  fclose(fp);
}

void NormalizeAccum(int macro_type, HMMSetNodeName nodeName,
                   void *data, void *userData) {
  int i, j, size;
  FLOAT *vector = NULL;

  if(macro_type == mt_mean || macro_type == mt_variance) {
    XformStatAccum *xfsa = NULL;
    int nxfsa = 0;

    if(macro_type == mt_mean) {
      xfsa   = ((Mean *)data)->xformStatAccum;
      nxfsa  = ((Mean *)data)->nxformStatAccums;
      size   = ((Mean *)data)->vec_size;
      vector = ((Mean *)data)->vector+size;
      size   = size + 1;
    } else if(macro_type == mt_variance) {
      xfsa   = ((Variance *)data)->xformStatAccum;
      nxfsa  = ((Variance *)data)->nxformStatAccums;
      size   = ((Variance *)data)->vec_size;
      vector = ((Variance *)data)->vector+size;
      size   = size * 2 + 1;
    }

    for(i=0; i < size; i++) vector[i] /= vector[size-1];

    for(i = 0; i < nxfsa; i++) {
      size = xfsa[i].xform->in_size;
      size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;

      for(j=0; j < size; j++) xfsa[i].stats[j] /= xfsa[i].norm;
      xfsa[i].norm = 1.0;
    }
  } else if(macro_type == mt_state) {
    State *state = (State *) data;
    if(state->outPDF_kind == KID_DiagC) {
      FLOAT accum_sum = 0.0;

      for(i = 0; i < state->num_mixtures; i++)
        accum_sum += state->mixture[i].weight_accum;

      if(accum_sum > 0.0) {
        for(i = 0; i < state->num_mixtures; i++)
          state->mixture[i].weight_accum /= accum_sum;
      }
    }
  } else if(macro_type == mt_transition) {
    int nstates = ((Transition *) data)->nstates;
    vector = ((Transition *) data)->matrix + SQR(nstates);

    for(i=0; i < nstates; i++) {
      FLOAT nrm = LOG_0;
      for(j=0; j < nstates; j++) {
        LOG_INC(nrm, vector[i * nstates + j]);
      }
      if(nrm < LOG_MIN) nrm = 0.0;
      for(j=0; j < nstates; j++) {
        vector[i * nstates + j] -= nrm;
      }
    }
  }
}

void NormalizeAccums(HMMSet *hmm_set)
{
  ScanHMMSet(hmm_set, mtm_all & ~(mtm_XformInstance|mtm_Xform), NULL,
             NormalizeAccum, NULL);
}



int faddfloat(FLOAT *vec, size_t size, float mul_const, FILE *fp) {
  int i;
  FLOAT f;

  for(i = 0; i < size; i++) {
    if(fread(&f, sizeof(FLOAT), 1, fp) != 1) break;
    vec[i] += f * mul_const;
  }
  return i;
}

typedef struct {
  FILE   *fp;
  char   *fn;
  HMMSet *hmm_set;
  float  weight;
  int    mmi;
} ReadAccumUserData;

void ReadAccum(int macro_type, HMMSetNodeName nodeName,
                void *data, void *userData) {
  int i, j, c, size;
  FLOAT *vector = NULL;
  Macro *macro;
  ReadAccumUserData *ud = (ReadAccumUserData *) userData;
//  FILE *fp =        ((ReadAccumUserData *) userData)->fp;
//  char *fn =        ((ReadAccumUserData *) userData)->fn;
//  HMMSet *hmm_set = ((ReadAccumUserData *) userData)->hmm_set;
//  float weight =    ((ReadAccumUserData *) userData)->weight;
  char xfName[128];

  xfName[sizeof(xfName)-1] = '\0';

  if(macro_type == mt_mean || macro_type == mt_variance) {
    XformStatAccum *xfsa = NULL;
    int size_inf, nxfsa_inf, nxfsa = 0;

    if(macro_type == mt_mean) {
      xfsa   = ((Mean *)data)->xformStatAccum;
      nxfsa  = ((Mean *)data)->nxformStatAccums;
      size   = ((Mean *)data)->vec_size;
      vector = ((Mean *)data)->vector+size;
      size   = size + 1;
    } else if(macro_type == mt_variance) {
      xfsa   = ((Variance *)data)->xformStatAccum;
      nxfsa  = ((Variance *)data)->nxformStatAccums;
      size   = ((Variance *)data)->vec_size;
      vector = ((Variance *)data)->vector+size;
      size   = size * 2 + 1;
    }

    if(ud->mmi) vector += size;

    if(faddfloat(vector, size, ud->weight,     ud->fp) != size ||
       fread(&nxfsa_inf, sizeof(nxfsa_inf), 1, ud->fp) != 1) {
      Error("Incompatible accumulator file: '%s'", ud->fn);
    }

    if(!ud->mmi) { // MMI estimation of Xform statistics has not been implemented yet
      for(i = 0; i < nxfsa_inf; i++) {
        if(getc(ud->fp) != '"') {
          Error("Incompatible accumulator file: '%s'", ud->fn);
        }

        for(j=0; (c=getc(ud->fp)) != EOF && c != '"' && j < sizeof(xfName)-1; j++) {
          xfName[j] = c;
        }

        xfName[j] = '\0';
        if(c == EOF) {
          Error("Incompatible accumulator file: '%s'", ud->fn);
        }

        macro = FindMacro(&ud->hmm_set->Xform_hash, xfName);

        if(fread(&size_inf, sizeof(int), 1, ud->fp) != 1) {
          Error("Incompatible accumulator file: '%s'", ud->fn);
        }

        if(macro != NULL) {
          size = ((LinearXform *) macro->data)->in_size;
          size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;

          if(size != size_inf) {
            Error("Incompatible accumulator file: '%s'", ud->fn);
          }

          for(j = 0; j < nxfsa && xfsa[j].xform != macro->data; j++);
          if(j < nxfsa) {
            if(faddfloat(xfsa[j].stats, size, ud->weight, ud->fp) != size  ||
              faddfloat(&xfsa[j].norm,    1, ud->weight, ud->fp) != 1) {
              Error("Invalid accumulator file: '%s'", ud->fn);
            }
          } else {
            macro = NULL;
          }
        }

        if(macro == NULL) { // Skip Xform accumulator
          FLOAT f;
          for(j = 0; j < size_inf+1; j++) fread(&f, sizeof(f), 1, ud->fp);
        }
      }
    }
  } else if(macro_type == mt_state) {
    State *state = (State *) data;
    if(state->outPDF_kind == KID_DiagC) {
      FLOAT junk;
      for(i = 0; i < state->num_mixtures; i++) {
        if(ud->mmi == 1) {
          if(faddfloat(&state->mixture[i].weight_accum_den, 1, ud->weight, ud->fp) != 1 ||
             faddfloat(&junk,                               1, ud->weight, ud->fp) != 1) {
            Error("Incompatible accumulator file: '%s'", ud->fn);
          }
        } else if(ud->mmi == 2) {
          if(faddfloat(&junk,                               1, ud->weight, ud->fp) != 1 ||
             faddfloat(&state->mixture[i].weight_accum_den, 1, ud->weight, ud->fp) != 1) {
            Error("Incompatible accumulator file: '%s'", ud->fn);
          }
        } else {
          if(faddfloat(&state->mixture[i].weight_accum,     1, ud->weight, ud->fp) != 1 ||
             faddfloat(&state->mixture[i].weight_accum_den, 1, ud->weight, ud->fp) != 1) {
            Error("Incompatible accumulator file: '%s'", ud->fn);
          }
        }
      }
    }
  } else if(macro_type == mt_transition) {
    FLOAT f;
    size   = SQR(((Transition *) data)->nstates);
    vector =     ((Transition *) data)->matrix + size;

    for(i = 0; i < size; i++) {
      if(fread(&f, sizeof(FLOAT), 1, ud->fp) != 1) {
       Error("Incompatible accumulator file: '%s'", ud->fn);
      }
      if(!ud->mmi) { // MMI estimation of transition probabilities has not been implemented yet
        f += log(ud->weight);
        LOG_INC(vector[i], f);
      }
    }
  }
}

void ReadAccums(char *fileName, float weight, HMMSet *hmm_set,
                long *totFrames, FLOAT *totLogLike, int MMI_denominator_accums) {
  FILE *fp;
  char macroName[128];
  struct my_hsearch_data *hash;
  int i, t = 0, c, skip_accum = 0;
  long occurances;
  int mtm = mtm_prescan | mtm_state | mtm_mean | mtm_variance | mtm_transition;
  ReadAccumUserData ud;
  Macro *macro;

  macroName[sizeof(macroName)-1] = '\0';


  if(!strcmp(fileName, "-")) {
    fp = stdin;
  } else if((fp = fopen(fileName, "rb")) == NULL) {
    Error("Cannot open input accumulator file: '%s'", fileName);
  }

  if(fread(totFrames,  sizeof(long),  1, fp) != 1 ||
    fread(totLogLike, sizeof(FLOAT), 1, fp) != 1) {
    Error("Invalid accumulator file: '%s'", fileName);
  }

  *totFrames  *= weight;
  *totLogLike *= weight;

  ud.fn      = fileName;
  ud.fp      = fp;
  ud.hmm_set = hmm_set;
  ud.weight  = weight;
  ud.mmi     = MMI_denominator_accums;

  for(;;) {
    if(skip_accum) { // Skip to the begining of the next macro accumulator
      for(;;) {
        while((c = getc(fp)) != '~' && c != EOF);
        if(c == EOF) break;
        if(strchr("hsmuvt", t = c = getc(fp)) &&
          (c = getc(fp)) == ' ' && (c = getc(fp)) == '"')
          break;
        ungetc(c, fp);
      }
      if(c == EOF) break;
    } else {
      if((c = getc(fp)) == EOF) break;
      if(c != '~'       || !strchr("hsmuvt", t = getc(fp)) ||
        getc(fp) != ' ' || getc(fp) != '"') {
        Error("Incomatible accumulator file: '%s'", fileName);
      }
    }

    for(i=0; (c = getc(fp))!=EOF && c!='"' && i<sizeof(macroName)-1; i++) {
      macroName[i] = c;
    }
    macroName[i] = '\0';

    hash = t == 'h' ? &hmm_set->hmm_hash :
           t == 's' ? &hmm_set->state_hash :
           t == 'm' ? &hmm_set->mixture_hash :
           t == 'u' ? &hmm_set->mean_hash :
           t == 'v' ? &hmm_set->variance_hash :
           t == 't' ? &hmm_set->transition_hash : NULL;

    assert(hash);
    if((macro = FindMacro(hash, macroName)) == NULL) {
      skip_accum = 1;
      continue;
    }

    skip_accum = 0;
    if(fread(&occurances, sizeof(occurances), 1, fp) != 1) {
      Error("Invalid accumulator file: '%s'", fileName);
    }

    if(!MMI_denominator_accums) macro->occurances += occurances;

    switch(t) {
      case 'h': ScanHMM((HMM *)macro->data, mtm, NULL, ReadAccum, &ud);    break;
      case 's': ScanState((State *)macro->data, mtm, NULL, ReadAccum, &ud);break;
      case 'm': ScanMixture((Mixture *)macro->data,mtm,NULL,ReadAccum,&ud);break;
      case 'u': ReadAccum(mt_mean, NULL, macro->data, &ud);                break;
      case 'v': ReadAccum(mt_variance, NULL, macro->data, &ud);            break;
      case 't': ReadAccum(mt_transition, NULL, macro->data, &ud);          break;
      default:  assert(0);
    }
  }

  fclose(fp);
}

struct my_hsearch_data MakeCIPhoneHash(HMMSet *hmms)
{
  int i, nCIphns = 0;
  struct my_hsearch_data tmpHash, retHash;
  ENTRY e, *ep;

  if(!my_hcreate_r(100, &tmpHash)) Error("Insufficient memory");

  // Find CI HMMs and put them into hash
  for(i = 0; i < hmms->hmm_hash.nentries; i++) {
    Macro *macro = (Macro *) hmms->hmm_hash.entry[i]->data;
    if(strpbrk(macro->name, "+-")) continue;

    e.key  = macro->name;
    e.data = macro->data;
    if(!my_hsearch_r(e, ENTER, &ep, &tmpHash)) Error("Insufficient memory");
    nCIphns++;
  }

  // Find CD HMMs and mark corresponding CI HMMs in the hash
  for(i = 0; i < hmms->hmm_hash.nentries; i++) {
    char *ciname, chr;
    int cinlen;
    Macro *macro = (Macro *) hmms->hmm_hash.entry[i]->data;
    if(!strpbrk(macro->name, "+-")) continue;

    ciname = strrchr(macro->name, '-');
    if(ciname == NULL) ciname = macro->name;
    else ciname++;
    cinlen = strcspn(ciname, "+");
    chr = ciname[cinlen];
    ciname[cinlen] = '\0';
    e.key  = ciname;
    my_hsearch_r(e, FIND, &ep, &tmpHash);
    ciname[cinlen] = chr;

    if(ep != NULL && ep->data != 0) {
      ep->data = NULL;
      nCIphns--;
    }
  }
  if(!my_hcreate_r(nCIphns, &retHash)) Error("Insufficient memory");

  // To the new hash, collect only those CI HMMs from tmpHash not having CD versions
  for(i = 0; i < tmpHash.nentries; i++) {
    if(tmpHash.entry[i]->data != NULL) {
      HMM *hmm  = (HMM *) tmpHash.entry[i]->data;
      int isTee = hmm->transition->matrix[hmm->nstates - 1] > LOG_MIN;
      e.key  = tmpHash.entry[i]->key;
      e.data = (void *) isTee;
      if(!my_hsearch_r(e, ENTER, &ep, &retHash)) Error("Insufficient memory");
    }
  }
  my_hdestroy_r(&tmpHash, 0);
  return retHash;
}

const char *hlist_filter;

void ReadHMMList(
  HMMSet *hmms,
  const char *file_name,
  const char *in_mmf_dir,
  const char *in_mmf_ext)
{
  struct readline_data rld = {0};
  char *lhmm, *fhmm, *chptr;
  Macro *macro, *macro2;
  int line_no = 0;
  char mmfile[1024];
  FILE *fp;

  if((fp = my_fopen(file_name, "rt", hlist_filter)) == NULL) {
    Error("Cannot open file: '%s'", file_name);
  }
  while((lhmm = fhmm = readline(fp, &rld)) != NULL) {
    line_no++;
    if(getHTKstr(lhmm, &chptr)) {
      Error("%s (%s:%d)", chptr, file_name, line_no);
    }
    if(*chptr && getHTKstr(fhmm = chptr, &chptr)) {
      Error("%s (%s:%d)", chptr, file_name, line_no);
    }
    if((macro = FindMacro(&hmms->hmm_hash, fhmm)) == NULL) {
      mmfile[0] = '\0';
      if(in_mmf_dir) strcat(strcat(mmfile, in_mmf_dir), "/");
      strcat(mmfile, fhmm);
      if(in_mmf_ext) strcat(strcat(mmfile, "."), in_mmf_ext);
      ReadHMMSet(mmfile, hmms, fhmm);

      if((macro = FindMacro(&hmms->hmm_hash, fhmm)) == NULL) {
        Error("Definition of model '%s' not found in file '%s'", fhmm, mmfile);
      }
    }
    if(lhmm != fhmm) {
      current_mmf_name = NULL; // Global variable; macro will not be written to any output MMF
      macro2 = AddMacroToHMMSet('h', lhmm, hmms);
      assert(macro2 != NULL);
      if(macro2->data != NULL) {
        if(hmms_ignore_macro_redefinition == 0) {
          Error("Redefinition of HMM %s (%s:%d)", lhmm, file_name, line_no);
        } else {
          Warning("Redefinition of HMM %s (%s:%d) is ignored",
                  lhmm, file_name, line_no);
        }
      } else {
        macro2->data = macro->data;
      }
    }
  }
  if(ferror(fp) || my_fclose(fp)) {
    Error("Cannot read HMM list file %s", file_name);
  }
}

/*struct my_hsearch_data ReadHMMList(HMMSet *hmms, HMMSet *hmmsToUpdate,
                                   char *hmmListFileName)
{
  char line[1024];
  FILE *fp;
  int nlines=0;
  struct my_hsearch_data hash;
  HMM *hmm;

  if((fp = fopen(hmmListFileName, "rt")) == NULL) {
      Error("Cannot open file: '%s'", hmmListFileName);
  }

  if(!my_hcreate_r(100, &hash)) Error("Insufficient memory");

  while(fgets(line, sizeof(line), fp)) {
    char *hmmName = NULL;
    char termChar = '\0';
    int i = strlen(line);
    Macro *macro;
    ENTRY e, *ep;

    nlines++;
    if(line[i-1] != '\n' && getc(fp) != EOF) {
      Error("Line %d is too long in file: %s",
            nlines, hmmListFileName);
    }

    for(; i > 0 && isspace(line[i-1]); i--) line[i-1] = '\0';

    if(i == 0) continue;

    for(i = 0; isspace(line[i]); i++);

    if(line[i] == '\"' || line[i] == '\'') {
      termChar = line[i++];
    }

    hmmName = &line[i];
    while(line[i] != '\0' && line[i] != termChar &&
          (termChar != '\0' || !isspace(line[i]))) i++;

    if(termChar != '\0') {
      if(line[i] != termChar) {
        Error("Terminanting %c expected at line %d in file %s",
              termChar, nlines, hmmListFileName);
      }
      line[i++] = '\0';
    }

    if(line[i] != '\0') { // shell command follows
      for(; isspace(line[i]); i++) line[i] = '\0';
      if(line[i] != '\0') {
        Error("Only one term is expected at line %d in file %s",
              nlines, hmmListFileName);
      }
    }

    macro = FindMacro(&hmms->hmm_hash, hmmName);
    if(macro == NULL) {
    //!!! I should try to open file with HMM named
      Error("Undefined HMM '%s' at line %d in file %s",
            hmmName, nlines, hmmListFileName);
    }

    hmm = (HMM *) macro->data;

    if(hmmsToUpdate != NULL && hmmsToUpdate != hmms) {
      macro = FindMacro(&hmmsToUpdate->hmm_hash, hmmName);
      if(macro == NULL) {
        Error("Model '%s' not found in target HMM (%s:%d)",
              hmmName, nlines, hmmListFileName);
      }

      hmm->hmmToUpdate = (HMM *) macro->data;
    }

    e.key  = macro->name;
    e.data = (void *) hmm;

    if(!my_hsearch_r(e, ENTER, &ep, &hash)) {
      Error("Insufficient memory");
    }
  }

  fclose(fp);

  return hash;
}*/

void ReadXformList(HMMSet *hmm_set, const char *xformListFileName)
{
  char line[1024];
  FILE *fp;
  int nlines=0;

  if((fp = fopen(xformListFileName, "rt")) == NULL) {
      Error("ReadXformList: Cannot open file: '%s'", xformListFileName);
  }


  while(fgets(line, sizeof(line), fp)) {
    char *xformName, *makeXformShellCommand = NULL;
    char termChar = '\0';
    int i = strlen(line);
    Macro *macro;
    MakeXformCommand *mxfc;

    nlines++;
    if(line[i-1] != '\n' && getc(fp) != EOF) {
      Error("ReadXformList: Line %d is too long in file: %s",
            nlines, xformListFileName);
    }

    for(; i > 0 && isspace(line[i-1]); i--) line[i-1] = '\0';

    if(i == 0) continue;

    for(i = 0; isspace(line[i]); i++);

    if(line[i] == '\"' || line[i] == '\'') {
      termChar = line[i++];
    }

    xformName = &line[i];
    while(line[i] != '\0' && line[i] != termChar &&
          (termChar != '\0' || !isspace(line[i]))) i++;

    if(termChar != '\0') {
      if(line[i] != termChar) {
        Error("ReadXformList: Terminanting %c expected at line %d in file %s",
              termChar, nlines, xformListFileName);
      }
      line[i++] = '\0';
    }

    if(line[i] != '\0') { // shell command follows
      for(; isspace(line[i]); i++) line[i] = '\0';
      makeXformShellCommand = &line[i];
    }

    macro = FindMacro(&hmm_set->Xform_hash, xformName);
    if(macro == NULL) {
      Error("ReadXformList: Undefined Xform '%s' at line %d in file %s",
            xformName, nlines, xformListFileName);
    }

    hmm_set->xformToUpdate = (MakeXformCommand*)
      realloc(hmm_set->xformToUpdate,
              sizeof(MakeXformCommand) * ++hmm_set->nxformsToUpdate);

    if(hmm_set->xformToUpdate == NULL) {
      Error("ReadXformList: Insufficient memory");
    }

    mxfc = &hmm_set->xformToUpdate[hmm_set->nxformsToUpdate-1];
    mxfc->xform = (Xform *) macro->data;
    mxfc->shellCommand = NULL;

    if(makeXformShellCommand) {
      if((mxfc->shellCommand = strdup(makeXformShellCommand)) == NULL) {
        Error("ReadXformList: Insufficient memory");
      }
    }
  }
}
