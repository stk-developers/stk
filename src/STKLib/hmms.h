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

#ifndef HMMS_H
#define HMMS_H

#include "common.h"

#define __USE_GNU
#include <search.h>


enum MacroType {
  mt_hmm             = 'h',
  mt_state           = 's',
  mt_mixture         = 'm',
  mt_mean            = 'u',
  mt_variance        = 'v',
//mt_variance        = 'i',
  mt_transition      = 't',
  mt_XformInstance   = 'j',
  mt_Xform           = 'x',
};

typedef enum _MacroTypeMask {
  mtm_hmm             = 0x0001,
  mtm_state           = 0x0002,
  mtm_mixture         = 0x0004,
  mtm_mean            = 0x0008,
  mtm_variance        = 0x0010,
  mtm_transition      = 0x0020,
  mtm_XformInstance   = 0x0040,
  mtm_Xform           = 0x0080,
  mtm_all             = 0x00ff,
  mtm_revpass         = 0x4000, // Process last macro first
  mtm_prescan         = 0x8000  // Process HMMs then states then mixtures, ...
} MacroTypeMask;


typedef struct _Macro Macro;
typedef struct _HMMSet HMMSet;
typedef struct _HMM HMM;
typedef struct _State State;
typedef struct _Mixture Mixture;
typedef struct _Mean Mean;
typedef struct _Variance Variance;
typedef struct _Transition Transition;
typedef struct _XformInstance XformInstance;
typedef struct _Xform Xform;
typedef struct _CompositeXform CompositeXform;
typedef struct _LinearXform LinearXform;
typedef struct _CopyXform CopyXform;
typedef struct _BiasXform BiasXform;
typedef struct _FuncXform FuncXform;
typedef struct _StackingXform StackingXform;
typedef struct _XformStatCache XformStatCache;
typedef struct _XformStatAccum XformStatAccum;


struct _Macro {
  char  *name;
  char  *file;
  void  *data;
  long  occurances;
//  Macro *next;
  int   type;
  Macro *nextAll;
  Macro *prevAll;
};

typedef enum {
  UM_TRANSITION = 1,
  UM_MEAN       = 2,
  UM_VARIANCE   = 4,
  UM_WEIGHT     = 8,
  UM_OLDMEANVAR = 16,
  UM_XFSTATS    = 32,
  UM_XFORM      = 64
} UpdateMask;

typedef enum {
  KID_BeginHMM,    KID_Use,        KID_EndHMM,    KID_NumMixes,  KID_NumStates,
  KID_StreamInfo,  KID_VecSize,    KID_NullD,     KID_PoissonD,  KID_GammaD,
  KID_RelD,        KID_GenD,       KID_DiagC,     KID_FullC,     KID_XformC,
  KID_State,       KID_TMix,       KID_Mixture,   KID_Stream,    KID_SWeights,
  KID_Mean,        KID_Variance,   KID_InvCovar,  KID_Xform,     KID_GConst,
  KID_Duration,    KID_InvDiagC,   KID_TransP,    KID_DProb,     KID_LLTC,
  KID_LLTCovar,
  KID_XformKind=90,KID_ParentXform,KID_NumXforms, KID_XformSet,  KID_LinXform,
  KID_Offset,      KID_Bias,       KID_BlockInfo, KID_Block,     KID_BaseClass,
  KID_Class,       KID_XformWgtSet,KID_ClassXform,KID_MMFIDMask, KID_Parameters,
  KID_NumClasses,  KID_AdaptKind,  KID_Prequal,   KID_InputXform,
  KID_RClass  =110,KID_RegTree,    KID_Node,      KID_TNode,
  KID_HMMSetID=119,KID_ParmKind,

  /* Non-HTK keywords */
  KID_FrmExt  =200,KID_PDFObsVec,  KID_ObsCoef,    KID_Input,    KID_NumLayers,
  KID_NumBlocks,   KID_Layer,      KID_Copy,       KID_Stacking,

  /* Numeric functions - FuncXform*/
  KID_Sigmoid,     KID_Log,        KID_Exp,        KID_Sqrt,     KID_SoftMax,

  KID_MaxKwdID
} KeywordID;

extern char *Kwds[KID_MaxKwdID];

typedef struct _MakeXformCommand {
  Xform *xform;
  char  *shellCommand;
} MakeXformCommand;

#define DEFAULT_XFORM_NAME "defaultInputXform"

struct _HMMSet {
  Macro *first_macro;
  Macro *last_macro;
  struct my_hsearch_data hmm_hash;
  struct my_hsearch_data state_hash;
  struct my_hsearch_data mixture_hash;
  struct my_hsearch_data mean_hash;
  struct my_hsearch_data variance_hash;
  struct my_hsearch_data transition_hash;
  struct my_hsearch_data Xform_instance_hash;
  struct my_hsearch_data Xform_hash;
  int   in_vec_size;
  int   param_kind;
  long  nmixtures;
  long  nstates;
  int   alloc_accums;
  int   totalDelay;
  int   allMixuresUpdatableFromStatAccums;
  int   isHTKCopatible; // Models use no extension with respecto to HTK
  KeywordID outPDF_kind;
  KeywordID dur_kind;
  XformInstance *inputXform;
//XformInstance *linXform;
  XformInstance *Xform_instances;
  
  //Reestimation params
  UpdateMask updateMask;
  FLOAT minMixWeight;
  Variance *varFloor;
  long minOccurances;
  MakeXformCommand *xformToUpdate;
  int nxformsToUpdate;
  int gaussLvl2ModelReest;
  int MMIUpdate;
  FLOAT MMI_E;
  FLOAT MMI_h;
  FLOAT MMI_tauI;
};

struct _HMM {
  Macro      *macro;
  int        nstates;
  Transition *transition;
  State      *state[1];
};

struct _State {
  Macro      *macro;
  long       state_id;

  KeywordID outPDF_kind;
  union {
    int num_mixtures;
    int PDF_obs_coef;
  };

  struct {
    Mixture *estimates;
    FLOAT   weight;
    FLOAT   weight_accum; //used for reestimation
    FLOAT   weight_accum_den;
  } mixture[1];
};

struct _Mixture {
  Macro    *macro;
  long     mixture_id;
  Mean     *mean;
  Variance *variance;
  FLOAT    g_const;
  XformInstance *inputXform;
//  XformInstance *linXform;
//  FLOAT *linXformStats;
};

struct _XformStatAccum {
  Xform *xform;
  FLOAT norm;
  FLOAT *stats;
};

struct _Mean {
  Macro          *macro;
  int            vec_size;
  XformStatAccum *xformStatAccum;
  int            nxformStatAccums;
  int            updatableFromStatAccums;
  FLOAT          vector[1];
};

struct _Variance {
  Macro          *macro;
//  BOOL         diagonal;
  int            vec_size;
  XformStatAccum *xformStatAccum;
  int            nxformStatAccums;
  int            updatableFromStatAccums;
  FLOAT          vector[1];
};

struct _Transition {
  Macro *macro;
  int   nstates;
  FLOAT matrix[1];
};

struct _XformStatCache {
  XformStatCache *upperLevelStats;
  Xform          *xform;
  int            norm;
  FLOAT          *stats;
};

struct _XformInstance {
  Macro *macro;
  XformInstance *input;
  Xform *xform;
  int time;
  XformInstance  *next; // Chain of all instances
  XformStatCache *xformStatCache;
  int nxformStatCaches;
  int out_size;
  int statCacheTime;
  char *memory;
  int totalDelay;
  
  FLOAT out_vec[1]; //stackSize * (xform ? xform->out_size : hmm_set->in_vec_size)
};

typedef enum {
  XT_LINEAR,
  XT_COPY,
  XT_BIAS,
  XT_FUNC,
  XT_STACKING,
  XT_COMPOSITE,
} XformType;


struct _Xform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;
  };

struct _CompositeXform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;

  int   nlayers;
  struct {
    FLOAT *out_vec;
    int   nblocks;
    Xform **block;
  } layer[1];
};

struct _LinearXform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;

  FLOAT matrix[1];
};

struct _BiasXform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;

  FLOAT vector[1];
};

struct _FuncXform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;

  int funcId;
};

struct _CopyXform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;

  int indices[1];
};


struct _StackingXform {
  Macro *macro;
  XformType xform_type;
  int in_size;
  int out_size;
  int memorySize;
  int delay;  
  
  int horiz_stack;
};


#ifdef __cplusplus
  extern "C" {
#endif

void InitHMMSet(HMMSet *hmm_set, int reest);
void ReadHMMSet(const char *mmFileName, HMMSet *hmm_set, char *expectHMM);
void WriteHMMSet(const char *mmfName, const char *out_mmf_dir,
                 const char *out_mmf_ext, int binary, HMMSet *hmm_set);
void ReadAccums(char *fileName, float weight, HMMSet *hmm_set,
                long *totFrames, FLOAT *totLogLike, int MMI_denominator_accums);
void WriteAccums(const char *accfName, const char *out_dir, HMMSet *hmm_set,
                 long totFrames, FLOAT totLogLike);
void NormalizeAccums(HMMSet *hmm_set);
void ReleaseHMMSet(HMMSet *hmm_set);
void ResetAccumsForHMMSet(HMMSet *hmm_set);
void GlobalStatsForHMMSet(HMMSet *hmm_set, FLOAT *observation, int time);
void UpdateHMMSetFromAccums(const char *out_dir, HMMSet *hmm_set);
void DistributeMacroOccurances(HMMSet *hmm_set);
typedef char HMMSetNodeName[128];

typedef void (*ScanAction)(int type, HMMSetNodeName nodeName,
                           void *data, void *userData);

void ScanHMMSet(HMMSet *hmm_set, MacroTypeMask mask,
                HMMSetNodeName nodeNameBuffer,
                ScanAction action, void *userData);

void ScanHMM(HMM *hmm, MacroTypeMask mask, HMMSetNodeName nodeName,
             ScanAction action, void *userData);

void ScanState(State *state, MacroTypeMask mask, HMMSetNodeName nodeName,
               ScanAction action, void *userData);

void ScanMixture(Mixture *mixture, MacroTypeMask mask,
                 HMMSetNodeName nodeName, ScanAction action, void *userData);

void ScanXformInstance(XformInstance *xformInstance, MacroTypeMask mask,
                       HMMSetNodeName nodeName, ScanAction action,
                       void *userData);

void ScanXform(Xform *xform, MacroTypeMask mask, HMMSetNodeName nodeName,
               ScanAction action, void *userData);

void WriteHMMStats(const char *stat_file, HMMSet *hmm_set);
void WriteXformStatsAndRunCommands(const char *out_dir, int binary, HMMSet *hmm_set);
void ReadXformStats(const char *out_dir, int binary, HMMSet *hmm_set);

/*void UpdateXformsAndModels(HMMSet *hmm_set,
                           const char *out_mmf_dir,
                           char *make_xform_command, int binary);*/


Macro *FindMacro(struct my_hsearch_data *macro_hash, const char *name);
FLOAT *XformPass(XformInstance *xformInstance, FLOAT *in_vec, int time, PropagDir dir);
void ResetXformInstances(HMMSet *hmm_set);
void AllocateAccumulatorsForXformStats(HMMSet *hmm_set);
void UpdateStacks(HMMSet *hmm_set, FLOAT *obs, int time,  PropagDir dir);
struct my_hsearch_data MakeCIPhoneHash(HMMSet *hmms);
//struct my_hsearch_data ReadHMMList(HMMSet *hmms, HMMSet *hmmsToUpdate,
//                                   char *hmmListFileName);

void ReadHMMList(
  HMMSet *hmms,
  const char *file_name,
  const char *in_mmf_dir,
  const char *in_mmf_ext);

void ReadXformList(HMMSet *hmm_set, const char *xformListFileName);

void NormalizeStatsForXform(int macro_type, HMMSetNodeName nodeName,
                            void *data, void *userData);
#ifdef __cplusplus
}
#endif

extern int hmms_ignore_macro_redefinition;
extern const char *hlist_filter;


#endif // HMMS_H
