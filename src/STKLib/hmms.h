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

  
#define DEFAULT_XFORM_NAME "defaultInputXForm"  
  
  class MacroHash;  
  class MacroData;
  class Macro;
  class Hmm;  
  class State;
  class Mixture;
  class MixtureLink;
  class Mean;
  class Variance;
  class Transition;
  class XFormInstance;
  class XForm;
  class CompositeXForm;
  class CopyXForm;
  class LinearXForm;
  class BiasXForm;
  class FuncXForm;
  class StackingXForm;
  class XFormStatCache;
  class XFormStatAccum;

  enum MacroType 
  {
    mt_hmm             = 'h',
    mt_state           = 's',
    mt_mixture         = 'm',
    mt_mean            = 'u',
    mt_variance        = 'v',
  //mt_variance        = 'i',
    mt_transition      = 't',
    mt_XFormInstance   = 'j',
    mt_XForm           = 'x',
  };

  typedef enum _MacroTypeMask 
  {
    mtm_hmm             = 0x0001,
    mtm_state           = 0x0002,
    mtm_mixture         = 0x0004,
    mtm_mean            = 0x0008,
    mtm_variance        = 0x0010,
    mtm_transition      = 0x0020,
    mtm_XFormInstance   = 0x0040,
    mtm_XForm           = 0x0080,
    mtm_all             = 0x00ff,
    mtm_revpass         = 0x4000, // Process last macro first
    mtm_prescan         = 0x8000  // Process HMMs then states then mixtures, ...
  } MacroTypeMask;

  
  /// model set flags
  /// @{  
  const FlagType MODEL_SET_WITH_ACCUM   = 1; ///< The set should allocate space for
                                             ///< accumulators 
  /// @}
  

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Macro data base class
   */
  class MacroData
  {
  public:
    Macro *                     mpMacro;
  
  private:
    long int                    mLinksCount;
    
  public:  
    /// Default constructor
    MacroData() : mpMacro(NULL), mLinksCount(0) {}
    
    const long int
    Links() const
    {
      return mLinksCount;
    }    
    
    // we'll make the Macro class a frien as it will increment the mLinksCount
    friend class Macro;
  }; // class MacroData

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Macro representation
   */
  class Macro 
  {
  public:
    char *        mpName;
    char *        mpFileName;
    MacroData *   mpData;
    long          mOccurances;
    int           type;
    Macro *       nextAll;
    Macro *       prevAll;
  };
  
          
  typedef enum 
  {
    UM_TRANSITION = 1,
    UM_MEAN       = 2,
    UM_VARIANCE   = 4,
    UM_WEIGHT     = 8,
    UM_OLDMEANVAR = 16,
    UM_XFSTATS    = 32,
    UM_XFORM      = 64
  } UpdateMask;
  
  typedef enum {
    KID_UNSET=-1,
    KID_BeginHMM,    KID_Use,        KID_EndHMM,    KID_NumMixes,  KID_NumStates,
    KID_StreamInfo,  KID_VecSize,    KID_NullD,     KID_PoissonD,  KID_GammaD,
    KID_RelD,        KID_GenD,       KID_DiagC,     KID_FullC,     KID_XFormC,
    KID_State,       KID_TMix,       KID_Mixture,   KID_Stream,    KID_SWeights,
    KID_Mean,        KID_Variance,   KID_InvCovar,  KID_XForm,     KID_GConst,
    KID_Duration,    KID_InvDiagC,   KID_TransP,    KID_DProb,     KID_LLTC,
    KID_LLTCovar,
    KID_XFormKind=90,KID_ParentXForm,KID_NumXForms, KID_XFormSet,  KID_LinXForm,
    KID_Offset,      KID_Bias,       KID_BlockInfo, KID_Block,     KID_BaseClass,
    KID_Class,       KID_XFormWgtSet,KID_ClassXForm,KID_MMFIDMask, KID_Parameters,
    KID_NumClasses,  KID_AdaptKind,  KID_Prequal,   KID_InputXForm,
    KID_RClass  =110,KID_RegTree,    KID_Node,      KID_TNode,
    KID_HMMSetID=119,KID_ParmKind,
  
    /* Non-HTK keywords */
    KID_FrmExt  =200,KID_PDFObsVec,  KID_ObsCoef,    KID_Input,    KID_NumLayers,
    KID_NumBlocks,   KID_Layer,      KID_Copy,       KID_Stacking,
  
    /* Numeric functions - FuncXForm*/
    KID_Sigmoid,     KID_Log,        KID_Exp,        KID_Sqrt,     KID_SoftMax,
  
    KID_MaxKwdID
  } KeywordID;
  
  extern char *Kwds[KID_MaxKwdID];
  
  
  class MakeXFormCommand 
  {
  public:
    XForm *xform;
    char  *shellCommand;
  };


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Master model class
   * 
   *  This class encapsulates STK model manipulation
   */
  class ModelSet 
  {
  public:
    Macro *                   first_macro;
    Macro *                   last_macro;
    
    struct my_hsearch_data    mHmmHash;
    struct my_hsearch_data    mStateHash;
    struct my_hsearch_data    mMixtureHash;
    struct my_hsearch_data    mMeanHash;
    struct my_hsearch_data    mVarianceHash;
    struct my_hsearch_data    mTransitionHash;
    struct my_hsearch_data    mXFormInstanceHash;
    struct my_hsearch_data    mXFormHash;
    
    int                       mInputVectorSize;
    int                       mParamKind;
    long                      mNMixtures;
    long                      mNStates;
    bool                      mAllocAccums;
    int                       mTotalDelay;
    int                       allMixuresUpdatableFromStatAccums;
    bool                      isHTKCopatible; // Models use no extension with respecto to HTK
    
    KeywordID                 outPDF_kind;
    KeywordID                 dur_kind;
    XFormInstance *           mpInputXForm;
  //XFormInstance *           linXForm;
    XFormInstance *           mpXFormInstances;
    
    //Reestimation params
    int                       updateMask;
    FLOAT                     minMixWeight;
    Variance *                varFloor;
    long                      minOccurances;
    MakeXFormCommand *        xformToUpdate;
    int                       nxformsToUpdate;
    int                       gaussLvl2ModelReest;
    int                       MMIUpdate;
    FLOAT                     MMI_E;
    FLOAT                     MMI_h;
    FLOAT                     MMI_tauI;
    
    
    /**
     *  @brief Reads definition from the stream and parses
     *  @param rName filename to parse
     *  @return @c true on success
     *  
     *  The entire stream is read and propriate structures are constructed.
     */  
    void
    ParseMmf(const std::string & rName, char * expectHMM);
    
    
    /**
     *  @brief Initializes the Model set
     *  @param modelSetType tells whether the models should allocate space for accumulators
     */  
    void
    Init(FlagType flags = 0);
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Basic class defining HMM 
   * 
   *  This class represents the Hidden Markov Model structure and defines
   *  supporting procedures.
   */
  class Hmm : public MacroData 
  {
  public:
    size_t                    mNStates;
    Transition *              mpTransition;
    State **                  mpState;
    State      *              state[1];
    
  public:
    /// Constructor
    Hmm(size_t nStates);
    
    ~Hmm();
  };
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Defines HMM State representation
   */
  class State : public MacroData  
  {
  public:
    long                      mID;
  
    KeywordID outPDF_kind;
    union 
    {
      unsigned int            num_mixtures;
      int                     PDF_obs_coef;
    };
  
    struct 
    {
    public:
      Mixture *               estimates;
      FLOAT                   weight;
      FLOAT                   weight_accum; //used for reestimation
      FLOAT                   weight_accum_den;
    } mixture[1];
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Defines HMM Mixture
   *
   */
  class Mixture : public MacroData 
  {
  public:
    long                      mID;
    Mean *                    mpMean;
    Variance *                mpVariance;
    FLOAT                     mGConst;
    XFormInstance *           mpInputXForm;
  
    /// Returns the GConst constant
    const FLOAT
    GConst() const
    {
      return mGConst;
    }
    
    /// Computes the GConst constant
    void
    ComputeGConst();
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief XForm statistics accumulator
   */
  class XFormStatAccum
  {
  public:
    XForm *                 xform;
    FLOAT                   norm;
    FLOAT *                 stats;
  };
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Mean representation
   */
  class Mean : public MacroData 
  {
  public:
    int                     vec_size;
    XFormStatAccum *        xformStatAccum;
    int                     nxformStatAccums;
    int                     updatableFromStatAccums;
    FLOAT                   vector[1];
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Variance representation
   */
  class Variance : public MacroData 
  {
  public:
  //  BOOL         diagonal;
    int                     vec_size;
    XFormStatAccum *        xformStatAccum;
    int                     nxformStatAccums;
    int                     updatableFromStatAccums;
    FLOAT                   vector[1];
  };
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Defines HMM Transition representation (matrix)
   *
   *
   */
  class Transition : public MacroData 
  {
  public:
    int                     mNStates;
    FLOAT                   matrix[1];
  };


  class XFormStatCache {
  public:
    XFormStatCache *upperLevelStats;
    XForm          *xform;
    int            norm;
    FLOAT          *stats;
  };

  class XFormInstance : public MacroData 
  {
  public:
    XFormInstance *       input;
    XForm *               xform;
    int                   time;
    XFormInstance  *      next; // Chain of all instances
    XFormStatCache *      xformStatCache;
    int                   nxformStatCaches;
    int                   out_size;
    int                   statCacheTime;
    char *                memory;
    int                   mTotalDelay;
    
    FLOAT                 out_vec[1]; 
    //stackSize * (xform ? xform->out_size : hmm_set->mInputVectorSize)
  };

  typedef enum {
    XT_LINEAR,
    XT_COPY,
    XT_BIAS,
    XT_FUNC,
    XT_STACKING,
    XT_COMPOSITE,
  } XFormType;
  

  class XForm : public MacroData 
  {
  public:
    XFormType xform_type;
    int in_size;
    int out_size;
    int memorySize;
    int delay;
  };

  class CompositeXForm : public XForm 
  {
  public:
    int   nlayers;
    
    struct
    {
      FLOAT *out_vec;
      int   nblocks;
      XForm **block;
    } layer[1];
  };
  
  class LinearXForm : public XForm 
  {
  public:
    FLOAT matrix[1];
  };
  
  class BiasXForm : public XForm 
  {
  public:
    FLOAT vector[1];
  };
  
  class FuncXForm : public XForm 
  {
  public:
    int funcId;
  };
  
  class CopyXForm : public XForm 
  {
  public:
    int indices[1];
  };
  
  
  class StackingXForm : public XForm 
  {
  public:
    int horiz_stack;
  };
  


void initKwdTab();
void ReadHMMSet(const char *mmFileName, ModelSet *hmm_set, char *expectHMM);
void WriteHMMSet(const char *mmfName, const char *out_mmf_dir,
                 const char *out_mmf_ext, int binary, ModelSet *hmm_set);
void ReadAccums(char *fileName, float weight, ModelSet *hmm_set,
                long *totFrames, FLOAT *totLogLike, int MMI_denominator_accums);
void WriteAccums(const char *accfName, const char *out_dir, ModelSet *hmm_set,
                 long totFrames, FLOAT totLogLike);
void NormalizeAccums(ModelSet *hmm_set);
void ReleaseHMMSet(ModelSet *hmm_set);
void ResetAccumsForHMMSet(ModelSet *hmm_set);
void GlobalStatsForHMMSet(ModelSet *hmm_set, FLOAT *observation, int time);
void UpdateHMMSetFromAccums(const char *out_dir, ModelSet *hmm_set);
void DistributeMacroOccurances(ModelSet *hmm_set);
typedef char HMMSetNodeName[128];

typedef void (*ScanAction)(int type, HMMSetNodeName nodeName,
                           void *data, void *userData);

void ScanHMMSet(ModelSet *hmm_set, int mask,
                HMMSetNodeName nodeNameBuffer,
                ScanAction action, void *userData);

void ScanHMM(Hmm *hmm, int mask, HMMSetNodeName nodeName,
             ScanAction action, void *userData);

void ScanState(State *state, int mask, HMMSetNodeName nodeName,
               ScanAction action, void *userData);

void ScanMixture(Mixture *mixture, int mask,
                 HMMSetNodeName nodeName, ScanAction action, void *userData);

void ScanXFormInstance(XFormInstance *xformInstance, int mask,
                       HMMSetNodeName nodeName, ScanAction action,
                       void *userData);

void ScanXForm(XForm *xform, int mask, HMMSetNodeName nodeName,
               ScanAction action, void *userData);

void WriteHMMStats(const char *stat_file, ModelSet *hmm_set);
void WriteXFormStatsAndRunCommands(const char *out_dir, int binary, ModelSet *hmm_set);
void ReadXFormStats(const char *out_dir, int binary, ModelSet *hmm_set);

/*void UpdateXFormsAndModels(ModelSet *hmm_set,
                           const char *out_mmf_dir,
                           char *make_xform_command, int binary);*/


Macro *FindMacro(struct my_hsearch_data *macro_hash, const char *name);
FLOAT *XFormPass(XFormInstance *xformInstance, FLOAT *in_vec, int time, PropagDir dir);
void ResetXFormInstances(ModelSet *hmm_set);
void AllocateAccumulatorsForXFormStats(ModelSet *hmm_set);
void UpdateStacks(ModelSet *hmm_set, FLOAT *obs, int time,  PropagDir dir);
struct my_hsearch_data MakeCIPhoneHash(ModelSet *hmms);
//struct my_hsearch_data ReadHMMList(ModelSet *hmms, ModelSet *hmmsToUpdate,
//                                   char *hmmListFileName);

void ReadHMMList(
  ModelSet *hmms,
  const char *file_name,
  const char *in_mmf_dir,
  const char *in_mmf_ext);

void ReadXFormList(ModelSet *hmm_set, const char *xformListFileName);

void NormalizeStatsForXForm(int macro_type, HMMSetNodeName nodeName,
                            void *data, void *userData);

extern int hmms_ignore_macro_redefinition;
extern const char *hlist_filter;


#endif // HMMS_H
