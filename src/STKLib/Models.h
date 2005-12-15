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

#ifndef STK_Models_h
#define STK_Models_h

#include "common.h"
#include "stkstream.h"
#include <search.h>


#define SQR(x) ((x) * (x))

#define WARN_FEW_EXAMPLES(type, name, exs) \
  Warning(type" %s is not updated (%s%ld example%s)", \
          name, exs == 0 ? "" : "only ", exs, exs == 1 ? "" : "s")
  
#define DEFAULT_XFORM_NAME "defaultInputXForm"  



namespace STK
{  
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

  
  /**
   * @name MacroTypeMask
   * These constants define the flags used by Scan(...) methods. The method 
   * uses these flags to decide which macros are processed by the
   * action procedure.
   */
  //@{
  const int MTM_HMM             = 0x0001;
  const int MTM_STATE           = 0x0002;
  const int MTM_MIXTURE         = 0x0004;
  const int MTM_MEAN            = 0x0008;
  const int MTM_VARIANCE        = 0x0010;
  const int MTM_TRANSITION      = 0x0020;
  const int MTM_XFORM_INSTANCE  = 0x0040;
  const int MTM_XFORM           = 0x0080;
  const int MTM_ALL             = 0x00ff;
  const int MTM_REVERSE_PASS    = 0x4000; ///< Process last macro first
  const int MTM_PRESCAN         = 0x8000; ///< Process HMMs then states then mixtures, ...
  //@}

  
  typedef char HMMSetNodeName[128];

  /**
   * @brief Scan action procedure type
   * 
   * This function type is used by the scan functions to perform the desired
   * action on them.
   */  
  typedef void (*ScanAction)( int             type,           ///< Type of data to perform the action on
                              HMMSetNodeName  nodeName,       ///< node name
                              MacroData *     pData,          ///< Macro data class 
                              void *          pUserData);     ///< User data to process

  
  /// model set flags
  //@{  
  const FlagType MODEL_SET_WITH_ACCUM   = 1; ///< The set should allocate space for
                                             ///< accumulators 
  //@}
  
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief MMF file input operation encapsulation
   */
   class MmfReader
   {
   private:
      istkstream        mStream;   ///< stream to read from
      int               mCurLine;  ///< current line
      int               mCurCol;   ///< current column
      
   public:      
      /**
       * @name File operation functions
       */
      //@{
      FLOAT       GetFloat();
      int         GetInt();
      char *      GetString(int eofNotExpected);
      void        SkipSpaces();
      void        UngetString();
      //@}
      
      const std::string
      FileName() const
      {
        return mStream.name();
      }
      
   };
  
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Macro data base class
   */
  class MacroData
  {
  public:
    Macro *                     mpMacro;
  
  public:  
    /// Default constructor
    MacroData() : mpMacro(NULL) {}
    
    // we'll make the Macro class a frien as it will increment the mLinksCount
    friend class Macro;
    
    virtual void
    Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *userData) = 0;
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
    int           mType;
    Macro *       nextAll;
    Macro *       prevAll;
  };
  
          
  enum UpdateMask
  {
    UM_TRANSITION = 1,
    UM_MEAN       = 2,
    UM_VARIANCE   = 4,
    UM_WEIGHT     = 8,
    UM_OLDMEANVAR = 16,
    UM_XFSTATS    = 32,
    UM_XFORM      = 64
  };
  
  enum KeywordID
  {
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
  };
 
  
  class MakeXFormCommand 
  {
  public:
    XForm *mpXForm;
    char  *mpShellCommand;
  };
  
  struct FunctionTable
  {
    void (*funcPtr)(FLOAT *, FLOAT*, int);
    KeywordID KID;
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

  private:
  
    
    /// This group of methods does exactly what their names state
    /// @{ 
    Hmm *           ReadHMM           (FILE * fp, Macro * macro);
    State *         ReadState         (FILE * fp, Macro * macro);
    Mixture *       ReadMixture       (FILE * fp, Macro * macro);
    Mean *          ReadMean          (FILE * fp, Macro * macro);
    Variance *      ReadVariance      (FILE * fp, Macro * macro);
    Transition *    ReadTransition    (FILE * fp, Macro * macro);
    XFormInstance * ReadXFormInstance (FILE * fp, Macro * macro);
    XForm *         ReadXForm         (FILE * fp, Macro * macro);
    CompositeXForm *ReadCompositeXForm(FILE * fp, Macro * macro);
    LinearXForm *   ReadLinearXForm   (FILE * fp, Macro * macro);
    CopyXForm *     ReadCopyXForm     (FILE * fp, Macro * macro);
    BiasXForm *     ReadBiasXForm     (FILE * fp, Macro * macro);
    FuncXForm *     ReadFuncXForm     (FILE * fp, Macro * macro, int funcId);
    StackingXForm * ReadStackingXForm (FILE * fp, Macro * macro);
    int             ReadGlobalOptions (FILE * fp);
    
    void   WriteHMM          (FILE * fp, bool binary, Hmm * hmm);
    void   WriteState        (FILE * fp, bool binary, State * state);
    void   WriteMixture      (FILE * fp, bool binary, Mixture * mixture);
    void   WriteMean         (FILE * fp, bool binary, Mean * mean);
    void   WriteVariance     (FILE * fp, bool binary, Variance * variance);
    void   WriteTransition   (FILE * fp, bool binary, Transition * transition);
    void   WriteXFormInstance(FILE * fp, bool binary, XFormInstance * xformInstance);
    void   WriteXForm        (FILE * fp, bool binary,          XForm * xform);
    void   WriteCompositeXForm(FILE *fp, bool binary, CompositeXForm * xform);
    void   WriteLinearXForm  (FILE * fp, bool binary,    LinearXForm * xform);
    void   WriteCopyXForm    (FILE * fp, bool binary,      CopyXForm * xform);
    void   WriteFuncXForm    (FILE * fp, bool binary,      FuncXForm * xform);
    void   WriteBiasXForm    (FILE * fp, bool binary,      BiasXForm * xform);
    void   WriteStackingXForm(FILE * fp, bool binary,  StackingXForm * xform);
    void   WriteGlobalOptions(FILE * fp, bool binary);  
    /// @}
    
    
  public:
    Macro *                   mpFirstMacro;
    Macro *                   mpLastMacro;
    
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
    size_t                    mNStates;
    bool                      mAllocAccums;
    int                       mTotalDelay;
    bool                      mAllMixuresUpdatableFromStatAccums;
    bool                      mIsHTKCopatible;                       ///< Models use no extension with respecto to HTK
    
    KeywordID                 mOutPdfKind;
    KeywordID                 mDurKind;
    XFormInstance *           mpInputXForm;
    XFormInstance *           mpXFormInstances;
    
    //Reestimation params
    int                       mUpdateMask;
    FLOAT                     mMinMixWeight;
    Variance *                mpVarFloor;
    long                      mMinOccurances;
    MakeXFormCommand *        mpXFormToUpdate;
    size_t                    mNumberOfXFormsToUpdate;
    int                       mGaussLvl2ModelReest;
    int                       mMmiUpdate;
    FLOAT                     MMI_E;
    FLOAT                     MMI_h;
    FLOAT                     MMI_tauI;
    
    
    /**
     * @name MMF output functions
     */
    //@{
    void 
    /**
     * @brief Writes the complete Model set to a file
     * @param rFileName file name to write to
     * @param rOutputDir output path
     * @param rOutputExt output file's implicit extension
     * @param binary write in binary mode if true
     */
    WriteMmf(const std::string & rFileName, 
             const std::string & rOutputDir,
             const std::string & rOutputExt, bool binary);
    
    void 
    /**
     * @brief Writes HMM statistics to a file
     * @param rFileName file to write to
     */
    WriteHMMStats(const std::string & rFileName);             
    
    
    /**
     * @brief Writes accumulators to a file
     * @param rFileName 
     * @param rOutputDir 
     * @param totFrames 
     * @param totLogLike 
     */
    void 
    WriteAccums(const std::string & rFileName, 
                const std::string & rOutputDir,
                long                totFrames, 
                FLOAT               totLogLike);
    
    void 
    /**
     * 
     * @param rOutDir 
     * @param binary 
     */
    WriteXFormStatsAndRunCommands(const std::string & rOutDir, bool binary);
    //@}
    
    
    /**
     *  @brief Initializes the Model set
     *  @param modelSetType tells whether the models should allocate space for accumulators
     */  
    void
    Init(FlagType flags = 0);
    
    /**
     *  @brief Releases memory occupied by this object's structures
     */
    void
    Release();
      
    
    /**
     * @name MMF input functions
     */
    //@{
    void
    /**
     *  @brief Reads definition from the stream and parses
     *  @param rName filename to parse
     *  @return @c true on success
     *  
     *  The entire stream is read and propriate structures are constructed.
     */  
    ParseMmf(const std::string & rName, char * expectHMM);
    
    void 
    /**
     * @brief Loads an HMM list from file
     * @param rFileName  file to open
     *
     */
    ReadHMMList(const std::string & rFileName, 
                const std::string & rInMmfDir, 
                const std::string & rInMmfExt);
    
    /**
     * @brief Reads accumulators from file
     * @param rFileName file to read from
     * @param weight 
     * @param totFrames 
     * @param totLogLike 
     * @param MMI_denominator_accums 
     */
    void
    ReadAccums(const std::string & rFileName, 
               float               weight,
               long *              totFrames, 
               FLOAT *             totLogLike, 
               int                 mmiDenominatorAccums);
    
    
    void
    ReadXFormStats(const std::string & rOutDir, bool binary);
    
    void
    ReadXFormList(const std::string & rFileName);
    //@}
    
    
    void
    AllocateAccumulatorsForXFormStats();
                    
    void 
    NormalizeAccums();
    
    /**
     * @brief Resets the accumulators to 0
     */
    void 
    ResetAccums();
    
    void 
    ComputeGlobalStats(FLOAT *observation, int time);
    
    void 
    UpdateFromAccums(const std::string & rOutputDir);
    
    void 
    DistributeMacroOccurances();
    
    void 
    ResetXFormInstances();
  
    void
    UpdateStacks(FLOAT *obs, int time,  PropagDir dir);
    
    Macro *
    pAddMacro(const char type, const std::string & rNewName);    
    
    struct my_hsearch_data 
    MakeCIPhoneHash();
    
    
    void 
    /**
     * @brief Performs desired @c action on the set's data which are chosen by @mask
     * @param mask bit mask of flags that tells which data to process by @c action
     * @param nodeNameBuffer 
     * @param action 
     * @param pUserData 
     */
    Scan( int             mask,
          HMMSetNodeName  nodeNameBuffer,
          ScanAction      action, 
          void *          pUserData);
  }; // ModelSet

  
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
    State *                   state[1];
    
  public:
    /// Constructor
    Hmm(size_t nStates);
    
    /// Destructor
    ~Hmm();
    
    void
    /**
     * @brief Updates the object from the accumulators
     * @param rModelSet ModelSet object which holds the accumulator configuration
     */
    UpdateFromAccums(const ModelSet * pModelSet);
    
    
    void
    /**
     * @brief Performs desired @c action on the HMM's data which are chosen by @mask
     * @param mask bit mask of flags that tells which data to process by @c action
     * @param nodeNameBuffer 
     * @param action routine to apply on pUserData
     * @param pUserData the data to process by @c action
     */
    Scan( int             mask,
          HMMSetNodeName  nodeNameBuffer,
          ScanAction      action, 
          void *          pUserData);
  };
  
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Defines HMM State representation
   */
  class State : public MacroData  
  {
  public:
    long                      mID;
  
    KeywordID                 mOutPdfKind;
    union 
    {
      size_t                  mNumberOfMixtures;
      int                     PDF_obs_coef;
    };
  
    struct 
    {
    public:
      Mixture *               estimates;
      FLOAT                   weight;
      FLOAT                   weight_accum; //used for reestimation
      FLOAT                   weight_accum_den;
    } mpMixture[1];
    
    void
    /**
     * @brief Updates the object from the accumulators
     * @param rModelSet ModelSet object which holds the accumulator configuration
     * @param rHmm parrent Hmm of this state
     */
    UpdateFromAccums(const ModelSet * pModelSet, const Hmm * pHmm);
  
  
    void
    /**
     * @brief Performs desired @c action on the HMM's data which are chosen by @mask
     * @param mask bit mask of flags that tells which data to process by @c action
     * @param nodeNameBuffer 
     * @param action routine to apply on pUserData
     * @param pUserData the data to process by @c action
     */
    Scan( int             mask,
          HMMSetNodeName  nodeNameBuffer,
          ScanAction      action, 
          void *          pUserData);
  
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
    
    void
    /**
     * @brief Updates the object from the accumulators
     * @param rModelSet ModelSet object which holds the accumulator configuration
     */
    UpdateFromAccums(const ModelSet * pModelSet);
    
    void
    /**
     * @brief Performs desired @c action on the HMM's data which are chosen by @mask
     * @param mask bit mask of flags that tells which data to process by @c action
     * @param nodeNameBuffer 
     * @param action routine to apply on pUserData
     * @param pUserData the data to process by @c action
     */
    Scan( int             mask,
          HMMSetNodeName  nodeNameBuffer,
          ScanAction      action, 
          void *          pUserData);
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief XForm statistics accumulator
   */
  class XFormStatAccum
  {
  public:
    XForm *                 mpXForm;
    FLOAT                   norm;
    FLOAT *                 mpStats;
  };
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Mean representation
   */
  class Mean : public MacroData 
  {
  public:
    size_t                  mVectorSize;
    XFormStatAccum *        mpXFormStatAccum;
    size_t                  mNumberOfXFormStatAccums;
    bool                    mUpdatableFromStatAccums;
    FLOAT                   mVector[1];
    
    
    void
    /**
     * @brief Updates the object from the accumulators
     * @param rModelSet ModelSet object which holds the accumulator configuration
     */
    UpdateFromAccums(const ModelSet * pModelSet);
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Variance representation
   */
  class Variance : public MacroData 
  {
  public:
  //  BOOL         diagonal;
    size_t                  mVectorSize;
    XFormStatAccum *        mpXFormStatAccum;
    size_t                  mNumberOfXFormStatAccums;
    bool                    mUpdatableFromStatAccums;
    FLOAT                   mVector[1];
    
    
    void
    /**
     * @brief Updates the object from the accumulators
     * @param rModelSet ModelSet object which holds the accumulator configuration
     */
    UpdateFromAccums(const ModelSet * pModelSet);
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
    size_t                  mNStates;
    FLOAT                   matrix[1];
  
    void
    /**
     * @brief Updates the object from the accumulators
     * @param rModelSet ModelSet object which holds the accumulator configuration
     */
    UpdateFromAccums(const ModelSet * pModelSet);
  };


  /** *************************************************************************
   ** *************************************************************************
   *  @brief XForm statisctics cache
   */
  class XFormStatCache 
  {
  public:
    XFormStatCache *        mpUpperLevelStats;
    XForm *                 mpXForm;
    int                     norm;
    FLOAT *                 mpStats;
  };

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief XForm Instance 
   */
  class XFormInstance : public MacroData 
  {
  public:
    XFormInstance *       mpInput;
    XForm *               mpXForm;
    int                   time;
    XFormInstance  *      next; // Chain of all instances
    XFormStatCache *      mpXFormStatCache;
    size_t                mNumberOfXFormStatCaches;
    size_t                mOutSize;
    int                   statCacheTime;
    char *                memory;
    int                   mTotalDelay;
    
    FLOAT                 out_vec[1]; 
    
    FLOAT *
    XFormPass(FLOAT *in_vec, int time, PropagDir dir);
  
  
    void
    /**
     * @brief Performs desired @c action on the HMM's data which are chosen by @mask
     * @param mask bit mask of flags that tells which data to process by @c action
     * @param nodeNameBuffer 
     * @param action routine to apply on pUserData
     * @param pUserData the data to process by @c action
     */
    Scan( int             mask,
          HMMSetNodeName  nodeNameBuffer,
          ScanAction      action, 
          void *          pUserData);
  
  };

  typedef enum 
  {
    XT_LINEAR,
    XT_COPY,
    XT_BIAS,
    XT_FUNC,
    XT_STACKING,
    XT_COMPOSITE,
  } XFormType;
  

  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief XForm interface (abstract class) to all kinds of transforms
   */
  class XForm : public MacroData 
  {
  public:
    XFormType           mXFormType;
    size_t              mInSize;
    size_t              mOutSize;
    size_t              mMemorySize;
    int                 mDelay;
    
    /**
     * @brief Interface to the evaluation procedure of a concrete XForm
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction) = 0;
  
    void
    /**
     * @brief Performs desired @c action on the HMM's data which are chosen by @mask
     * @param mask bit mask of flags that tells which data to process by @c action
     * @param nodeNameBuffer 
     * @param action routine to apply on pUserData
     * @param pUserData the data to process by @c action
     */
    Scan( int             mask,
          HMMSetNodeName  nodeNameBuffer,
          ScanAction      action, 
          void *          pUserData);
  
  };

  
  class XFormLayer
  {
  public:
    FLOAT *             out_vec;
    size_t              mNBlocks;
    XForm **            block;
  }; // class XFormLayer
  
  
  class CompositeXForm : public XForm 
  {
  public:
    size_t              mNLayers;    
    XFormLayer          layer[1];
    
    /**
     * @brief Composite XForm evaluation 
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction);
  };
  
  class LinearXForm : public XForm 
  {
  public:
    FLOAT               matrix[1];
    
    /**
     * @brief Linear XForm evaluation 
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction);
  };
  
  class BiasXForm : public XForm 
  {
  public:
    FLOAT               mVector[1];
    
    
    /**
     * @brief Bias XForm evaluation 
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction);
  };
  
  
  
  class FuncXForm : public XForm 
  {
  public:
    int                 mFuncId;
  
  
    /**
     * @brief Function XForm evaluation 
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction);
  };
  
  
  
  class CopyXForm : public XForm 
  {
  public:
    int                 indices[1];
  
  
    /**
     * @brief Copy XForm evaluation 
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction);
  
  };
  
  
  class StackingXForm : public XForm 
  {
  public:
    int                 horiz_stack;
  
    /**
     * @brief Stacking XForm evaluation 
     * @param pInputVector pointer to the input vector
     * @param pOutputVector pointer to the output vector
     * @param pMemory pointer to the extra memory needed by the operation
     * @param direction propagation direction (forward/backward)
     */
    virtual FLOAT * 
    Evaluate(FLOAT *    pInputVector, 
             FLOAT *    pOutputVector,
             char *     pMemory,
             PropagDir  direction);
  };
  
  
  class XFormStatsFileNames
  {
  public:
    FILE *      mpStatsP;
    char *      mpStatsN;
    FILE *      mpOccupP;
    char *      mpOccupN;
  };
  
  ///
  /// @{
  class ReadAccumUserData 
  {
  public:
    FILE   *      fp;
    char   *      mpFileName;
    ModelSet *    mpModelSet;
    float         mWeight;
    int           mMmi;
  };
    
  class WriteAccumUserData
  {
  public:
    FILE *        fp;
    char *        mpFileName;
    int           mMmi;
  };

  class ReplaceItemUserData
  {
  public:
    MacroData *        mpOldData;
    MacroData *        mpNewData;
    int           mType;
  };

  struct GlobalStatsUserData
  {
    FLOAT *       observation;
    int           time;
  } ;
  
  class WriteStatsForXFormUserData 
  {
  public:
    LinearXForm *             mpXForm;
    XFormStatsFileNames       mMeanFile;
    XFormStatsFileNames       mCovFile;
    bool                      mBinary;
  };
  /// @}

    
  
                              
  extern bool           gHmmsIgnoreMacroRedefinition;         ///< Controls macro redefinition behavior
  extern const char *   gpHListFilter;                        ///< HMM list Filter command
  extern FunctionTable  gFuncTable[];                         ///< FuncXForm function table
  extern char *         gpKwds[KID_MaxKwdID];                 ///< MMF keyword table
                                                           

                  
  /**
   * @brief Passes the vector through XFormInstance
   */
  FLOAT *     XFormPass(XFormInstance *xformInstance, FLOAT *in_vec, int time, PropagDir dir);
  
  void        PutFlt(FILE *fp, bool binary, FLOAT f);
  void        PutInt(FILE *fp, bool binary, int i);
  void        PutKwd(FILE *fp, bool binary, KeywordID kwdID);
  void        PutNLn(FILE *fp, bool binary);
  
  FLOAT       GetFloat(FILE *fp);
  int         GetInt(FILE *fp);
  char *      GetString(FILE *fp, int eofNotExpected);
  void        RemoveSpaces(FILE *fp);
  void        UngetString(void);
  
  void        ReleaseMacroHash(struct my_hsearch_data *macro_hash);
  Macro *     FindMacro(struct my_hsearch_data *macro_hash, const char *name);
  
  int         CheckKwd(const char *str, KeywordID kwdID);
  void        InitKwdTable();
  KeywordID   ReadDurKind(char *str);
  KeywordID   ReadOutPDFKind(char *str);
  
  /// Action procedures
  /// @{
  void        AllocateXFormStatCachesAndAccums(int, HMMSetNodeName, MacroData * pData, void *);
  void        GlobalStats(int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        NormalizeAccum(int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        NormalizeStatsForXForm(int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        ReadAccum  (int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        ReadStatsForXForm(int macro_type, HMMSetNodeName, void *data, void *userData);
  void        ReleaseItem(int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        ReplaceItem(int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        ResetAccum (int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        WriteAccum (int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  void        WriteStatsForXForm(int macro_type, HMMSetNodeName, MacroData * pData, void *userData);
  /// @}  
  
}; //namespace STK

#endif // STK_Models_h
