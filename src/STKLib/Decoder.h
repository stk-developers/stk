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

#ifndef STK_Decoder_h
#define STK_Decoder_h

#define bordel_staff

#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "Net.h"
#include "DecoderNetwork.h"


#include <list>

//#define DEBUG_MSGS
//#define TRACE_TOKENS

#define IS_ACTIVE(token) ((token).mLike > LOG_MIN)



namespace STK
{

  //###########################################################################
  //###########################################################################
  // CLASS OVERVIEW
  //###########################################################################
  //###########################################################################
  class Cache;
  class WordLinkRecord;
  class Decoder;
  class ActiveNodeRecord;
  class Token;
  class Lattice;
  
  //###########################################################################
  //###########################################################################
  // GENERAL CONSTS
  //###########################################################################
  //###########################################################################
  
  
  //###########################################################################
  //###########################################################################
  // GENERAL ENUMS
  //###########################################################################
  //###########################################################################
  typedef enum 
  {
    SP_ALL, 
    SP_TRUE_ONLY
  } SearchPathsType;
  
  
  typedef enum 
  {
    NO_ALIGNMENT    = 0,
    WORD_ALIGNMENT  = 1,
    MODEL_ALIGNMENT = 2,
    STATE_ALIGNMENT = 4,
    FRAME_ALIGNMENT = 8
  } AlignmentType;
  
  
  typedef enum 
  {
    AT_ML=0,
    AT_MPE,
    AT_MFE,
    AT_MCE,
    AT_MMI
  } AccumType;



  
  
  //###########################################################################
  //###########################################################################
  // CLASS DEFINITIONS
  //###########################################################################
  //###########################################################################
  class Cache 
  {
  public:
    FLOAT mValue;
    long  mTime;
  };
  
  
  /** 
   * @brief Description of the NBest mixtures structure
   */
  struct NBestRecord 
  {
    FLOAT    like;
    FLOAT    P;
    size_t   index;
  };

  
  

  /** *************************************************************************
   ** *************************************************************************
   *  @brief Decoder representation
   */
  class Decoder 
  {
  public:
    typedef DecoderNetwork  NetworkType;

  private:
    /**
     * @brief Forward-Backward method return type
     */
    struct FWBWRet 
    {
      double          totLike;            ///< total likelihood
      FloatInLog      avgAccuracy;        ///< average accuracy
    }; // struct FWBWRet

    
    /**
     * @brief Sorts nodes in this network
     */
//    void 
//    SortNodes();
    
    /**
     * @brief Checks for a cyclic structure in network
     */
    int 
    HasCycle();
    

    /** 
     * @brief State reestimation procedure
     * 
     * @param pNode 
     * @param stateIndex 
     * @param logPriorProb 
     * @param updateDir 
     * @param obs 
     * @param obs2 
     */
    void 
    ReestState(NetworkType::NodeType  *    pNode,
               int       stateIndex, 
               FLOAT     logPriorProb, 
               FLOAT     updateDir,
               FLOAT *   obs, 
               FLOAT *   obs2);
      
    /**
     * @brief  Computes likelihood using the forward/backward algorithm
     * @param  pObsMx 
     * @param  nFrames number of frames to use
     * @return @c FWBWRet structure containing total likelihood and average
     *         accuracy
     */
    Decoder::FWBWRet
    ForwardBackward(FLOAT* pObsMx, int nFrames);
    
    Decoder::FWBWRet
    ForwardBackward(const Matrix<FLOAT>& rFeatureMatrix, size_t nFrames);
    
    
    /**
     * @brief Frees 
     */
    void 
    FreeFWBWRecords();
  
    /**
     * @brief Returns true if propagation is forward
     */
    const bool
    InForwardPass() const {return mPropagDir == FORWARD;}
    
    NetworkType::NodeType*
    pActivateWordNodesLeadingFrom(NetworkType::NodeType* pNode);
    
    void
    ActivateNode(NetworkType::NodeType* pNode);

    void
    ActivateModel(NetworkType::NodeType* pNode);

    void
    DeactivateModel(NetworkType::NodeType* pNode);
    
    void
    DeactivateWordNodesLeadingFrom(NetworkType::NodeType* pNode);
    
    void 
    MarkWordNodesLeadingFrom(NetworkType::NodeType* node);
    
    bool
    AllWordSuccessorsAreActive();
    
    void
    GenerateToken(NetworkType::NodeType* pNode);
  
    void 
    TokenPropagationInit();
    
    void 
    TokenPropagationInNetwork();

    void
    TokenPropagationInModels(FLOAT* observation);
        
    void 
    TokenPropagationDone();

    FLOAT 
    DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs);
  
    FLOAT 
    DiagCGaussianMixtureDensity(State* pState, FLOAT* pObs);
    
    FLOAT 
    FromObservationAtStateId(State* pState, FLOAT* pObs);

    void 
    AddLinkToLattice(NetworkType::NodeType* from, NetworkType::NodeType* to, 
        FLOAT lmLike);

    int  
    BackwardPruning(int time, NetworkType::NodeType* pNode, int state);
    
    void 
    WriteAlpha(int time, NetworkType::NodeType* pNode, int state, Token *token);
    
    void 
    WriteBeta(int time, NetworkType::NodeType* pNode, int state, Token *token);

    WordLinkRecord*
    TimePruning(int frame_delay);

  public:

    // the constructor .........................................................
    Decoder()
    { mpNetwork = new NetworkType(); }

    // the destructor ..........................................................
    ~Decoder()
    { delete mpNetwork; }

    
    // init and release operators ..............................................
    /**
     * @brief Initializes the network
     * @param net 
     * @param first 
     * @param hmms 
     * @param hmmsToUptade 
     */
    void 
    Init(ModelSet* pHmms, ModelSet* pHmmsToUpdate, 
        bool compactRepresentation = false);
    
    /**
     * @brief Initializes the network
     * @param net 
     * @param first 
     * @param hmms 
     * @param hmmsToUptade 
     */
    void 
    Init(NetworkType::NodeType* pFirstNode, ModelSet* pHmms, 
        ModelSet* pHmmsToUpdate, bool compactRepresentation = false);
    
    /**
     * @brief Releases memory occupied by the network resources
     */
    void 
    Clear();


    // accessors ...............................................................
    NetworkType&
    rNetwork()
    { return *mpNetwork; }
    

    NetworkType::NodeType*
    pFirst() 
    { return mpNetwork->pFirst(); }

    NetworkType::NodeType*
    pLast()
    { return mpNetwork->pLast(); }


    
    const NetworkType&
    Net() const
    { return *mpNetwork; }


    FLOAT   
    (*OutputProbability) (State *state, FLOAT *observation, Decoder *network);

    int     
    (*PassTokenInNetwork)(Token* pFrom, Token* pTo, FLOAT addLogLike, FLOAT acousticLike);

    int     
    (*PassTokenInModel)  (Token* pFrom, Token* pTo, FLOAT addLogLike, FLOAT acousticLike);
    
    // decoding functions ......................................................
    void              
    ViterbiInit();  
    
    void              
    ViterbiStep(FLOAT* pObservation);
    
    FLOAT             
    ViterbiDone(Label** pLabels, Lattice* pNetwork = NULL);
    
    //FLOAT             
    //ViterbiDone(Label** pLabels, Node ** pLattice = NULL);
    
    FLOAT 
    MCEReest(FLOAT* pObsMx, FLOAT * pObsMx2, int nFrames, FLOAT weight, 
        FLOAT sigSlope);
    
    FLOAT 
    MCEReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, 
        int nFrames, FLOAT weight, FLOAT sigSlope);
            
    FLOAT 
    ViterbiReest(FLOAT* pObsMx, FLOAT* pObsMx2, int nFrames, FLOAT weight);
    
    FLOAT
    ViterbiReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, int nFrames, FLOAT weight);
    
    FLOAT
    BaumWelchReest(FLOAT* pObsMx, FLOAT* pObsMx2, int nFrames, FLOAT weight);

    FLOAT
    BaumWelchReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, int nFrames, FLOAT weight);


    // public atributes ........................................................
    NetworkType*                mpNetwork;
     
    NetworkType::NodeType*     mpActiveModels;
    NetworkType::NodeType*     mpActiveNodes;
    int                         mActiveTokens;
                            
    int                         mNumberOfNetStates; // Not used for anything important
    Token*                      mpAuxTokens;
    Cache*                      mpOutPCache;
    Cache*                      mpMixPCache;
  
    long                        mTime;
    Token*                      mpBestToken;
    NetworkType::NodeType*     mpBestNode;
    
    // Passing parameters
    State*                      mpThreshState;
    FLOAT                       mBeamThresh;
    FLOAT                       mWordThresh;
                                
    FLOAT                       mWPenalty;  
    FLOAT                       mMPenalty;
    FLOAT                       mPronScale;
    FLOAT                       mLmScale;
    FLOAT                       mTranScale;
    FLOAT                       mOutpScale;
    FLOAT                       mOcpScale;
    FLOAT                       mPruningThresh;
    SearchPathsType             mSearchPaths;

    bool                        mLatticeGeneration;
    bool                        mCompactRepresentation;

    bool                        mContinuousTokenGeneration;
    bool                        mKeepExitToken;    
    //    MyHSearchData               mLatticeNodeHash;
    //    Node *                      mpLatticeLastNode;

                                  
    PropagDirectionType         mPropagDir;
    int                         mAlignment;
    int                         mCollectAlphaBeta;
    
    //  int                         mmi_den_pass;
    AccumType                   mAccumType;
    ModelSet*                   mpModelSet;
    ModelSet*                   mpModelSetToUpdate;

  private:
    /** 
     * @brief Maps network phone nodes to hmm models in physical memory
     * 
     * @param pHmms self descriptive
     * @param pHmmsToUpdate self descriptive
     * @param maxStatesInModel self descriptive
     * 
     * @return Pointer to the last node in the mNetwork node set
     */
    NetworkType::NodeType* 
    PhoneNodesToModelNodes(ModelSet * pHmms, ModelSet *pHmmsToUpdate, 
        int& maxStatesInModel);
  }; // class Decoder
  //***************************************************************************
  //***************************************************************************
  

  /** *************************************************************************
   ** *************************************************************************
   *  @brief Word link record reference
   *
   *  This class provides a way that Token connects to a WLR. Besides the 
   *  pointer to the WLR, it also holds various probabilities
   *
   *  In fact, I would be happy if this class dissappeared or was somehow
   *  integreted to WordLinkRecord itself.
   */
  class WlrReference
  {
  public:
    typedef double           LikeType;


    WordLinkRecord*          mpWlr;            ///< Associated word link record
    LikeType                 mLike;            ///< Total likelihood
    LikeType                 mAcousticLike;    ///< Acoustic likelihood


    // CONSTRUCTORS
    // copy constructor
    WlrReference(const WlrReference& rWlrRef) : mpWlr(rWlrRef.mpWlr), 
      mLike(rWlrRef.mLike), mAcousticLike(rWlrRef.mAcousticLike) {}

    // create using separate values 
    WlrReference(WordLinkRecord* pWlr, double  like, double acousticLike) : 
      mpWlr(pWlr), mLike(like), mAcousticLike(acousticLike) {}
    
  }; // class WlrReference


  /** 
   * @brief Data structure used for alternative hpothesis list
   */
  typedef std::list<WlrReference>   AltHypList;


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Token representation
   */
  class Token 
  {
  public:
    /// Token's likelihood precision type
    typedef double                            LikeType;
    typedef Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>  NodeType;

    LikeType                mLike;            ///< Total likelihood
    LikeType                mAcousticLike;    ///< Acoustic likelihood
    FloatInLog              mAccuracy;        ///< Accuracy

#   ifdef bordel_staff
    WordLinkRecord*         mpTWlr;
    FLOAT                   mBestLike;
#   endif

    WordLinkRecord*         mpWlr;            ///< Associated word link record
    AltHypList*             mpAltHyps;        ///< Associated alternative hypotheses
    
  
    /**
     * @brief Returns true if token is active
     */
    Token() : mLike(LOG_0), mpWlr(NULL), mpAltHyps(NULL), mpTWlr(NULL) {}
     
    inline const bool
    IsActive() const 
    { return mLike > LOG_MIN; }
    
    /**
     * @brief Returns pointer to an array of this token's labels
     */
    Label*
    pGetLabels();
    
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>*
    pGetLattice();
    
    void
    AddWordLinkRecord(Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* pNode, int stateIdx, int time);
    
    void 
    AddAlternativeHypothesis(WordLinkRecord* pWlr);
    
    void 
    AddAlternativeHypothesis(WordLinkRecord* pWlr, LikeType like, 
        LikeType acousticLike);
  };
  // class Token 
  //***************************************************************************
  

  /** *************************************************************************
   ** *************************************************************************
   * @brief 
   */
  class WordLinkRecord 
  {
  public:
    typedef double     LikeType;

    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>*              mpNode;
    int                mStateIdx;
    int                mAux;             
    LikeType           mLike;
    LikeType           mAcousticLike;

    long               mTime;
    WordLinkRecord*    mpNext;
    AltHypList*        mpAltHyps;
    int                mNReferences;
  #ifdef DEBUG_MSGS
    WordLinkRecord*    mpTmpNext;
    bool               mIsFreed;
  #endif

    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>*
    pNode() const
    { return mpNode; }
  }; 
  // class WordLinkRecord
  //***************************************************************************
  

  
  /** *************************************************************************
   ** *************************************************************************
   * @brief Active node record
   *
   * Token passing procedure is applied only to active nodes (models).Therefore
   * it is necessary that each node keeps a record of its tokens. 
   * ActiveNodeRecord serves as a token pool.
   */
  class ActiveNodeRecord
  {
  public:
    typedef Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> NodeType;


    NodeType*             mpNode;
    NodeType*             mpNextActiveModel;
    NodeType*             mpPrevActiveModel;
    NodeType*             mpNextActiveNode;
    NodeType*             mpPrevActiveNode;
    
    bool                  mIsActiveModel;
    int                   mActiveNodeFlag;
    int                   mAux;

    Token*                mpExitToken;
    Token*                mpTokens;


    ActiveNodeRecord(NodeType* pNode) : mpNode(pNode), mIsActiveModel(false), mActiveNodeFlag(0)
    {
      int numOfTokens = mpNode->mType & NT_MODEL ? mpNode->mpHmm->mNStates : 1;
      mpTokens        = new Token[numOfTokens];
      mpExitToken     = &mpTokens[numOfTokens-1];
    }
    
    ~ActiveNodeRecord() 
    { delete [] mpTokens; }
  };
  // class ActiveNodeRecord
  //***************************************************************************



  /** *************************************************************************
   ** *************************************************************************
   *  @brief Alpha and Beta pair for forward/backward propagation
   */
  struct AlphaBeta
  {
    typedef  FLOAT   LikeType;

    AlphaBeta() : mAlpha(LOG_0), mBeta(LOG_0)
    {}

    ~AlphaBeta()
    {}

    LikeType             mAlpha;
    LikeType             mBeta;
  };
    
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Alpha and Beta pair for forward/backward propagation augmented
   *  with accuracies for MPE training
   */
  struct AlphaBetaMPE : public AlphaBeta
  {
    FloatInLog        mAlphaAccuracy;
    FloatInLog        mBetaAccuracy;
  };
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Forward-Backward record
   */
  struct FWBWR 
  {
    FWBWR*            mpNext;
    int               mTime;
    AlphaBetaMPE      mpState[1];
  };
  
  
  
  
  void              
  KillToken(Token* token);
  
  void
  FindNBestMixtures(State* pState, FLOAT* pObs, NBestRecord* pNBest, 
      size_t nBest, int time);

  FLOAT             
  DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs, Decoder* net);

  FLOAT             
  DiagCGaussianMixtureDensity(State* state, FLOAT* obs, Decoder* network);

  FLOAT             
  FromObservationAtStateId(State* state, FLOAT* obs, Decoder* network);


#ifdef USE_OLD_TOKEN_PASSING
  int               
  PassTokenMaxForLattices(Token* from, Token* to, FLOAT addLogLike);

  int               
  PassTokenMax(Token* from, Token* to, FLOAT addLogLike);

  int               
  PassTokenSum(Token* from, Token* to, FLOAT addLogLike);
#else
  int               
  PassTokenMaxForLattices(Token* from, Token* to, FLOAT addLogLike, FLOAT acousticLike);

  int               
  PassTokenMax(Token* from, Token* to, FLOAT addLogLike, FLOAT acousticLike);

  int               
  PassTokenSum(Token* from, Token* to, FLOAT addLogLike, FLOAT acousticLike);
#endif


  WordLinkRecord*  
  TimePruning(Decoder* pNetwork, int frameDelay);
  
  void 
  LoadRecognitionNetwork(
      char*         netFileName,
      ModelSet*     hmms,
      Decoder*      net);
  
  void 
  ReadRecognitionNetwork(
      FILE*         fp,
      ModelSet*     hmms,
      Decoder*      network,
      LabelFormat   labelFormat,
      long          sampPeriod,
      const char*   file_name,
      const char*   in_MLF);
  
  FLOAT*
  StateOccupationProbability(
      Decoder*      network,
      FLOAT*        observationMx,
      ModelSet*     hmmset,
      int           nFrames,
      FLOAT**       outProbOrMahDist,
      int           getMahalDist);
    
  /*
  FLOAT MCEReest(
    Decoder *     net,
    FLOAT *       obsMx,
    FLOAT *       obsMx2,
    int           nFrames,
    FLOAT         weight,
    FLOAT         sigSlope);
  
  FLOAT BaumWelchReest(
    Decoder *     net,
    FLOAT *       obsMx,
    FLOAT *       obsMx2,
    int           nFrames,
    FLOAT         weight);
  
  FLOAT ViterbiReest(
    Decoder *     net,
    FLOAT *       observationMx,
    FLOAT *       observationMx2,
    int           nFrames,
    FLOAT         weight);
  */
}; // namespace STK

#endif  // #ifndef STK_Decoder_h
  
