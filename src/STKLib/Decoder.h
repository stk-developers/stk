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
#include "Lattice.h"


#include <list>
#include <vector>

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
//  class Lattice;

  template<typename _NodeType>
    class WordLinkRecord;

  template<typename _NodeType>
    class ActiveNodeRecord;

  template<typename _NetworkType>
    class Decoder;

  template<typename _NodeType>
    class Token;
  
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
  template <typename _NetworkType>
    class Decoder 
    {
    public:

      typedef _NetworkType                            NetworkType;
      typedef typename NetworkType::NodeType          NetworkNode;
      typedef STK::WordLinkRecord<NetworkNode>        WordLinkRecord;
      typedef STK::Token<NetworkNode>                 Token;
      typedef typename STK::ActiveNodeRecord<NetworkNode> 
                                                      ActiveNodeRecord;
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
      void 
      SortNodes();
      
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
      ReestState(NetworkNode  *    pNode,
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
      FWBWRet
      ForwardBackward(FLOAT* pObsMx, int nFrames);
      
      FWBWRet
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
      
      NetworkNode*
      pActivateWordNodesLeadingFrom(NetworkNode* pNode);
      
      void
      ActivateNode(NetworkNode* pNode);

      void
      ActivateModel(NetworkNode* pNode);

      void
      DeactivateModel(NetworkNode* pNode);
      
      void
      DeactivateWordNodesLeadingFrom(NetworkNode* pNode);
      
      void 
      MarkWordNodesLeadingFrom(NetworkNode* node);
      
      bool
      AllWordSuccessorsAreActive();
      
      void
      GenerateToken(NetworkNode* pNode);
    
      void 
      TokenPropagationInit();
      
      void 
      TokenPropagationInNetwork();

      void
      TokenPropagationInModels(FLOAT* observation);
          
      void 
      TokenPropagationDone();

      void 
      AddLinkToLattice(NetworkNode* from, NetworkNode* to, 
          FLOAT lmLike);

      int  
      BackwardPruning(int time, NetworkNode* pNode, int state);
      
      void 
      WriteAlpha(int time, NetworkNode* pNode, int state, Token *token);
      
      void 
      WriteBeta(int time, NetworkNode* pNode, int state, Token *token);



      void
      Wlr2Lattice(WordLinkRecord* pWlr, STK::Lattice& rLattice);

      void 
      Wlr2Lattice_AllocateNodes(WordLinkRecord* pWlr, 
          typename STK::Lattice::NodeType*& rpNode);

      void 
      Wlr2Lattice_EstablishLinks(WordLinkRecord* pWlr);
      
      std::vector<FLOAT> mActiveModelsBestLikes;

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
      Init(ModelSet* pHmms, ModelSet* pHmmsToUpdate/*, bool compactRepresentation = false*/);
      
//      /**
//       * @brief Initializes the network
//       * @param net 
//       * @param first 
//       * @param hmms 
//       * @param hmmsToUptade 
//       */
//
//      void 
//      Init(typename NetworkType::NodeType* pFirstNode, ModelSet* pHmms, 
//          ModelSet* pHmmsToUpdate/*, bool compactRepresentation = false*/);
      
      /**
       * @brief Releases memory occupied by the network resources
       */
      void 
      Clear();


      // accessors ...............................................................
      NetworkType&
      rNetwork()
      { return *mpNetwork; }
      

      NetworkNode*
      pFirst() 
      { return mpNetwork->pFirst(); }

      NetworkNode*
      pLast()
      { return mpNetwork->pLast(); }


      
      const NetworkType&
      Net() const
      { return *mpNetwork; }


      FLOAT   
      (*OutputProbability) (State *state, FLOAT *observation, Decoder *network);

      int     
      (*PassTokenInNetwork)(Token* pFrom, Token* pTo, FLOAT addLogLike, 
          FLOAT acousticLike);

      int     
      (*PassTokenInModel)  (Token* pFrom, Token* pTo, FLOAT addLogLike, 
          FLOAT acousticLike);
      
      // decoding functions ......................................................
      void              
      ViterbiInit();  
      
      void              
      ViterbiStep(FLOAT* pObservation);
      
      FLOAT             
      ViterbiDone(Label** pLabels, Lattice* pNetwork = NULL, bool getTimesFromNetwork = false);
      
      WordLinkRecord*
      TimePruning(int frame_delay);
      
      FLOAT 
      MCEReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, 
          int nFrames, FLOAT weight, FLOAT sigSlope);
              
      FLOAT
      ViterbiReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, 
          int nFrames, FLOAT weight, BasicVector<FLOAT>* pWeightVector);
      
      FLOAT
      BaumWelchReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, 
          int nFrames, FLOAT weight, BasicVector<FLOAT>* pWeightVector);

      FLOAT 
      GetMpeGamma(const Matrix<FLOAT>& rObsMx, Matrix<FLOAT>& rGamma, FLOAT& avgAcc, 
          int nFrames, FLOAT weight, BasicVector<FLOAT>* pWeightVector);



      static int               
      PassTokenMaxForLattices(Token* from, Token* to, FLOAT addLogLike, 
          FLOAT acousticLike);

      static int               
      PassTokenMax(Token* from, Token* to, FLOAT addLogLike, FLOAT acousticLike);

      static int               
      PassTokenSum(Token* from, Token* to, FLOAT addLogLike, FLOAT acousticLike);


      static void              
      KillToken(Token* token);
    
      static FLOAT             
      DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs, Decoder* net);

      static FLOAT             
      DiagCGaussianMixtureDensity(State* state, FLOAT* obs, Decoder* network);

      static FLOAT             
      FromObservationAtStateId(State* state, FLOAT* obs, Decoder* network);


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
      NetworkNode* 
      PhoneNodesToModelNodes(ModelSet * pHmms, ModelSet *pHmmsToUpdate, 
          int& maxStatesInModel);


    public:
      // public atributes ........................................................
      NetworkType*                mpNetwork;
       
      NetworkNode*                mpActiveModels;
      NetworkNode*                mpActiveNodes;
      
      int                         mNActiveTokensForUtterance;
      int                         mNActiveModelsForUtterance;
      int                         mNActiveTokensForObservation;
      int                         mNActiveModelsForObservation;
                              
      int                         mNumberOfNetStates; // Not used for anything important
      Token*                      mpAuxTokens;
      Cache*                      mpOutPCache;
      Cache*                      mpMixPCache;
    
      long                        mTime;
      Token*                      mpBestToken;
      NetworkNode*                mpBestNode;
      
      // Passing parameters
      State*                      mpThreshState;
      FLOAT                       mBeamThresh;
      FLOAT                       mMaxThreshold;
      FLOAT                       mWordThresh;
                                  
      FLOAT                       mWPenalty;  
      FLOAT                       mMPenalty;
      FLOAT                       mPronScale;
      FLOAT                       mLmScale;
      FLOAT                       mPosteriorScale;
      FLOAT                       mTranScale;
      FLOAT                       mOutpScale;
      FLOAT                       mOcpScale;
      FLOAT                       mPruningThresh;
      SearchPathsType             mSearchPaths;
      int                         mMaxActiveModels;
      int                         mMinActiveModels;

      bool                        mLatticeGeneration;
      bool                        mTimePruning;
//      bool                        mCompactRepresentation;

      bool                        mContinuousTokenGeneration;
      bool                        mKeepExitToken;    
      //    MyHSearchData               mLatticeNodeHash;
      //    NodeType *                      mpLatticeLastNode;

                                    
      PropagDirectionType         mPropagDir;
      int                         mAlignment;
      int                         mCollectAlphaBeta;
      
      //  int                         mmi_den_pass;
      AccumType                   mAccumType;
      ModelSet*                   mpModelSet;
      ModelSet*                   mpModelSetToUpdate;

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
  template<typename _Wlr>
  class WlrReference
  {
  public:
    typedef double           LikeType;

    typedef _Wlr             WordLinkRecord;

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
  //typedef std::list<WlrReference>   AltHypList;


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Token representation
   */
  template <typename _NodeType>
    class Token 
    {
    public:
      /// Token's likelihood precision type
      typedef double                            LikeType;
      //typedef Node<NodeContent, LinkContent, LinkArray>  NodeType;
      typedef _NodeType  NodeType;
      typedef STK::WordLinkRecord<NodeType>     WordLinkRecord;
      typedef STK::WlrReference<WordLinkRecord>           WlrReference;
      typedef std::list<WlrReference >    AltHypList;

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
      Token() : mLike(LOG_0), mpTWlr(NULL), mpWlr(NULL), mpAltHyps(NULL) {}
       
      inline const bool
      IsActive() const 
      { return mLike > LOG_MIN; }
      
      /**
       * @brief Returns pointer to an array of this token's labels
       */
      Label*
      pGetLabels(bool getTimesFromNetwork = false, FLOAT total_like = LOG_0);
      
      NodeType*
      pGetLattice();
      
      void
      AddWordLinkRecord(NodeType* pNode, int stateIdx, int time);
      
      void 
      AddAlternativeHypothesis(WordLinkRecord* pWlr);
      
      void 
      AddAlternativeHypothesis(WordLinkRecord* pWlr, LikeType like, 
          LikeType acousticLike);

      void
      Penalize(LikeType penalty);
    };
  // class Token 
  //***************************************************************************
  

  /** *************************************************************************
   ** *************************************************************************
   * @brief 
   */
  template <typename _NodeType>
  class WordLinkRecord 
  {
  public:
    typedef double     LikeType;
    //typedef Node<NodeContent, LinkContent, LinkArray>  NodeType;
    typedef _NodeType  NodeType;
    typedef WordLinkRecord Self;
    typedef std::list<WlrReference<Self> >   AltHypList;

    NodeType*          mpNode;
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

    NodeType*
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
  template<typename _NodeType>
    class ActiveNodeRecord
    {
    public:
      typedef _NodeType                                   NodeType;
      typedef STK::Token<NodeType>                        Token;     


      ActiveNodeRecord(NodeType* pNode) : mpNode(pNode), mIsActiveModel(false), mActiveNodeFlag(0)
      {
        int numOfTokens = mpNode->mC.mType & NT_MODEL ? mpNode->mC.mpHmm->mNStates : 1;
        mpTokens        = new Token[numOfTokens];
        mpExitToken     = &mpTokens[numOfTokens-1];
#ifdef REPORT_TOKEN_ACTIVITY
      mTokensEntered =  mTokensExited = 0;
#endif
      }
      
      ~ActiveNodeRecord() 
      { delete [] mpTokens; }

      NodeType*             mpNode;
      NodeType*             mpNextActiveModel;
      NodeType*             mpPrevActiveModel;
      NodeType*             mpNextActiveNode;
      NodeType*             mpPrevActiveNode;
      
      bool              mIsActiveModel;
      int               mActiveNodeFlag;
      int               mAux;

      Token*            mpExitToken;
      Token*            mpTokens;
      Token*            mpBestToken;
#ifdef REPORT_TOKEN_ACTIVITY
      int               mTokensEntered;
      int               mTokensExited;
      int               mNodeNumber;
      int               mActivationTime;
      void ReportActivity(int time) { 
        printf("Type:%d Name:%s Num:%d Start:%d End:%d Tokens:%d/%d\n", 
               mpNode->mC.mType,
               mpNode->mC.mType == 1 ? ( mpNode->mC.mpPronun == NULL ? "NULL" : mpNode->mC.mpPronun->mpWord->mpName) :
               mpNode->mC.mType == 2 ? mpNode->mC.mpHmm->mpMacro->mpName :
               "<unk>",
               mNodeNumber, mActivationTime, time, mTokensEntered, mTokensExited); }
#endif
    };
  // class ActiveNodeRecord
  //***************************************************************************


  
  
  void
  FindNBestMixtures(State* pState, FLOAT* pObs, NBestRecord* pNBest, 
      size_t nBest, int time);

  template<typename _NetworkType>
    WordLinkRecord<typename _NetworkType::NodeType> *  
    TimePruning(Decoder<_NetworkType>* pNetwork, int frameDelay);
  
  void 
  UpdateXformInstanceStatCaches(XformInstance* xformInstance, FLOAT* pObservation, int time);


/*  FLOAT 
  compute_diag_c_gaussian_density(
    const FLOAT*  pObs, 
    const FLOAT   gConst,
    const FLOAT*  pMean,
    const FLOAT*  pVar,
    const size_t  vSize); */
  

  template<typename _NodeType>
    void 
    FreeWordLinkRecords(WordLinkRecord<_NodeType>* wlr);


  /*
  FLOAT*
  StateOccupationProbability(
      Decoder*      network,
      FLOAT*        observationMx,
      ModelSet*     hmmset,
      int           nFrames,
      FLOAT**       outProbOrMahDist,
      int           getMahalDist);
  */
    
}; // namespace STK

#include "Decoder.tcc"

#endif  // #ifndef STK_Decoder_h
