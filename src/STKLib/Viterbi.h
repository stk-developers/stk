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

#ifndef STK_Viterbi_h
#define STK_Viterbi_h

#define bordel_staff

#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "Net.h"

//#define DEBUG_MSGS
//#define TRACE_TOKENS

#define IS_ACTIVE(token) ((token).mLike > LOG_MIN)


namespace STK
{
  class Cache;
  class WLR;
  class SubNet;
  class Network;
  
  
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
    AT_MCE
  } AccumType;
  
  
  class Cache 
  {
  public:
    FLOAT mValue;
    long  mTime;
  };
  
  
  class SubNet 
  {
  public:
    int           mNumberOfNodes;
    Node  *       mpFirst;
    Node  *       mpLast;
    Node  *       mpNodes;
  };
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network representation
   */
  class Network 
  {
  public:
    //Subnet part
    //  int   nnodes;
    Node  *                 mpFirst;
    Node  *                 mpLast;
    //  Node  *nodes;
  
    Node  *                 mpActiveModels;
    Node  *                 mpActiveNodes;
    int                     mActiveTokens;
                            
    int                     mNumberOfNetStates;
    Token *                 mpAuxTokens;
    Cache *                 mpOutPCache;
    Cache *                 mpMixPCache;
  
    long                    mTime;
    Token *                 mpBestToken;
    Node  *                 mpBestNode;
    
    // struct my_hsearch_data   mWordHash;
    
    // Passing parameters
    State *                 mpThreshState;
    FLOAT                   mBeamThresh;
    FLOAT                   mWordThresh;
                            
    FLOAT                   mWPenalty;  
    FLOAT                   mMPenalty;
    FLOAT                   mPronScale;
    FLOAT                   mLmScale;
    FLOAT                   mTranScale;
    FLOAT                   mOutpScale;
    FLOAT                   mOcpScale;
    FLOAT                   mPruningThresh;
    SearchPathsType         mSearchPaths;
                              
    FLOAT   (*OutputProbability) (State *state, FLOAT *observation, Network *network);
    int     (*PassTokenInNetwork)(Token *from, Token *to, FLOAT addLogLike);
    int     (*PassTokenInModel)  (Token *from, Token *to, FLOAT addLogLike);
    
    PropagDirectionType     mPropagDir;
    int                     mAlignment;
    int                     mCollectAlphaBeta;
    
  //  int                     mmi_den_pass;
    AccumType               mAccumType;
    ModelSet *              mpModelSet;
    ModelSet *              mpModelSetToUpdate;
    
    
    //*************************************************************************
    void 
    /**
     * @brief Initializes the network
     * @param net 
     * @param first 
     * @param hmms 
     * @param hmmsToUptade 
     */
    Init(Node * pFirst, ModelSet * pHmms, ModelSet * pHmmsToUptade);
    
    void
    /**
     * @brief Releases memory occupied the network resources
     */
    Release();
  }; // class Network
  //***************************************************************************
  //***************************************************************************
  
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Token representation
   */
  class Token 
  {
  public:
    double                  mLike;              ///< Likelihood
    WLR   *                 wlr;                ///<
    FloatInLog              mAccuracy;          ///< Accuracy
#   ifdef bordel_staff
    WLR   *                 twlr;
    FLOAT                   bestlike;
#   endif
  
    /**
     * @brief Returns true if token is active
     */
    bool
    IsActive() const {return mLike > LOG_MIN;}
    
    /**
     * @brief Returns pointer to an array of this token's labels
     */
    Label *
    pGetLabels();
    
    void 
    /**
     * @brief Kills this token
     */
    Kill();
  };
    
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Alpha and Beta pair for forward/backward propagation
   */
  struct AlphaBeta 
  {
    FLOAT       alpha;
    FLOAT       beta;
    FloatInLog  alphaAccuracy;
    FloatInLog  betaAccuracy;
  };
  
  
  class FWBWR 
  {
  public:
    FWBWR *     mpNext;
    int         mTime;
    AlphaBeta   mpState[1];
  };
  
  
  class WLR 
  {
  public:
    Node  *     mpNode;
    int         mStateIdx;
    FLOAT       mLike;
    long        mTime;
    WLR   *     mpNext;
    int         refs;
  #ifdef DEBUG_MSGS
    WLR   *     mpTmpNext;
    BOOL        freed;
  #endif
  };
  
  
  void LoadRecognitionNetwork(
    char *netFileName,
    ModelSet *hmms,
    Network *net);
  
  void ReadRecognitionNetwork(
    FILE *fp,
    ModelSet *hmms,
    Network *network,
    LabelFormat labelFormat,
    long sampPeriod,
    const char *file_name,
    const char *in_MLF);
  
  void InitNetwork(Network *net, Node *first, ModelSet *hmms, ModelSet *hmmsToUptade);
  void ReleaseNetwork(Network *network);
  void ViterbiInit(Network *network);
  void ViterbiStep(Network *network, FLOAT *observation);
  FLOAT ViterbiDone(Network *network, Label **labels);
  Label *GetLabels(Token *token);
  void KillToken(Token *token);
  
  FLOAT *StateOccupationProbability(
    Network *network,
    FLOAT *observationMx,
    ModelSet *hmmset,
    int nFrames,
    FLOAT **outProbOrMahDist,
    int getMahalDist);
    
  FLOAT MCEReest(
    Network *net,
    FLOAT *obsMx,
    FLOAT *obsMx2,
    int nFrames,
    FLOAT weight,
    FLOAT sigSlope);
  
  FLOAT BaumWelchReest(
    Network *net,
    FLOAT *obsMx,
    FLOAT *obsMx2,
    int nFrames,
    FLOAT weight);
  
  FLOAT ViterbiReest(
    Network *net,
    FLOAT *observationMx,
    FLOAT *observationMx2,
    int nFrames,
    FLOAT weight);
  
  FLOAT DiagCGaussianMixtureDensity(State *state, FLOAT *obs, Network *network);
  FLOAT FromObservationAtStateId(State *state, FLOAT *obs, Network *network);
  int PassTokenMax(Token *from, Token *to, FLOAT addLogLike);
  int PassTokenSum(Token *from, Token *to, FLOAT addLogLike);
  int PassTokenSumUnlogLikes(Token *from, Token *to, FLOAT addLogLike);
  
  WLR *TimePruning(Network *network, int frame_delay);
  
}; // namespace STK

#endif  // #ifndef STK_Viterbi_h
  
