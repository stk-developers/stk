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

#ifndef STK_DecoderNetwork_h
#define STK_DecoderNetwork_h

#include "Net.h"
#include <map>

namespace STK
{
  template<typename _NodeType>
    class ActiveNodeRecord;


  class AlphaBeta;
  class DecoderNetwork;
  class LinkContent;

  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network expansion options definitions
   */
  class ExpansionOptions 
  {
  public:
    unsigned mNoOptimization    : 1;
    unsigned mNoWordExpansion   : 1;
    unsigned mRemoveNulls       : 1;
    unsigned mRespectPronunVar  : 1;
    unsigned mRemoveWordsNodes  : 1;
    unsigned mCDPhoneExpansion  : 1;
    unsigned mStrictTiming      : 1;
    unsigned mNoWeightPushing   : 1;
    unsigned mTraceFlag;
  };
  // class ExpansionOptions 
  //****************************************************************************


  // Enums
  //
  enum NodeKind 
  {
    NT_UNDEF  = 0x00,
    NT_WORD   = 0x01,
    NT_MODEL  = 0x02,
    NT_PHONE  = 0x04,
    NT_SUBNET = 0x08,
    NT_TEE    = 0x10,
    NT_STICKY = 0x20,
    NT_TRUE   = 0x40,
    NT_LATTIMAGIC = 0x80   // this flag is set only for lattices
                           // When set, model nodes act as word nodes
  }; // NodeKind


  enum NotInDictActionType
  {
    WORD_NOT_IN_DIC_UNSET = 0,
    WORD_NOT_IN_DIC_ERROR = 1,
    WORD_NOT_IN_DIC_WARN  = 2,
    PRON_NOT_IN_DIC_ERROR = 4
  }; // NotInDictActionType;

  enum PhoneCorrectnessApproximationType
  {
    MPE_ApproximateAccuracy,
//    MPE_ApproximateError,
    MPE_FrameAccuracy,
    MPE_FrameError
  };
  
  /** *************************************************************************
   ** *************************************************************************
   *  @brief Alpha and Beta pair for forward/backward propagation
   */
  struct AlphaBeta
  {
    typedef  double   LinkType;

    AlphaBeta() : mAlpha(LOG_0), mBeta(LOG_0)
    {}

    ~AlphaBeta()
    {}

    LinkType             mAlpha;
    LinkType             mBeta;
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
  
  

  /** *************************************************************************
   ** *************************************************************************
   * @brief 
   */
  //template <typename _Parent>
    class NodeBasicContent
    {
    public:
      typedef       long long  TimingType;
      typedef ActiveNodeRecord< Node<NodeBasicContent, LinkContent, LinkArray> > ActiveNodeRecordType;

      NodeBasicContent() : mpName(NULL), mType(NT_UNDEF)//, mpAnr(NULL)
      {}

      ~NodeBasicContent()
      { 
        if (NULL != mpAlphaBeta);
          //delete mpAlphaBeta;
      }


      //time range when model can be active - apply only for model type
      void
      SetStart(const TimingType& start)
      {}

      void
      SetStop(const TimingType& stop)
      {}

      TimingType
      Start() const
      { return UNDEF_TIME; }

      TimingType
      Stop() const
      { return UNDEF_TIME; }

      void
      SetPhoneAccuracy(FLOAT pAcc)
      {  }

      FLOAT
      PhoneAccuracy() const
      { return 1.0; }

      Hmm*&
      rpHmmToUpdate()
      { assert(false); }

      FWBWR*&
      rpAlphaBetaList()
      { assert(false); }
    
      FWBWR*&
      rpAlphaBetaListReverse()
      { assert(false); }

    public:
      union 
      {
        char*         mpName;
        Hmm*          mpHmm;
        Pronun*       mpPronun;
      };
    
      int             mType;
      
#   ifndef NDEBUG
      //id of first emiting state - apply only for model type
      int           mEmittingStateId;
      int           mAux2;
#   endif
#   ifndef EXPANDNET_ONLY
      struct //union { We need both in SExpant - node posteriors for one best string
      {
        ActiveNodeRecordType*          mpAnr;
        AlphaBeta*                 mpAlphaBeta; 
      };
#   endif
    };


  /** *************************************************************************
   ** *************************************************************************
   * @brief 
   */
  //template <typename _Parent>
    class NodeContent 
    //: public NodeBasicContent<_Parent>
    {
    public:
      typedef       long long  TimingType;
      typedef ActiveNodeRecord<Node<NodeContent, LinkContent, LinkArray> > ActiveNodeRecordType;


      void
      Init()
      { mStart = mStop = UNDEF_TIME; 
        mpAlphaBeta = NULL; }

      //time range when model can be active - apply only for model type
      void
      SetStart(const TimingType& start)
      { mStart = start; }

      void
      SetStop(const TimingType& stop)
      { mStop = stop; }

      TimingType
      Start() const
      { return mStart; }

      TimingType
      Stop() const
      { return mStop; }

      void
      SetPhoneAccuracy(FLOAT pAcc)
      { mPhoneAccuracy = pAcc; }

      FLOAT
      PhoneAccuracy() const
      { return mPhoneAccuracy; }


      Hmm*&
      rpHmmToUpdate()
      { return mpHmmToUpdate; }

      FWBWR*&
      rpAlphaBetaList()
      { return mpAlphaBetaList; }
    
      FWBWR*&
      rpAlphaBetaListReverse()
      { return mpAlphaBetaListReverse; }

    public:
#   ifndef EXPANDNET_ONLY    
      union 
      {
        char*         mpName;
        Hmm*          mpHmm;
        Pronun*       mpPronun;
      };
    
      int             mType;
      
#   ifndef NDEBUG
      //id of first emiting state - apply only for model type
      int           mEmittingStateId;
      int           mAux2;
#   endif
      struct //union { We need both in SExpant - node posteriors for one best string; AlphaBeta also initialized in Init()
      {
        ActiveNodeRecordType*      mpAnr;
        AlphaBeta*                 mpAlphaBeta;
      };
#   endif

    protected:
#   ifndef EXPANDNET_ONLY
      Hmm*               mpHmmToUpdate;
      FWBWR*             mpAlphaBetaList;
      FWBWR*             mpAlphaBetaListReverse;
#   endif

      FLOAT          mPhoneAccuracy;
      TimingType     mStart;
      TimingType     mStop;
    };


  /** *************************************************************************
   ** *************************************************************************
   * @brief 
   */
  class LinkContent
  {
  public:
    typedef FLOAT  LikeType;
    
    LinkContent()
    : mAcousticLike(0.0), mLmLike(0.0)
    {}

    ~LinkContent()
    {}

    /** 
     * @brief Initializes the link probabilities
     */
    void
    Init()
    { mLmLike = mAcousticLike = 0.0; }

    LikeType
    Like(FLOAT lmScale) const
    { return mLmLike*lmScale + mAcousticLike; }

    /** 
     * @brief Returns link's acoustic likelihood
     */
    const LikeType&
    AcousticLike() const
    { return mAcousticLike; }

    /** 
     * @brief Returns link's LM likelihood
     */
    const LikeType&
    LmLike() const
    { return mLmLike; }

    /** 
     * @brief Sets link's acoustic likelihood
     */
    void
    SetAcousticLike(const LikeType& like)
    { mAcousticLike = like; }
    

    /** 
     * @brief Returns link's LM likelihood
     */
    void
    SetLmLike(const LikeType& like) 
    { mLmLike = like; }

    /** 
     * @brief Sets link's acoustic likelihood
     */
    void
    AddAcousticLike(const LikeType& like)
    { mAcousticLike += like; }
    

    /** 
     * @brief Returns link's LM likelihood
     */
    void
    AddLmLike(const LikeType& like) 
    { mLmLike += like; }



    FLOAT mAcousticLike;
    FLOAT mLmLike;
  };



  //###########################################################################
  //###########################################################################
  // Explicit instantiation of the network types used in decoder
  //###########################################################################
  //###########################################################################
  
  typedef Node<NodeContent, LinkContent, LinkArray>             DecoderNode;
  typedef Link<DecoderNode, LinkContent, LinkArray>             DecoderLink;

  typedef Node<NodeBasicContent, LinkContent, LinkArray>    CompactDecoderNode;
  typedef Link<CompactDecoderNode, LinkContent, LinkArray>      CompactDecoderLink;

  template 
    class Network<NodeBasicContent, LinkContent, NodeArray, LinkArray>;

  typedef Network<NodeBasicContent, LinkContent, NodeArray, LinkArray>
    _CompactDecoderNetwork;


  template 
    class Network<NodeContent, LinkContent, ListStorage, LinkArray>;

  typedef Network<NodeContent, LinkContent, ListStorage, LinkArray>
    _DecoderNetwork;


  /** *************************************************************************
   ** *************************************************************************
   * @brief 
   */
  class DecoderNetwork 
  : public _DecoderNetwork 
  {
  public:


    static int 
    lnkcmp(const void *a, const void *b);

    // Construcotrs ..........................................................
    /// Basic constructor
    DecoderNetwork()
    : _DecoderNetwork()
    {}

    /** 
     * @brief Takes control over an existing network specified by @c pNode.
     * 
     * In this case, no deallocation will take place when destructor is called
     */
    DecoderNetwork(NodeType* pNode) 
    : _DecoderNetwork(pNode)
    { }

    /// Builds a linear network from labels 
    DecoderNetwork(const Label* pLabels, STK::NodeKind  nodeKind)
    { 
      BuildFromLabels(pLabels, nodeKind); 
    }

    // Destructor ............................................................
    ~DecoderNetwork()
    {}
    

    // creation and destruction functions ....................................
    void
    BuildFromLabels(const Label* pLabels, NodeKind nodeType);

    /** 
     * @brief Self explanative
     * 
     * @param dict 
     * @param keep_word_nodes 
     * @param multiple_pronun 
     */
    void
    ExpandByDictionary(MyHSearchData* pDict, bool keepWordNodes, 
        bool multiplePronun);

    /** 
     * @brief Performs various optimizations on lattice
     * 
     * @param strictTiming 
     * @param trace_flag 
     */
    void
    LatticeLocalOptimization(const ExpansionOptions &expOptions,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale);


    /** 
     * @brief Self explanative
     * 
     * @param nonCDphones 
     * @param CDphones 
     */
    void
    ExpandMonophonesToTriphones(MyHSearchData *nonCDphones, 
        MyHSearchData *CDphones);


    /** 
     * @brief Discards unwanted information in network records
     * 
     * @param rFormat 
     *
     * The function discard the information in network records that is not to be
     * saved to the output. This should allow for more effective network
     * optimization, which will be run after calling this function and before
     * saving network to file.
     */
    void
    DiscardUnwantedInfo(const STKNetworkOutputFormat& rFormat);


    /** 
     * @brief 
     * 
     * @param expOptions 
     * @param out_net_fmt 
     * @param wordHash 
     * @param nonCDphHash 
     * @param triphHash 
     */
    void 
    ExpansionsAndOptimizations(
      ExpansionOptions        expOptions,
      const STKNetworkOutputFormat&  rFormat,
      MyHSearchData *         wordHash,
      MyHSearchData *         nonCDphHash,
      MyHSearchData *         triphHash,
      FLOAT                   wordPenalty = 0.0,
      FLOAT                   modelPenalty = 0.0,
      FLOAT                   lmScale = 1.0,
      FLOAT                   posteriorScale = 1.0);

    /** 
     * @brief Remove null nones having less than three predecessors or less 
     * than three successors
     */
    int
    RemoveRedundantNullNodes(bool removeAllNullNodes,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale);

    void 
    SelfLinksToNullNodes();

    void 
    ComputePhoneCorrectnes(PhoneCorrectnessApproximationType approxType, MyHSearchData *silencePhones);
//    ComputeApproximatePhoneAccuracy(int type);
    
    /** 
     * @brief Performs forward-backward to get posterior probabilities
     * 
     * @return 
     */
    FLOAT
    ForwardBackward(
      FLOAT wordPenalty = 0.0,
      FLOAT modelPenalty = 0.0 ,
      FLOAT lmScale = 1.0,
      FLOAT posteriorScale = 1.0,
      bool viterbi = true);

    
    /** 
     * @brief Free records with posterior probabilities
     * 
     */
    void
    FreePosteriors();

    /** 
     * @brief Prunes the lattice based on the posterior probabilities
     * 
     * @param thresh 
     */
    void
    PosteriorPrune(
      const FLOAT& thresh,
      FLOAT wordPenalty = 0.0,
      FLOAT modelPenalty = 0.0,
      FLOAT lmScale = 1.0,
      FLOAT posteriorScale = 1.0);
      
      
    void
    PosteriorExpectedCounts(
      std::map<char *,float> &countMap);

    void
    DensityPrune(
      const int maxDensity,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale);
    
    void 
    TopologicalSort();

    void 
    SortNodes();

  private:
    int 
    LatticeLocalOptimization_ForwardPass(const ExpansionOptions &expOptions,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale);

    int 
    LatticeLocalOptimization_BackwardPass(const ExpansionOptions &expOptions,
      FLOAT wordPenalty,
      FLOAT modelPenalty,
      FLOAT lmScale,
      FLOAT posteriorScale);

  };
  // class DecoderNetwork
  //***************************************************************************





  int 
  getInteger(char *str, char **endPtr, const char *file_name, int line_no);

  float 
  getFloat(char *str, char **endPtr, const char *file_name, int line_no);

  int 
  RemoveCommentLines(FILE *fp);


  template<class _NetworkType>
    void WriteSTKNetwork(
      FILE*                     flp,
      _NetworkType&             rNetwork,
      STKNetworkOutputFormat    format,
      long                      sampPeriod,
      const char*               label_file,
      const char*               out_MNF,
      const FLOAT&              wordPenalty,
      const FLOAT&              modelPenalty,
      const FLOAT&              lmScale,
      const FLOAT               posteriorScale = 1.0);
  

  template<class _NetworkType>
    void WriteSTKNetworkInOldFormat(
      FILE* flp,
      _NetworkType&             rNetwork,
//      LabelFormat labelFormat,
      STKNetworkOutputFormat    format,
      long sampPeriod,
      const char* label_file,
      const char* out_MNF);


  template <class _NetworkType>
    void
    ReadSTKNetwork(
      FILE*                     lfp,
      struct MyHSearchData *    word_hash,
      struct MyHSearchData *    phone_hash,
      int                       notInDict,
//      LabelFormat               labelFormat,
      STKNetworkOutputFormat    format,
      long                      sampPeriod,
      const char *              file_name,
      const char *              in_MLF,
      bool                      compactRepresentation,
      _NetworkType&             rNetwork);


  template <class _NetworkType>
    void
    ReadSTKNetworkInOldFormat(
      FILE*                     lfp,
      MyHSearchData*            word_hash,
      MyHSearchData*            phone_hash,
//      LabelFormat               labelFormat,
      STKNetworkOutputFormat    format,
      long                      sampPeriod,
      const char*               file_name,
      const char*               in_MLF,
      _NetworkType&             rNetwork);

    
  
//  void FreeNetwork(DecoderNetwork::Node *node, bool compactRepresentation = false);
  

  //DecoderNetwork::Node* 
  //find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, DecoderNetwork::Node **last);

  

  int fprintBase62(FILE *fp, int v);
  
}
// namespace STK
// ****************************************************************************


#include "DecoderNetwork_IO.tcc"
  
#endif // #define STK_Decoder_Network_h

