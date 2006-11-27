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

namespace STK
{

  /** **************************************************************************
   ** **************************************************************************
   *  @brief STK network specific output format options
   */
  class STKNetworkOutputFormat 
  {
  public:
    unsigned mNoLMLikes              : 1 ;
    unsigned mNoTimes                : 1 ;
    unsigned mStartTimes             : 1 ;
    unsigned mNoWordNodes            : 1 ;
    unsigned mNoModelNodes           : 1 ;
    unsigned mNoPronunVars           : 1 ;
    unsigned mNoDefaults             : 1 ;
    unsigned mAllFieldNames          : 1 ;
    unsigned mArcDefsToEnd           : 1 ;
    unsigned mArcDefsWithJ           : 1 ;
    unsigned mBase62Labels           : 1 ;
    unsigned mAproxAccuracy          : 1 ;
    unsigned mNoAcousticLikes        : 1 ;
                                     
    //Have no effect yet             
    unsigned mStripTriphones         : 1 ;
    unsigned mLinNodeSeqs            : 1 ;
  };
  // class STKNetworkOutputFormat 
  //****************************************************************************


  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network expansion options definitions
   */
  class ExpansionOptions 
  {
  public:
    unsigned mNoOptimization    : 1;
    unsigned mNoWordExpansion   : 1;
    unsigned mRespectPronunVar  : 1;
    unsigned mRemoveWordsNodes  : 1;
    unsigned mCDPhoneExpansion  : 1;
    unsigned mStrictTiming      : 1;
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
  }; // NodeType


  enum NotInDictActionType
  {
    WORD_NOT_IN_DIC_UNSET = 0,
    WORD_NOT_IN_DIC_ERROR = 1,
    WORD_NOT_IN_DIC_WARN  = 2,
    PRON_NOT_IN_DIC_ERROR = 4
  }; // NotInDictActionType;


  struct NodeBasicContent
  {
    NodeBasicContent() : mpName(NULL), mType(NT_UNDEF), mpAnr(NULL)
    {}

    ~NodeBasicContent()
    { 
      if (NULL != mpAlphaBeta);
        //delete mpAlphaBeta;
    }

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
    union
    {
      ActiveNodeRecord* mpAnr;
      AlphaBeta*        mpAlphaBeta;
    };
#   endif
  };

  class NodeContent : public NodeBasicContent
  {
    typedef       long long  TimingType;
    
    int           mAux;

    void
    Init()
    { mStart = mStop = UNDEF_TIME; }

    //time range when model can be active - apply only for model type
    void
    SetStart(const TimingType& start)
    { mStart = start; }

    void
    SetStop(const TimingType& stop)
    { mStop = stop; }

    const TimingType&
    Start() const
    { return mStart; }

    const TimingType&
    Stop() const
    { return mStop; }


    FLOAT         mPhoneAccuracy;
  
#   ifndef EXPANDNET_ONLY
    Hmm*               mpHmmToUpdate;
    FWBWR*             mpAlphaBetaList;
    FWBWR*             mpAlphaBetaListReverse;
#   endif        

  protected:
    TimingType     mStart;
    TimingType     mStop;
  };


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
    Like() const
    { return mLmLike + mAcousticLike; }

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
  
  template 
    class Network<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR,
          NETWORK_REGULAR, ListStorage>;

  typedef Network<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR,
          NETWORK_REGULAR, ListStorage>
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
    LatticeLocalOptimization(int strictTiming, int trace_flag);


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
      MyHSearchData *         triphHash);

    /** 
     * @brief Remove null nones having less than three predecessors or less 
     * than three successors
     */
    int
    RemoveRedundantNullNodes();

    void 
    SelfLinksToNullNodes();

    private:

      int 
      LatticeLocalOptimization_ForwardPass(int strictTiming);

      int 
      LatticeLocalOptimization_BackwardPass(int strictTiming);

  };
  // class DecoderNetwork
  //***************************************************************************




  Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* 
  MakeNetworkFromLabels(Label* labels, enum NodeKind nodeKind);
  
  void ExpandWordNetworkByDictionary(
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* first,
    MyHSearchData* dict,
    int keep_word_nodes,
    int multiple_pronun);
  
  void ExpandMonophoneNetworkToTriphones(
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* first,
    MyHSearchData* nonCDphones,
    MyHSearchData* CDphones);
  
  void LatticeLocalOptimization(
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* first,
    int strictTiming,
    int trace_flag);
  
  Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* DiscardUnwantedInfoInNetwork(
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* first,
    STKNetworkOutputFormat format);
  

  int 
  getInteger(char *str, char **endPtr, const char *file_name, int line_no);

  float 
  getFloat(char *str, char **endPtr, const char *file_name, int line_no);

  template<class _NetworkType>
    void WriteSTKNetwork(
      FILE*                     flp,
      _NetworkType&             rNetwork,
      STKNetworkOutputFormat    format,
      long                      sampPeriod,
      const char*               label_file,
      const char*               out_MNF,
      const FLOAT&              wPenalty,
      const FLOAT&              lmScale);
  

  template <class _NetworkType>
    void
    ReadSTKNetwork(
      FILE*                     lfp,
      struct MyHSearchData *    word_hash,
      struct MyHSearchData *    phone_hash,
      int                       notInDict,
      LabelFormat               labelFormat,
      long                      sampPeriod,
      const char *              file_name,
      const char *              in_MLF,
      bool                      compactRepresentation,
      _NetworkType&             rNetwork);


    
  void WriteSTKNetwork(
    FILE* flp,
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* node,
    STKNetworkOutputFormat format,
    long sampPeriod,
    const char* label_file,
    const char* out_MNF);
  
  void WriteSTKNetworkInOldFormat(
    FILE* flp,
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* node,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* label_file,
    const char* out_MNF);
  
  void FreeNetwork(Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *node, bool compactRepresentation = false);
  
  Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>*
  ReadSTKNetwork(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    int notInDict,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name,
    const char* in_MLF,
    bool compactRepresentation = false);
  
  Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>*
  ReadSTKNetworkInOldFormat(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name,
    const char* in_MLF);

  Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* 
  find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> **last);

  
  Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* 
  ReadHTKLattice(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name);
  
  void ComputeAproximatePhoneAccuracy(
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *first,
    int type);
  
  void SelfLinksToNullNodes(Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *first);
  int RemoveRedundantNullNodes(Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *first);
  
  void NetworkExpansionsAndOptimizations(
    Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *node,
    ExpansionOptions expOptions,
    STKNetworkOutputFormat out_net_fmt,
    MyHSearchData *dictHash,
    MyHSearchData *nonCDphHash,
    MyHSearchData *triphHash);

  int fprintBase62(FILE *fp, int v);
  
}
// namespace STK
// ****************************************************************************


#include "DecoderNetwork_IO.tcc"
  
#endif // #define STK_Decoder_Network_h

