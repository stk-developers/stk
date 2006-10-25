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

#ifndef STK_Net_h
#define STK_Net_h

#include "Models.h"
#include "labels.h"
#include "dict.h"
#include "common.h"

namespace STK
{
  // Class list (overview)
  //
  class Node;


  class FWBWR;
  class ActiveNodeRecord;
  
  
  // Enums
  //
  enum NodeType 
  {
    NT_UNDEF  = 0x00,
    NT_WORD   = 0x01,
    NT_MODEL  = 0x02,
    NT_PHONE  = 0x04,
    NT_SUBNET = 0x08,
    NT_TEE    = 0x10,
    NT_STICKY = 0x20,
    NT_TRUE   = 0x40
  }; // NodeType

  
  enum NodeRepresentationType
  {
    NODE_WITH_STATIC                = 0
  }; // NodeRepresentationType


  enum LinkStorageType
  {
    LINKS_IN_ARRAY            = 0,
    LINKS_IN_LIST             = 1
  }; // enum LinkStorageType


  enum LinkRepresentationType
  {
    LINK_BASIC                = 0,
    LINK_WITH_ACOUSTIC_LIKE   = 1
  }; // LinkRepresentationType


  /** 
   * @brief Network representation type
   */
  enum NetworkStorageType
  {
    NETWORK_REGULAR               = 0,
    NETWORK_COMPACT               = 1
  }; // enum NetworkType


  enum NotInDictActionType
  {
    WORD_NOT_IN_DIC_UNSET = 0,
    WORD_NOT_IN_DIC_ERROR = 1,
    WORD_NOT_IN_DIC_WARN  = 2,
    PRON_NOT_IN_DIC_ERROR = 4
  }; // NotInDictActionType;
    

  template<LinkRepresentationType _LR>
    class Link;

  template<
    typename            _NodeType, 
    typename            _LinkType,
    NetworkStorageType  _NetworkType >
    class Network;

  typedef Network<STK::Node, STK::Link<LINK_BASIC>, NETWORK_REGULAR> RegularNetwork;


  // Class declarations
  //


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
  


  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network link representation class
   */
  template<LinkRepresentationType _LR>
    class Link 
    {
    public:
      typedef FLOAT LikeType;


      Node *        mpNode;
      LikeType      mAcousticLike;
      LikeType      mLmLike;


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

    }; 
  // Link
  //****************************************************************************


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  class NodeBasic
  {
  public:
    union {
      char*         mpName;
      Hmm*          mpHmm;
      Pronun*       mpPronun;
    };
  
    int                 mType;
    int                 mNLinks;
    Link<LINK_BASIC>*   mpLinks;
    
#   ifndef NDEBUG
    //id of first emiting state - apply only for model type
    int           mEmittingStateId;
    int           mAux2;
#   endif
#   ifndef EXPANDNET_ONLY    
    ActiveNodeRecord* mpAnr;
#   endif

    NodeBasic() : mpPronun(NULL), mType(NT_UNDEF), mNLinks(0), mpLinks(NULL), 
      mpAnr(NULL) {}
  }; 
  // class BasicNode
  //****************************************************************************
  


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  class Node : public NodeBasic
  {
  public:
    Node*         mpNext;
    Node*         mpBackNext;
    int           mNBackLinks;
    Link<LINK_BASIC>*         mpBackLinks;
    int           mAux;

    //time range when model can be active - apply only for model type
    long long     mStart;
    long long     mStop;
    FLOAT         mPhoneAccuracy;
  
#   ifndef EXPANDNET_ONLY
    Hmm*               mpHmmToUpdate;
    FWBWR*             mpAlphaBetaList;
    FWBWR*             mpAlphaBetaListReverse;
#   endif        
  }; 
  // class Node
  //****************************************************************************


  /** *************************************************************************
   ** *************************************************************************
   * @brief Node iterator class
   */
  template<class _NodeType, NetworkStorageType _StorageType>
    class NodeIterator 
    {
    public:
      typedef NodeIterator<_NodeType, _StorageType>           _Self;
      typedef _NodeType                                       _Node;

      typedef std::bidirectional_iterator_tag  iterator_category;
      typedef _NodeType                        value_type;
      typedef value_type &                     reference;
      typedef value_type *                     pointer;
      typedef ptrdiff_t                        difference_type;

    public:
      /// default constructor
      NodeIterator();

      /// copy constructor
      NodeIterator(value_type *pPtr);


      /// return storage type
      const NetworkStorageType&
      StorageType() const;


      /// return designated object
      reference 
      operator*() const 
      {
        return (*mpPtr);
      }

      /// return a pointer to designated object
      pointer 
      operator->() const 
      {
        return &(*mpPtr);
      }

      /// preincrement
      _Self& 
      operator++();

      /// predecrement
      _Self& 
      operator--();


      /// postincrement
      _Self operator++(int);

      /// postdecrement
      _Self operator--(int);

      /// test for iterator equality
      bool 
      operator==(const _Self& rRight) const;

      /// test for iterator inequality
      bool 
      operator!=(const _Self& rRight) const;

    private:
      pointer mpPtr; // pointer to container value type

    };
  // template<class _NodeType, NetworkStorageType _StorageType>
  //   class NodeIterator 
  //****************************************************************************



  // class specializations .....................................................

  /** *************************************************************************
   ** *************************************************************************
   * @brief Node iterator class (REGULAR_NETWORK specialization)
   */
  template<class _NodeType>
    class NodeIterator<_NodeType, NETWORK_REGULAR> 
    {
    public:
      typedef NodeIterator<_NodeType, NETWORK_REGULAR>        _Self;
      typedef _NodeType                                       _Node;

      typedef std::bidirectional_iterator_tag  iterator_category;
      typedef _NodeType                        value_type;
      typedef value_type &                     reference;
      typedef value_type *                     pointer;
      typedef ptrdiff_t                        difference_type;

    public:
      /// default constructor
      NodeIterator() 
      : mpPtr() { }

      /// copy constructor
      NodeIterator(value_type *pPtr) 
      : mpPtr(pPtr) { }


      /// return storage type
      const NetworkStorageType&
      StorageType() const
      { return NETWORK_REGULAR; }


      /// return designated object
      reference 
      operator*() const 
      {
        return (*mpPtr);
      }

      /// return a pointer to designated object
      pointer 
      operator->() const 
      {
        return &(*mpPtr);
      }

      /// preincrement
      _Self& 
      operator++()
      {
        mpPtr = mpPtr->mpNext;
        return *this;
      };

      /// predecrement
      _Self& 
      operator--()
      {
        mpPtr = mpPtr->mpPrev;
        return *this;
      };


      /// postincrement
      _Self operator++(int) 
      {
        _Self tmp = *this;
        ++*this;
        return (tmp);
      }

      /// postdecrement
      _Self operator--(int) 
      {
        _Self tmp = *this;
        --*this;
        return (tmp);
      }

      /// test for iterator equality
      bool 
      operator==(const _Self& rRight) const 
      {
        return (mpPtr == rRight.mpPtr);
      }

      /// test for iterator inequality
      bool 
      operator!=(const _Self& rRight) const 
      {
        return (!(mpPtr == rRight.mpPtr));
      }


    private:
      pointer mpPtr; // pointer to container value type

    };
  // template<class _NodeType>
  //   class NodeIterator<_NodeType, NETWORK_REGULAR> 
  //****************************************************************************



  /** *************************************************************************
   ** *************************************************************************
   * @brief Node iterator class (REGULAR_COMPACT specialization)
   */
  template<class _NodeType>
    class NodeIterator<_NodeType, NETWORK_COMPACT> 
    {
    public:
      typedef NodeIterator<_NodeType, NETWORK_COMPACT>        _Self;
      typedef _NodeType                                       _Node;

      typedef std::bidirectional_iterator_tag  iterator_category;
      typedef _NodeType                        value_type;
      typedef value_type &                     reference;
      typedef value_type *                     pointer;
      typedef ptrdiff_t                        difference_type;

    public:
      /// default constructor
      NodeIterator() 
      { }

      /// copy constructor
      NodeIterator(value_type *pPtr) : mpPtr(pPtr)
      { }


      /// return storage type
      const NetworkStorageType&
      StorageType() const
      { return NETWORK_COMPACT; }


      /// return designated object
      reference 
      operator*() const 
      {
        return (*mpPtr);
      }

      /// return a pointer to designated object
      pointer 
      operator->() const 
      {
        return &(*mpPtr);
      }

      /// preincrement
      _Self& 
      operator++()
      {
        mpPtr++;
        return *this;
      };

      /// predecrement
      _Self& 
      operator--()
      {
        mpPtr--;
        return *this;
      };


      /// postincrement
      _Self operator++(int) 
      {
        _Self tmp = *this;
        ++*this;
        return (tmp);
      }

      /// postdecrement
      _Self operator--(int) 
      {
        _Self tmp = *this;
        --*this;
        return (tmp);
      }

      /// test for iterator equality
      bool 
      operator==(const _Self& rRight) const 
      {
        return (mpPtr == rRight.mpPtr);
      }

      /// test for iterator inequality
      bool 
      operator!=(const _Self& rRight) const 
      {
        return (!(mpPtr == rRight.mpPtr));
      }


    private:
      pointer mpPtr; // pointer to container value type

    };
  // template<class _NodeType>
  //   class NodeIterator<_NodeType, NETWORK_COMPACT> 
  //****************************************************************************



  /** **************************************************************************
   ** **************************************************************************
   * @brief Network encapsulation class
   *
   * The Network class provides basic operations on graph structure. It 
   * encapsulates basic access to the network structure.
   */
  template <class _NodeType, class _LinkType, NetworkStorageType _NetworkType>
    class Network
    {
    public:
      typedef _NodeType       node_type;
      typedef _LinkType       link_type;

      static const NetworkStorageType StorageType=_NetworkType;

      typedef NodeIterator<node_type, _NetworkType>       iterator;
      typedef const NodeIterator<node_type, _NetworkType> const_iterator;

    public:
      // Construcotrs ..........................................................
      Network() : mpFirst(NULL), mpLast(NULL), mCompactRepresentation(false) 
      { }

      Network(node_type* pNode) : mpFirst(pNode), mpLast(NULL), 
        mCompactRepresentation(false) 
      { }

      /// Builds a linear network from labels 
      Network(const Label* pLabels, STK::NodeType  nodeType)
      { 
        BuildFromLabels(pLabels, nodeType); 
      }

      // Destructor ............................................................
      ~Network()
      { 
        Release(); 
      }


      // Iterators and STL-like interface.......................................
      /// return iterator for beginning of mutable sequence
      iterator begin()
      {	
        //return iterator(mBaseNode.mpNext);
        return iterator(mpFirst);
      }

      /// return iterator for end of mutable sequence
      iterator end()
      {	
        return iterator(NULL);
      }

      bool
      empty() const
      { return mBaseNode.mpNext == &mBaseNode; }


      // creation and destruction functions ....................................
      void
      BuildFromLabels(const Label* pLabels, NodeType nodeType);

      void
      Release();


      // accessors ............................................................. 
      node_type*
      pFirst()
      { 
        return mpFirst; 
      }

      node_type*
      pLast()
      { 
        return mpLast; 
      }

      const node_type*
      pFirst() const
      { 
        return mpFirst; 
      }


      node_type*
      SetFirst(node_type* pFirst)
      { 
        mpFirst = pFirst; 
      }

      node_type*
      SetLast(node_type* pLast)
      { 
        mpLast = pLast; 
      }



      // helping functions
      /// Returns true if network is empty
      bool
      IsEmpty() const
      { return NULL == mpFirst; }

      
      /** 
       * @brief Insert null node to self links
       * 
       * Insert null node to self links. Self links are not supported
       * in network manipulation functions  
       */
      void
      SelfLinksToNullNodes();


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
       * @param dict 
       * @param keep_word_nodes 
       * @param multiple_pronun 
       */
      void
      ExpandByDictionary(MyHSearchData* pDict, bool keepWordNodes, 
          bool multiplePronun);


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



    private:
      node_type*     mpFirst;   ///< self descriptive
      node_type*     mpLast;    ///< self descriptive

      /** 
       * @brief Implicit node
       *
       * Also defines the end marker
       */
      node_type      mBaseNode;    

      bool          mCompactRepresentation;
    }; // class Network



  
  // GLOBAL FUNCTIONS
  //
  
  Node* 
  MakeNetworkFromLabels(Label* labels, enum NodeType node_type);
  
  void ExpandWordNetworkByDictionary(
    Node* first,
    MyHSearchData* dict,
    int keep_word_nodes,
    int multiple_pronun);
  
  void ExpandMonophoneNetworkToTriphones(
    Node* first,
    MyHSearchData* nonCDphones,
    MyHSearchData* CDphones);
  
  void LatticeLocalOptimization(
    Node* first,
    int strictTiming,
    int trace_flag);
  
  Node* DiscardUnwantedInfoInNetwork(
    Node* first,
    STKNetworkOutputFormat format);
  

  static int 
  getInteger(char *str, char **endPtr, const char *file_name, int line_no);

  static float 
  getFloat(char *str, char **endPtr, const char *file_name, int line_no);

  template<class _NetworkType>
    void WriteSTKNetwork(
      FILE*                     flp,
      _NetworkType&             rNetwork,
      STKNetworkOutputFormat    format,
      long                      sampPeriod,
      const char*               label_file,
      const char*               out_MNF);
  

  template <class _NetworkType>
    void
    ReadSTKNetwork(
      FILE*                     lfp,
      _NetworkType&             rNetwork,
      struct MyHSearchData *    word_hash,
      struct MyHSearchData *    phone_hash,
      int                       notInDict,
      LabelFormat               labelFormat,
      long                      sampPeriod,
      const char *              file_name,
      const char *              in_MLF,
      bool                      compactRepresentation);

    
  void WriteSTKNetwork(
    FILE* flp,
    Node* node,
    STKNetworkOutputFormat format,
    long sampPeriod,
    const char* label_file,
    const char* out_MNF);
  
  void WriteSTKNetworkInOldFormat(
    FILE* flp,
    Node* node,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* label_file,
    const char* out_MNF);
  
  void FreeNetwork(Node *node, bool compactRepresentation = false);
  
  Node*
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
  
  Node*
  ReadSTKNetworkInOldFormat(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name,
    const char* in_MLF);

  Node* 
  find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, Node **last);

  
  Node* ReadHTKLattice(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name);
  
  void ComputeAproximatePhoneAccuracy(
    Node *first,
    int type);
  
  void SelfLinksToNullNodes(Node *first);
  int RemoveRedundantNullNodes(Node *first);
  
  void NetworkExpansionsAndOptimizations(
    Node *node,
    ExpansionOptions expOptions,
    STKNetworkOutputFormat out_net_fmt,
    MyHSearchData *dictHash,
    MyHSearchData *nonCDphHash,
    MyHSearchData *triphHash);

  int fprintBase62(FILE *fp, int v);


}; // namespace STK

#include "Net.tcc"
#include "Net_IO.tcc"

#endif // STK_Net_h
