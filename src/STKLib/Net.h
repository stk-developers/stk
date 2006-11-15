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

  class FWBWR;
  class ActiveNodeRecord;
  
  
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
    NT_TRUE   = 0x40
  }; // NodeType

  
  enum NodeRepresentationType
  {
    NODE_REGULAR              = 0,
    NODE_COMPACT              = 1
  }; // NodeRepresentationType


  enum LinkStorageType
  {
    LINKS_IN_ARRAY            = 0,
    LINKS_IN_LIST             = 1
  }; // enum LinkStorageType


  enum LinkRepresentationType
  {
    LINK_REGULAR              = 0,
    LINK_COMPACT              = 1
  }; // LinkRepresentationType


  /** 
   * @brief Network representation type
   */
  enum NetworkStorageType
  {
    NETWORK_REGULAR           = 0,
    NETWORK_COMPACT           = 1
  }; // enum NetworkType


  enum NotInDictActionType
  {
    WORD_NOT_IN_DIC_UNSET = 0,
    WORD_NOT_IN_DIC_ERROR = 1,
    WORD_NOT_IN_DIC_WARN  = 2,
    PRON_NOT_IN_DIC_ERROR = 4
  }; // NotInDictActionType;
    


  //############################################################################
  //############################################################################
  // Class list (overview)

  template<NodeRepresentationType _NR, LinkRepresentationType _LR>
    class Link;

  template<NodeRepresentationType _NR, LinkRepresentationType _LR>
    class NodeBasic;

  template<NodeRepresentationType _NR, LinkRepresentationType _LR>
    class Node;


  template< NodeRepresentationType    _NodeType, 
            LinkRepresentationType    _LinkType,
            NetworkStorageType        _NetworkType,
            template<class> class     _StorageType > 
    class Network;


  //############################################################################
  //############################################################################
  // CLASS DECLARATIONS
  //############################################################################
  //############################################################################


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
  template<NodeRepresentationType _NR>
    class Link<_NR, LINK_REGULAR> 
    {
    public:
      typedef FLOAT                     LikeType;
      typedef Node<_NR, LINK_REGULAR>   NodeType;

      // Constructors and destructor ...........................................
      Link() 
      : mAcousticLike(0.0), mLmLike(0.0)
      { }

      Link(NodeType* pNode)
      : mpNode(pNode), mAcousticLike(0.0), mLmLike(0.0)
      { }

      ~Link()
      { }

      /** 
       * @brief Initializes the link probabilities
       */
      void
      Init()
      { mLmLike = mAcousticLike = 0.0; }


      // Accessors ............................................................
      NodeType*
      pNode() const
      { return mpNode; }

      void
      SetNode(NodeType* pNode)
      { mpNode = pNode; }

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



    private:
      NodeType *    mpNode;
      LikeType      mAcousticLike;
      LikeType      mLmLike;
    }; 
  // Link
  //****************************************************************************

  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network link representation class
   */
  template<NodeRepresentationType _NR>
    class Link<_NR, LINK_COMPACT> 
    {
    public:
      typedef FLOAT                     LikeType;
      typedef Node<_NR, LINK_COMPACT>   NodeType;

      // Constructors and destructor ...........................................
      Link() 
      : mLmLike(0.0)
      { }

      Link(NodeType* pNode)
      : mpNode(pNode), mLmLike(0.0)
      { }

      ~Link()
      { }

      /** 
       * @brief Initializes the link probabilities
       */
      void
      Init()
      { mLmLike = 0.0; }


      // Accessors ............................................................
      NodeType*
      pNode() const
      { return mpNode; }

      void
      SetNode(NodeType* pNode)
      { mpNode = pNode; }


      /** 
       * @brief Returns link's acoustic likelihood
       */
      const LikeType&
      AcousticLike() const
      { return 0.0; }

      /** 
       * @brief Returns link's LM likelihood
       */
      const LikeType&
      LmLike() const
      { return mLmLike; }

      /** 
       * @brief Sets the link's acoustic likelihood
       *
       * This function really does'nt do anything as compact node does not 
       * contain any acoustic info. The method is implemented for compatibility
       * reasons only!
       */
      void
      SetAcousticLike(const LikeType& like)
      { }

      /** 
       * @brief Sets link's LM likelihood
       */
      void
      SetLmLike(const LikeType& like) 
      { mLmLike = like; }

      /** 
       * @brief Adds like to the link's acoustic likelihood
       *
       * This function really does'nt do anything as compact node does not 
       * contain any acoustic info. The method is implemented for compatibility
       * reasons only!
       */
      void
      AddAcousticLike(const LikeType& like)
      { }

      /** 
       * @brief Adds like to the link's LM likelihood
       */
      void
      AddLmLike(const LikeType& like) 
      { mLmLike += like; }


    private:
      NodeType *    mpNode;
      LikeType      mLmLike;
    }; 
  // Link
  //****************************************************************************

  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template<NodeRepresentationType _NR, LinkRepresentationType _LR>
    class NodeBasic
    {
    public:
      union 
      {
        char*         mpName;
        Hmm*          mpHmm;
        Pronun*       mpPronun;
      };
    
      int                 mType;
      int                 mNLinks;
      Link<NODE_REGULAR, LINK_REGULAR>*   mpLinks;
      
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
  template<NodeRepresentationType _NR, LinkRepresentationType _LR>
    class Node : public NodeBasic<_NR, _LR>
    {
    public:
      typedef       long long  TimingType;

      Node*         mpNext;
      Node*         mpBackNext;
      int           mNBackLinks;
      Link<_NR, _LR>*         mpBackLinks;
      int           mAux;

      //time range when model can be active - apply only for model type
      FLOAT         mPhoneAccuracy;

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

    
#   ifndef EXPANDNET_ONLY
      Hmm*               mpHmmToUpdate;
      FWBWR*             mpAlphaBetaList;
      FWBWR*             mpAlphaBetaListReverse;
#   endif        
    private:
      TimingType     mStart;
      TimingType     mStop;
      
    }; 
  // class Node
  //****************************************************************************


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template<LinkRepresentationType _LR>
    class Node<NODE_REGULAR, _LR> : public NodeBasic<NODE_REGULAR, _LR>
    {
    public:
      typedef       long long  TimingType;

      Node*         mpNext;
      Node*         mpBackNext;
      int           mNBackLinks;
      Link<NODE_REGULAR, _LR>*         mpBackLinks;
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
    private:
      TimingType     mStart;
      TimingType     mStop;
      
    }; 
  // class Node
  //****************************************************************************


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template<LinkRepresentationType _LR>
    class Node<NODE_COMPACT, _LR> : public NodeBasic<NODE_COMPACT, _LR>
    {
      typedef       long long  TimingType;
      
      //time range when model can be active - apply only for model type
      const TimingType&
      Start() const
      { return 0; }

      const TimingType&
      Stop() const
      { return 0; }

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
      typedef const reference                  const_reference;
      typedef value_type *                     pointer;
      typedef const pointer                    const_pointer;
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

      typedef std::bidirectional_iterator_tag   iterator_category;
      typedef _NodeType                         value_type;
      typedef value_type &                      reference;
      typedef const reference                   const_reference;
      typedef value_type *                      pointer;
      typedef const pointer                     const_pointer;
      typedef ptrdiff_t                         difference_type;

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



  /** 
   * @brief Represents a node storage for Network class
   *
   * It uses ListIterator class as its iterator, which assumes that the _Content
   * has mpNext and mpBackNext fields.
   *
   * TODO: Probably make it more general and leave the mpNext and mpBackNext 
   * fields for some parrent class...
   */
  template<typename _Content>
    class ListStorage
    {
    public:
      typedef NodeIterator<_Content, NETWORK_REGULAR>           iterator;
      typedef const iterator                                    const_iterator;
      typedef _Content                                          value_type;
      typedef typename iterator::pointer                        pointer;
      typedef typename iterator::const_pointer                  const_pointer;
      typedef typename iterator::reference                      reference;
      typedef typename iterator::const_reference                const_reference;


      ListStorage()
      : mpFirst(NULL), mpLast(NULL)
      { }

      ListStorage(pointer pInit)
      : mpFirst(pInit), mpLast(NULL)
      { }

      iterator 
      begin()
      { return mpFirst; }

      iterator
      end()
      { return iterator(NULL); }

      // element access
      /**
       *  Returns a read/write reference to the data at the first
       *  element of the %list.
       */
      reference
      front()
      { return *begin(); }

      /**
       *  Returns a read-only (constant) reference to the data at the first
       *  element of the %list.
       */
      const_reference
      front() const
      { return *begin(); }

      /**
       *  Returns a read/write reference to the data at the last element
       *  of the %list.
       */
      reference
      back()
      { 
        return *mpLast;
      }

      /**
       *  Returns a read-only (constant) reference to the data at the last
       *  element of the %list.
       */
      const_reference
      back() const
      { 
        return *mpLast;
      }


      bool
      empty() const
      { return NULL == mpFirst; }


    protected:
      _Content*     mpFirst;   ///< self descriptive
      _Content*     mpLast;    ///< self descriptive
        
    };
  // StorageList : public Storage<_Content>
  //****************************************************************************


  template<typename _Content>
    class ArrayStorage
    {
      typedef NodeIterator<_Content, NETWORK_COMPACT>           iterator;
      typedef const iterator                                    const_iterator;
      typedef _Content                                          value_type;
      typedef typename iterator::pointer                        pointer;
      typedef typename iterator::const_pointer                  const_pointer;
      typedef typename iterator::reference                      reference;
      typedef typename iterator::const_reference                const_reference;

      ArrayStorage()
      : mpFirst(NULL), mNElements(0), mpLast(NULL)
      { }

      ArrayStorage(pointer pInit)
      : mpFirst(pInit), mpLast(NULL)
      { }

      iterator 
      begin()
      { return mpFirst; }

      iterator
      end()
      { return iterator(NULL); }

      // element access
      /**
       *  Returns a read/write reference to the data at the first
       *  element of the %list.
       */
      reference
      front()
      { return *begin(); }

      /**
       *  Returns a read-only (constant) reference to the data at the first
       *  element of the %list.
       */
      const_reference
      front() const
      { return *begin(); }

      /**
       *  Returns a read/write reference to the data at the last element
       *  of the %list.
       */
      reference
      back()
      { 
        return *mpLast;
      }

      /**
       *  Returns a read-only (constant) reference to the data at the last
       *  element of the %list.
       */
      const_reference
      back() const
      { 
        return *mpLast;
      }


      bool
      empty() const
      { return !mNElements; }

    protected:
      _Content*     mpFirst;   ///< self descriptive
      size_t        mNElements;
      _Content*     mpLast;    ///< self descriptive
        
    };
  // StorageList : public Storage<_Content>
  //****************************************************************************



  /** **************************************************************************
   ** **************************************************************************
   * @brief Network encapsulation class
   *
   * The Network class provides basic operations on graph structure. It 
   * encapsulates basic access to the network structure.
   */
  template <NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, 
           NetworkStorageType _NetworkType, template<class> class _StorageType>
    class Network : public _StorageType<Node<_NodeType, _LinkType> >
    {
    public:
      typedef Node<_NodeType,_LinkType>                   NodeType;
      typedef Link<_NodeType,_LinkType>                   LinkType;
      typedef _StorageType<NodeType>                      StorageType;

      //typedef NodeIterator<NodeType, _NetworkType>       iterator;
      typedef typename StorageType::iterator             iterator;
      typedef const NodeIterator<NodeType, _NetworkType>  const_iterator;

      static const NetworkStorageType NetworkType       = _NetworkType;


    public:

      using StorageType::begin;
      using StorageType::end;
      using StorageType::empty;
      using StorageType::front;
      using StorageType::back;


      // Construcotrs ..........................................................
      Network() 
      : StorageType(), mCompactRepresentation(false), mIsExternal(false)
      { }

      /** 
       * @brief Takes control over an existing network specified by @c pNode.
       * 
       * In this case, no deallocation will take place when destructor is called
       */
      Network(NodeType* pNode) 
      : StorageType(pNode), mCompactRepresentation(false),
        mIsExternal(true) 
      { }

      /// Builds a linear network from labels 
      Network(const Label* pLabels, STK::NodeKind  nodeKind)
      { 
        BuildFromLabels(pLabels, nodeKind); 
      }

      // Destructor ............................................................
      ~Network()
      { 
        if (!mIsExternal)
          Clear(); 
      }

      const bool
      IsCompact() const
      { return mCompactRepresentation; }


      // creation and destruction functions ....................................
      void
      BuildFromLabels(const Label* pLabels, NodeKind nodeType);

      void
      Clear();


      // accessors ............................................................. 
      NodeType*
      pFirst()
      { 
        return StorageType::mpFirst; 
      }

      NodeType*
      pLast()
      { 
        return StorageType::mpLast; 
      }

      const NodeType*
      pFirst() const
      { 
        return StorageType::mpFirst; 
      }


      NodeType*
      SetFirst(NodeType* pFirst)
      { 
        // we don't want any memory leaks
        assert(IsEmpty());

        StorageType::mpFirst     = pFirst; 
        mIsExternal = true;
      }

      NodeType*
      SetLast(NodeType* pLast)
      { 
        // we don't want any memory leaks
        //assert(IsEmpty());

        StorageType::mpLast = pLast; 
        //mIsExternal = true;
      }



      // helping functions
      /// Returns true if network is empty
      bool
      IsEmpty() const
      { return this->empty(); }

      
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


      bool          mCompactRepresentation;

    private:

      int 
      LatticeLocalOptimization_ForwardPass(int strictTiming);

      int 
      LatticeLocalOptimization_BackwardPass(int strictTiming);

      Network&
      Reverse();


      /// Indicates, whether network was created from an external structure.
      /// It controls whether destruction of @c this calls the Clear()
      /// function.
      bool          mIsExternal;
    }; // class Network


  // Explicit instantiation of the mostly used network types
  template 
    class Network<NODE_REGULAR, LINK_REGULAR, NETWORK_REGULAR, ListStorage>;

  typedef Network<NODE_REGULAR, LINK_REGULAR, NETWORK_REGULAR, ListStorage>
    RegularNetwork;



  
  // GLOBAL FUNCTIONS
  //
  
  Node<NODE_REGULAR, LINK_REGULAR>* 
  MakeNetworkFromLabels(Label* labels, enum NodeKind nodeKind);
  
  void ExpandWordNetworkByDictionary(
    Node<NODE_REGULAR, LINK_REGULAR>* first,
    MyHSearchData* dict,
    int keep_word_nodes,
    int multiple_pronun);
  
  void ExpandMonophoneNetworkToTriphones(
    Node<NODE_REGULAR, LINK_REGULAR>* first,
    MyHSearchData* nonCDphones,
    MyHSearchData* CDphones);
  
  void LatticeLocalOptimization(
    Node<NODE_REGULAR, LINK_REGULAR>* first,
    int strictTiming,
    int trace_flag);
  
  Node<NODE_REGULAR, LINK_REGULAR>* DiscardUnwantedInfoInNetwork(
    Node<NODE_REGULAR, LINK_REGULAR>* first,
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
      const char*               out_MNF);
  

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
    Node<NODE_REGULAR, LINK_REGULAR>* node,
    STKNetworkOutputFormat format,
    long sampPeriod,
    const char* label_file,
    const char* out_MNF);
  
  void WriteSTKNetworkInOldFormat(
    FILE* flp,
    Node<NODE_REGULAR, LINK_REGULAR>* node,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* label_file,
    const char* out_MNF);
  
  void FreeNetwork(Node<NODE_REGULAR, LINK_REGULAR> *node, bool compactRepresentation = false);
  
  Node<NODE_REGULAR, LINK_REGULAR>*
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
  
  Node<NODE_REGULAR, LINK_REGULAR>*
  ReadSTKNetworkInOldFormat(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name,
    const char* in_MLF);

  Node<NODE_REGULAR, LINK_REGULAR>* 
  find_or_create_node(struct MyHSearchData *node_hash, const char *node_id, Node<NODE_REGULAR, LINK_REGULAR> **last);

  
  Node<NODE_REGULAR, LINK_REGULAR>* 
  ReadHTKLattice(
    FILE* lfp,
    MyHSearchData* word_hash,
    MyHSearchData* phone_hash,
    LabelFormat labelFormat,
    long sampPeriod,
    const char* file_name);
  
  void ComputeAproximatePhoneAccuracy(
    Node<NODE_REGULAR, LINK_REGULAR> *first,
    int type);
  
  void SelfLinksToNullNodes(Node<NODE_REGULAR, LINK_REGULAR> *first);
  int RemoveRedundantNullNodes(Node<NODE_REGULAR, LINK_REGULAR> *first);
  
  void NetworkExpansionsAndOptimizations(
    Node<NODE_REGULAR, LINK_REGULAR> *node,
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
