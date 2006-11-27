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
  class AlphaBeta;
  class NodeBasicContent;
  

  
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


    


  //############################################################################
  //############################################################################
  // Class list (overview)

  template<typename _NodeContent, typename _LinkContent, NodeRepresentationType _NR, LinkRepresentationType _LR>
    class Link;

  template<typename _NodeContent, typename _LinkContent, NodeRepresentationType _NR, LinkRepresentationType _LR>
    class NodeBasic;

  template<typename _NodeContent, typename _LinkContent, NodeRepresentationType _NR, LinkRepresentationType _LR>
    class Node;


  template< typename                  _NodeContent,
            typename                  _LinkContent,
            NodeRepresentationType    _NodeType, 
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
   *  @brief Network link representation class
   */
  template<class _NodeContent, typename _LinkContent, NodeRepresentationType _NR>
    class Link<_NodeContent, _LinkContent, _NR, LINK_REGULAR> : public _LinkContent
    {
    public:
      typedef FLOAT                     LikeType;
      typedef Node<_NodeContent, _LinkContent, _NR, LINK_REGULAR>   NodeType;

      // Constructors and destructor ...........................................
      Link() 
      : _LinkContent(), mpNode(NULL)
      { }

      Link(NodeType* pNode)
      : _LinkContent(), mpNode(pNode)
      { }

      ~Link()
      { }


      // Accessors ............................................................
      NodeType*
      pNode() const
      { return mpNode; }

      void
      SetNode(NodeType* pNode)
      { mpNode = pNode; }

      /** 
       * @brief Makes the link to point nowhere
       */
      void
      Detach()
      { mpNode = NULL; }

      bool
      PointsNowhere() const
      { return NULL == mpNode; }


    private:
      NodeType *    mpNode;
      //LikeType      mAcousticLike;
      //LikeType      mLmLike;
    }; 
  // Link
  //****************************************************************************

  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network link representation class
   */
  template<typename _NodeContent, typename _LinkContent, NodeRepresentationType _NR>
    class Link<_NodeContent, _LinkContent, _NR, LINK_COMPACT> : public _LinkContent 
    {
    public:
      typedef FLOAT                     LikeType;
      typedef Node<_NodeContent, _LinkContent, _NR, LINK_COMPACT>   NodeType;

      // Constructors and destructor ...........................................
      Link() 
      : _LinkContent()
      { }

      Link(NodeType* pNode)
      : _LinkContent(), mpNode(pNode)
      { }

      ~Link()
      { }

      /** 
       * @brief Initializes the link probabilities
       */
      void
      Init()
      { }


      // Accessors ............................................................
      NodeType*
      pNode() const
      { return mpNode; }

      void
      SetNode(NodeType* pNode)
      { mpNode = pNode; }


    private:
      NodeType *    mpNode;
    }; 
  // Link
  //****************************************************************************


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template<typename _NodeContent, typename _LinkContent, NodeRepresentationType _NR, LinkRepresentationType _LR>
    class NodeBasic : public _NodeContent
    {
    public:
      typedef Link<_NodeContent, _LinkContent, _NR, _LR> LinkType;
      typedef Node<_NodeContent, _LinkContent, _NR, _LR> NodeType;

      NodeBasic() : NodeBasicContent(), mNLinks(0), mpLinks(NULL) {}


      LinkType*
      pLinks() const
      { return mpLinks; }

      int
      NLinks() const
      { return mNLinks; }


      LinkType*
      pBackLinks() const
      { return NULL; }

      int
      NBackLinks() const
      { return 0; }


      /** 
       * @brief Returns number of nodes that really point to me
       * 
       * In this case, the number is not known so -1 is returned.
       */
      int
      NPredecessors() const
      { return -1; }

      /** 
       * @brief Returns number of nodes that I really point to
       * 
       */
      int
      NSuccessors()
      {
        int s(0);

        for (size_t i(0); i < mNLinks; i++)
        {
          if (! mpLinks[i].PointsNowhere())
            ++s;
        }

        return s;
      }

      /** 
       * @brief Finds a link to a given node
       * 
       * @param pNode 
       */
      LinkType*
      pFindLink(const NodeType* pNode)
      {
        LinkType* p_link;

        for (size_t i(0); i < mNLinks; i++)
        {
          p_link = &( mpLinks[i] );

          if (p_link->pNode() == pNode)
            return p_link;
        }

        return NULL;
      }

      LinkType*           mpLinks;
      int                 mNLinks;

    protected:

    }; 
  // class BasicNode
  //****************************************************************************
  


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template<typename _NodeContent, typename _LinkContent, NodeRepresentationType _NR, LinkRepresentationType _LR>
    class Node : public NodeBasic<_NodeContent, _LinkContent, _NR, _LR>
    {
    public:
      typedef       long long  TimingType;
      typedef       Link<_NodeContent, _LinkContent, _NR, _LR> LinkType;

      Node*         mpNext;
      Node*         mpBackNext;
      int           mNBackLinks;
      LinkType*     mpBackLinks;
      int           mAux;


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
  template<typename _NodeContent, typename _LinkContent, LinkRepresentationType _LR>
    class Node<_NodeContent, _LinkContent, NODE_REGULAR, _LR> : public NodeBasic<_NodeContent, _LinkContent, NODE_REGULAR, _LR>
    {
    public:
      typedef       long long  TimingType;
      typedef       Link<_NodeContent, _LinkContent, NODE_REGULAR, _LR> LinkType;
      typedef       Node<_NodeContent, _LinkContent, NODE_REGULAR, _LR> NodeType;

      Node*         mpNext;
      Node*         mpBackNext;
      int           mNBackLinks;
      LinkType*     mpBackLinks;
      int           mAux;

      LinkType*
      pBackLinks() const
      { return mpBackLinks; }

      int
      NBackLinks() const
      { return mNBackLinks; }

      /** 
       * @brief Returns number of nodes that I really point to
       * 
       */
      int
      NPredecessors()
      {
        int s(0);

        for (size_t i(0); i < mNBackLinks; i++)
        {
          if (! mpBackLinks[i].PointsNowhere())
            ++s;
        }

        return s;
      }


      /** 
       * @brief Makes the specified link point nowhere. The target node is
       * updated as well
       * 
       * @param pLink pointer to the link
       *
       * Use Link::PointNowhere() to check whether the link is detached
       */
      void
      DetachLink(LinkType* pLink)
      {
        LinkType* p_back_link;

        if (NULL != (p_back_link = pLink->pNode()->pFindBackLink(this)))
        {
          p_back_link->Detach();
        }

        pLink->Detach();
      }
      

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

      /** 
       * @brief Finds a link to a given node
       * 
       * @param pNode 
       */
      LinkType*
      pFindBackLink(const NodeType* pNode)
      {
        LinkType* p_link;

        for (size_t i(0); i < mNBackLinks; i++)
        {
          p_link = &( mpBackLinks[i] );

          if (p_link->pNode() == pNode)
            return p_link;
        }

        return NULL;
      }

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
  template<typename _NodeContent, typename _LinkContent, LinkRepresentationType _LR>
    class Node<_NodeContent, _LinkContent, NODE_COMPACT, _LR> : public NodeBasic<_NodeContent, _LinkContent, NODE_COMPACT, _LR>
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
        mpPtr = mpPtr->mpBackNext;
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

      iterator
      erase(iterator pos);

    protected:
      _Content*     mpFirst;   ///< self descriptive
      _Content*     mpLast;    ///< self descriptive
        
    };
  // StorageList : public Storage<_Content>
  //****************************************************************************


  template<typename _Content>
    class ArrayStorage
    {
    public:
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
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, 
           LinkRepresentationType _LinkType, NetworkStorageType _NetworkType, 
           template<class> class _StorageType>
    class Network : public _StorageType<Node<_NodeContent, _LinkContent, _NodeType, _LinkType> >
    {
    public:
      typedef Node<_NodeContent, _LinkContent, _NodeType,_LinkType>     NodeType;
      typedef Link<_NodeContent, _LinkContent, _NodeType,_LinkType>     LinkType;
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

      // Destructor ............................................................
      ~Network()
      { 
        if (!mIsExternal)
          Clear(); 
      }

      const bool
      IsCompact() const
      { return mCompactRepresentation; }


      void
      Clear();

      // Various ............................................................... 
      Network&
      TopologicalSort();


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

      void
      PruneNode(NodeType* pNode);

      /** 
       * @brief Safely remove node
       * 
       * @param pNode pointer to a node to be removed
       *
       * The function also delete all links and backlinks to preserve network
       * consistency
       */
      void
      RemoveNode(NodeType* pNode);


      /** 
       * @brief This function unlinks the node, i.e. removes all links
       * 
       * @param pNode 
       */
      void
      IsolateNode(NodeType* pNode);

      
      /** 
       * @brief Insert null node to self links
       * 
       * Insert null node to self links. Self links are not supported
       * in network manipulation functions  
       */
      void
      SelfLinksToNullNodes();


      bool          mCompactRepresentation;

      Network&
      Reverse();

    private:
      /// Indicates, whether network was created from an external structure.
      /// It controls whether destruction of @c this calls the Clear()
      /// function.
      bool          mIsExternal;
    }; // class Network
  
}; // namespace STK

#include "Net.tcc"

#endif // STK_Net_h
