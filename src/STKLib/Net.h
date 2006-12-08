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
  class NodeBasicContent;
  



  //############################################################################
  //############################################################################
  // Class list (overview)

  template< typename _NodeContent, 
            typename _LinkContent, 
            template<class> class     _LinkContainer > 
    class Link;


  template< typename _NodeContent, 
            typename _LinkContent, 
            template<class> class     _LinkContainer > 
    class NodeBasic;


  template< typename _NodeContent, 
            typename _LinkContent, 
            template<class> class     _LinkContainer > 
    class Node;


  template< typename                  _NodeContent,
            typename                  _LinkContent,
            template<class> class     _NodeContainer,
            template<class> class     _LinkContainer > 
    class Network;


  //############################################################################
  //############################################################################
  // CLASS DECLARATIONS
  //############################################################################
  //############################################################################


  /** *************************************************************************
   ** *************************************************************************
   * @brief Node iterator class (REGULAR_COMPACT specialization)
   */
  template<class _LinkType>
    class LinkIterator
    {
    public:
      typedef LinkIterator<_LinkType>           _Self;
      typedef _LinkType                         LinkType;

      typedef std::bidirectional_iterator_tag   iterator_category;
      typedef _LinkType                         value_type;
      typedef value_type &                      reference;
      typedef const reference                   const_reference;
      typedef value_type *                      pointer;
      typedef const pointer                     const_pointer;
      typedef ptrdiff_t                         difference_type;


    public:
      /// default constructor
      LinkIterator() 
      { }

      /// copy constructor
      LinkIterator(value_type *pPtr) : mpPtr(pPtr)
      { }


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

  template<class _LinkContent>
    class LinkArray
    {
    public:
      typedef  _LinkContent                 value_type;
      typedef value_type &                  reference;
      typedef const value_type &            const_reference;
      typedef LinkIterator<_LinkContent>    iterator;
      typedef const iterator                const_iterator;
      typedef ptrdiff_t                     difference_type;
      typedef int                           size_type;

      // CONSTRUCTORS AND DESTRUCTORS ..........................................
      LinkArray() 
      : mpPtr(NULL), mSize(0)
      { }

      LinkArray(size_type size)
      { 
        if (NULL == (mpPtr = (value_type*) calloc(size, sizeof(value_type))))
          Error("Cannot allocate memory");
      }

      ~LinkArray()
      {
        if (NULL != mpPtr)
          free(mpPtr);
      }


      void
      resize(size_type newSize);

      /// Returns number of elements
      size_type
      size() const
      { return mSize; }

      //
      /// Returns maximum number of links. 
      size_type
      max_size() const
      { return size_type(mSize); }


      /// True if empty
      bool
      empty() const
      { return mSize > 0 ? false : true; }



      /** 
       * @brief Returns an iterator pointing to the begining of the array
       */
      iterator
      begin()
      { return iterator(&mpPtr[0]); }

      /** 
       * @brief Returns an iterator pointing past the end of the array
       */
      iterator
      end()
      { return iterator(&mpPtr[mSize]); }

      
      /** 
       * @brief Clears the array
       */
      void
      clear()
      {
        if (NULL != mpPtr)
          free(mpPtr);
        
        mpPtr = NULL;
        mSize = 0;
      }

      /** 
       * @brief Defragments the storage space. 
       *
       * It is possilbe that some links point nowhere. These link records are
       * replaced by the next non-empty links, so no clusters occur in the 
       * array.
       */
      void
      defragment();

    public:

      value_type*  mpPtr;
      size_type    mSize;
    };
  


  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network link representation class
   */
  template<class _NodeContent, typename _LinkContent, 
    template<class> class _LinkContainer>
    class Link 
    : public _LinkContent
    {
    public:
      typedef FLOAT                     LikeType;
      typedef Node<_NodeContent, _LinkContent, _LinkContainer>   NodeType;

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

    }; 
  // Link
  //****************************************************************************



   

  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template< typename                    _NodeContent, 
            typename                    _LinkContent, 
            template<class> class       _LinkContainer>
    class NodeBasic : public _NodeContent
    {
    public:
      typedef Link<_NodeContent, _LinkContent, _LinkContainer>      LinkType;
      typedef Node<_NodeContent, _LinkContent, _LinkContainer>      NodeType;
      typedef _LinkContainer<LinkType>                              LinkContainer;


      NodeBasic() 
      : _NodeContent(), mLinks() {}


      // TODO: This will be somehow replaced
      LinkType*&
      rpLinks() 
      { return mLinks.mpPtr; }

      typename LinkContainer::size_type&
      rNLinks() 
      { return mLinks.mSize; }

      const typename LinkContainer::size_type&
      NLinks() const
      { return mLinks.mSize; }

      LinkType*&
      rpBackLinks() 
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

        for (size_t i(0); i < NLinks(); i++)
        {
          if (! rpLinks()[i].PointsNowhere())
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

        for (size_t i(0); i < NLinks(); i++)
        {
          p_link = &( rpLinks()[i] );

          if (p_link->pNode() == pNode)
            return p_link;
        }

        return NULL;
      }

      
      /** 
       * @brief Gives access to the forward links container
       * 
       * @return 
       */
      LinkContainer&
      Links() 
      { return mLinks; }


      LinkContainer       mLinks;
    }; 
  // class BasicNode
  //****************************************************************************
  




  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template< typename                _NodeContent, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    class Node
    : public NodeBasic<_NodeContent, _LinkContent, _LinkContainer>
    {
    public:
      typedef long long  TimingType;
      typedef Link<_NodeContent, _LinkContent, _LinkContainer> LinkType;
      typedef Node<_NodeContent, _LinkContent, _LinkContainer> NodeType;
      typedef ptrdiff_t                        difference_type;
      typedef _LinkContainer<LinkType>          LinkContainer;


      Node*         mpNext;
      Node*         mpBackNext;

    private:
      LinkContainer     mBackLinks;
     
    public:
      int           mAux;


      LinkType*&
      rpBackLinks() 
      { return mBackLinks.mpPtr; }

      LinkContainer&
      BackLinks()
      {return mBackLinks; }


      typename LinkContainer::size_type&
      rNBackLinks() 
      { return mBackLinks.mSize; }

      int
      NBackLinks() 
      { return rNBackLinks(); }

      /** 
       * @brief Returns number of nodes that I really point to
       * 
       */
      int
      NPredecessors()
      {
        int s(0);

        for (size_t i(0); i < NBackLinks(); i++)
        {
          if (! rpBackLinks()[i].PointsNowhere())
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
      DeleteLink(LinkType* pLink);

      /** 
       * @brief Makes the specified link point nowhere. The target node is
       * updated as well
       * 
       * @param pLink pointer to the link
       *
       * Use Link::PointNowhere() to check whether the link is detached
       */
      void
      DeleteBackLink(LinkType* pLink);


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
      
      /** 
       * @brief Makes the specified link point nowhere. The target node is
       * updated as well
       * 
       * @param pLink pointer to the link
       *
       * Use Link::PointNowhere() to check whether the link is detached
       */
      void
      RemoveLink(LinkType* pLink);

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

        for (size_t i(0); i < NBackLinks(); i++)
        {
          p_link = &( rpBackLinks()[i] );

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





  // class specializations .....................................................

  /** *************************************************************************
   ** *************************************************************************
   * @brief Node iterator class 
   */
  template<class _NodeType>
    class NodeListIterator
    {
    public:
      typedef NodeListIterator<_NodeType>        _Self;
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
      NodeListIterator() 
      : mpPtr() { }

      /// copy constructor
      NodeListIterator(value_type *pPtr) 
      : mpPtr(pPtr) { }


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
    class NodeArrayIterator
    {
    public:
      typedef NodeArrayIterator<_NodeType>        _Self;
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
      NodeArrayIterator() 
      { }

      /// copy constructor
      NodeArrayIterator(value_type *pPtr) : mpPtr(pPtr)
      { }



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
  //   class NodeArrayIterator<_NodeType> 
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
      typedef NodeListIterator<_Content>           iterator;
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

      void
      clear();

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
      typedef NodeArrayIterator<_Content>           iterator;
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
  template< typename                  _NodeContent, 
            typename                  _LinkContent, 
            template<class> class     _NodeContainer, 
            template<class> class     _LinkContainer>
    class Network 
    : public _NodeContainer< Node<_NodeContent, _LinkContent, _LinkContainer> >
    {
    public:
      typedef Node<_NodeContent, _LinkContent, _LinkContainer>     NodeType;
      typedef Link<_NodeContent, _LinkContent, _LinkContainer>     LinkType;
      typedef _NodeContainer<NodeType>                      NodeContainer;

      typedef typename NodeContainer::iterator             iterator;
      typedef const iterator  const_iterator;



    public:

      /*
      using NodeContainer::end;
      using NodeContainer::empty;
      using NodeContainer::front;
      using NodeContainer::back;
      */

      // Construcotrs ..........................................................
      Network() 
      : NodeContainer(), mCompactRepresentation(false), mIsExternal(false)
      { }

      /** 
       * @brief Takes control over an existing network specified by @c pNode.
       * 
       * In this case, no deallocation will take place when destructor is called
       */
      Network(NodeType* pNode) 
      : NodeContainer(pNode), mCompactRepresentation(false),
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
        return NodeContainer::mpFirst; 
      }

      NodeType*
      pLast()
      { 
        return NodeContainer::mpLast; 
      }

      const NodeType*
      pFirst() const
      { 
        return NodeContainer::mpFirst; 
      }


      NodeType*
      SetFirst(NodeType* pFirst)
      { 
        // we don't want any memory leaks
        assert(IsEmpty());

        NodeContainer::mpFirst     = pFirst; 
        mIsExternal = true;
      }

      NodeType*
      SetLast(NodeType* pLast)
      { 
        // we don't want any memory leaks
        //assert(IsEmpty());

        NodeContainer::mpLast = pLast; 
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
      iterator
      RemoveNode(iterator iNode);


      /** 
       * @brief Safely remove node
       * 
       * @param pNode pointer to a node to be removed
       *
       * The function also delete all links and backlinks to preserve network
       * consistency
       */
      iterator
      RemoveNode(NodeType* pNode);


      /** 
       * @brief This function unlinks the node, i.e. removes all links
       * 
       * @param pNode 
       */
      void
      IsolateNode(iterator iNode);


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



  template< typename                _NodeContent, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    void
    Node<_NodeContent, _LinkContent, _LinkContainer>::
    DeleteLink(LinkType* pLink)
    {
      NodeType* p_node = pLink->pNode();
      LinkType* p_back_link = p_node->pFindBackLink(this);

      assert (NULL != p_back_link);

      ptrdiff_t offset = p_back_link - p_node->rpBackLinks();

      memmove(p_back_link, p_back_link + 1, (p_node->NBackLinks() - offset - 1) * sizeof(LinkType)); 
      p_node->rNLinks()--;

      offset = pLink - NodeBasic<_NodeContent, _LinkContent, _LinkContainer>::rpLinks();

      memmove(pLink, pLink + 1, (NodeBasic<_NodeContent, _LinkContent, _LinkContainer>::NLinks() - offset) * sizeof(LinkType));
      NodeBasic<_NodeContent, _LinkContent, _LinkContainer>::rNLinks() --;
    }
  
  template< typename                _NodeContent, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    void
    Node<_NodeContent, _LinkContent, _LinkContainer>::
    DeleteBackLink(LinkType* pLink)
    {
      NodeType* p_back_node = pLink->pNode();
      LinkType* p_back_link = p_back_node->pFindLink(this);

      assert (NULL != p_back_link);

      ptrdiff_t offset = p_back_link - p_back_node->rpLinks();

      memmove(p_back_link, p_back_link + 1, (p_back_node->NLinks() - offset - 1) * sizeof(LinkType)); 
      p_back_node->rNLinks()--;

      offset = pLink - rpBackLinks();

      memmove(pLink, pLink + 1, (NBackLinks() - offset) * sizeof(LinkType));
      rNBackLinks() --;
    }
  
}; // namespace STK

#include "Net.tcc"

#endif // STK_Net_h
