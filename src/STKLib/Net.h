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

  //############################################################################
  //############################################################################
  // Class list (overview)

  template< typename                  _NodeContent, 
            typename _LinkContent, 
            template<class> class     _LinkContainer > 
    class Link;


//  template< typename                  _NodeContent, 
//            typename _LinkContent, 
//            template<class> class     _LinkContainer > 
//    class NodeBasic;


  template< typename                  _NodeContent, 
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
   * @brief Array iterator class
   */
  template<class _LinkType>
    class ArrayIterator
    {
    public:
      typedef ArrayIterator<_LinkType>           _Self;
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
      ArrayIterator() 
      { }

      /// copy constructor
      ArrayIterator(value_type *pPtr) : mpPtr(pPtr)
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

      // Random access iterator requirements
      reference
      operator[](const difference_type& n) const
      { return mpPtr[n]; }

      ArrayIterator&
      operator+=(const difference_type& n)
      { mpPtr += n; return *this; }

      ArrayIterator
      operator+(const difference_type& n) const
      { return ArrayIterator(mpPtr + n); }

      ArrayIterator&
      operator-=(const difference_type& n)
      { mpPtr -= n; return *this; }

      ArrayIterator
      operator-(const difference_type& n) const
      { return ArrayIterator(mpPtr - n); }

        
      pointer mpPtr; // pointer to container value type

    };
  // template<class _NodeType>
  //   class NodeIterator<_NodeType, NETWORK_COMPACT> 
  //****************************************************************************



  //****************************************************************************
  //****************************************************************************
  template<class _Content>
    class Array
    {
    public:
      typedef  _Content                     value_type;
      typedef value_type &                  reference;
      typedef const value_type &            const_reference;
      typedef ArrayIterator<_Content>       iterator;
      typedef const iterator                const_iterator;
      typedef ptrdiff_t                     difference_type;
      typedef int                           size_type;

      // CONSTRUCTORS AND DESTRUCTORS ..........................................
      Array() 
      : mpPtr(NULL), mSize(0)
      { }

      Array(size_type size)
      { 
        if (NULL == (mpPtr = (value_type*) calloc(size, sizeof(value_type))))
          Error("Cannot allocate memory");
      }

      ~Array()
      {
        clear();
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
      { return iterator(mpPtr); }

      /** 
       * @brief Returns an iterator pointing to the begining of the array
       */
      const_iterator
      begin() const
      { return iterator(mpPtr); }

      /** 
       * @brief Returns an iterator pointing past the end of the array
       */
      iterator
      end()
      { return iterator(&mpPtr[mSize]); }

      /** 
       * @brief Returns an iterator pointing past the end of the array
       */
      const_iterator
      end() const
      { return iterator(&mpPtr[mSize]); }

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
        iterator i_tmp = end();
        
        return *(--i_tmp);
      }

      /**
       *  Returns a read-only (constant) reference to the data at the last
       *  element of the %list.
       */
      const_reference
      back() const
      { 
        iterator i_tmp = end();
        
        return *(--i_tmp);
      }

      iterator
      erase(iterator pos)
      {
        assert(false);
      }

      
      /** 
       * @brief Clears the array
       */
      void
      clear()
      {
        if (NULL != mpPtr)
        {
          for (iterator i_tmp = begin(); i_tmp != end(); i_tmp++)
            i_tmp->~value_type();

          free(mpPtr);
        }
        
        mpPtr = NULL;
        mSize = 0;
      }


      // element access
      /**
       *  @brief  Subscript access to the data contained in the %vector.
       *  @param n The index of the element for which data should be
       *  accessed.
       *  @return  Read/write reference to data.
       *
       *  This operator allows for easy, array-style, data access.
       *  Note that data access with this operator is unchecked and
       *  out_of_range lookups are not defined. (For checked lookups
       *  see at().)
       */
      reference
      operator[](size_type n)
      { return *(begin() + n); }

      /**
       *  @brief  Subscript access to the data contained in the %vector.
       *  @param n The index of the element for which data should be
       *  accessed.
       *  @return  Read-only (constant) reference to data.
       *
       *  This operator allows for easy, array-style, data access.
       *  Note that data access with this operator is unchecked and
       *  out_of_range lookups are not defined. (For checked lookups
       *  see at().)
       */
      const_reference
      operator[](size_type n) const
      { return *(begin() + n); }

      /** 
       * @brief Splices the original array given by pointers with this
       * 
       * @param pFrom 
       * @param pTo 
       *
       * Basically this method takes over controll of the original array
       * given by pFrom. This method is temorary only. In future, we want real
       * iterators to be used
       */
      void
      splice(_Content* pFrom, _Content* pTo)
      {
        mpPtr = pFrom;
        mSize = size_type(pTo - pFrom ) + 1;
      }


    public:

      value_type*  mpPtr;
      size_type    mSize;
    };


  /** 
   * @brief 
   */
  template<class _LinkContent>
    class LinkArray : public STK::Array<_LinkContent>
    {
    public:
      typedef  _LinkContent                 value_type;
      typedef value_type &                  reference;
      typedef const value_type &            const_reference;
      typedef ArrayIterator<_LinkContent>   iterator;
      typedef const iterator                const_iterator;
      typedef ptrdiff_t                     difference_type;
      typedef int                           size_type;

      using STK::Array<_LinkContent>::mSize;
      using STK::Array<_LinkContent>::mpPtr;

      // CONSTRUCTORS AND DESTRUCTORS ..........................................

      LinkArray() : STK::Array<_LinkContent>() 
      { }

      LinkArray(size_type size) : STK::Array<_LinkContent>(size) 
      { }

      ~LinkArray()
      { }


      /** 
       * @brief Defragments the storage space. 
       *
       * It is possilbe that some links point nowhere. These link records are
       * replaced by the next non-empty links, so no clusters occur in the 
       * array.
       */
      void
      defragment();
    };
  


  /** **************************************************************************
   ** **************************************************************************
   *  @brief Network link representation class
   */
  template< typename                _NodeType, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    class Link 
    : public _LinkContent
    {
    public:
      typedef FLOAT                     LikeType;
      typedef _NodeType              Node;

      // Constructors and destructor ...........................................
      Link() 
      : _LinkContent(), mpNode(NULL)
      { }

      Link(Node* pNode)
      : _LinkContent(), mpNode(pNode)
      { }

      ~Link()
      { }


      // Accessors ............................................................
      Node*
      pNode() const
      { return mpNode; }

      void
      SetNode(Node* pNode)
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
      Node *    mpNode;

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
    class NodeBasic
    // : public _NodeContent

    {
    public:
      typedef _NodeContent                                          NodeContent;

      typedef NodeBasic                                             Self;
      typedef Link<NodeBasic, _LinkContent, _LinkContainer>         LinkType;
      typedef _LinkContainer<LinkType>                              LinkContainer;
      typedef ptrdiff_t                                             difference_type;
//      typedef NodeBasic<_NodeContent, _LinkContent, _LinkContainer> Node;

      NodeBasic() : NodeContent(), mLinks() 
      { }

      ~NodeBasic()
      { }


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
       * @brief Returns number of nodes that really point to me
       * 
       * In this case, the number is not known so -1 is returned.
       */
      int
      NPredecessors() const
      { return -1; }

      /** 
       * @brief Finds a link to a given node
       * 
       * @param pNode 
       */
      LinkType*
      pFindLink(const NodeBasic* pNode)
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

    public:
      LinkContainer       mLinks;
    }; 
  // class BasicNode
  //****************************************************************************



  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
/*  template<typename _Derivate>
    class NodeBase
    {
    public:
      typedef typename _Derivate::LinkContainer LinkContainer;


    public:
    };
*/


  /** *************************************************************************
   ** *************************************************************************
   *  @brief Network node representation class
   */   
  template< typename                _NodeContent, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    class Node
    //: public NodeBasic<_NodeContent, _LinkContent, _LinkContainer>
    //: public _NodeContent
    {
    public:
      typedef _NodeContent                                    NodeContent;
//      typedef long long                                       TimingType;
      typedef Node                                            Self;
      typedef Link<Node, _LinkContent, _LinkContainer>        LinkType;
      typedef _LinkContainer<LinkType>                        LinkContainer;
      typedef ptrdiff_t                                       difference_type;
      //NodeBase<Node<_NodeContent, _LinkContent, _LinkContainer> > NodeBase;

      Node() : NodeContent()
      { }

      ~Node()
      { }


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
      { return mBackLinks.mpPtr; }

      int
      NBackLinks() 
      { return rNBackLinks(); }
      
      /** 
       * @brief Returns number of nodes that I really point to
       */
      int
      NSuccessors()
      {
        int s(0);

        for (int i(0); i < NLinks(); i++)
        {
          if (! rpLinks()[i].PointsNowhere())
            ++s;
        }

        return s;
      }

      /** 
       * @brief Returns number of nodes that I really point to
       */
      int
      NPredecessors()
      {
        int s(0);

        for (int i(0); i < NBackLinks(); i++)
        {
          if (! rpBackLinks()[i].PointsNowhere())
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
      pFindLink(const Node* pNode)
      {
        LinkType* p_link;

        for (int i(0); i < NLinks(); i++)
        {
          p_link = &( rpLinks()[i] );

          if (p_link->pNode() == pNode)
            return p_link;
        }

        return NULL;
      }


      /** 
       * @brief Finds a link to a given node
       * 
       * @param pNode 
       */
      LinkType*
      pFindBackLink(const Node* pNode)
      {
        LinkType* p_link;

        for (int i(0); i < NBackLinks(); i++)
        {
          p_link = &( rpBackLinks()[i] );

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

      LinkContainer&
      BackLinks()
      {return mBackLinks; }


      typename LinkContainer::size_type&
      rNBackLinks() 
      { return mBackLinks.mSize; }

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


      const NodeContent&
      Content() const
      { return mC; }

      NodeContent&
      Content()
      { return mC; }

    public:
      LinkContainer     mLinks;
      LinkContainer     mBackLinks;
      Node*             mpNext;
      Node*             mpBackNext;
      int               mAux;

      NodeContent       mC;
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
      typedef int                                               size_type;


      ListStorage()
      : mpPtr(NULL), mpLast(NULL)
      { }

      ListStorage(pointer pInit)
      : mpPtr(pInit), mpLast(NULL)
      { }

      iterator 
      begin()
      { return mpPtr; }

      const_iterator 
      begin() const
      { return mpPtr; }

      iterator
      end()
      { return iterator(NULL); }

      const_iterator
      end() const
      { return const_iterator(NULL); }

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
      { return NULL == mpPtr; }

      iterator
      erase(iterator pos);

      void
      clear();

      /// Returns number of elements
      size_type
      size() const
      { int sz = 0; for(iterator i = begin(); i != end(); i++) sz++; return sz; }

      void
      splice(_Content* pFrom, _Content* pTo);

      void
      splice(iterator pos, _Content* pFrom, _Content* pTo);

    protected:
      _Content*     mpPtr;   ///< self descriptive
      _Content*     mpLast;    ///< self descriptive
        
    };
  // StorageList : public Storage<_Content>
  //****************************************************************************


  template<typename _Content>
    class NodeArray : public STK::Array<_Content>
    {
    public:
      /*
      typedef NodeArrayIterator<_Content>           iterator;
      typedef const iterator                                    const_iterator;
      */
      typedef _Content                                          value_type;
      typedef typename STK::Array<_Content>::iterator           iterator;
      typedef typename STK::Array<_Content>::const_iterator     const_iterator;
      typedef typename iterator::pointer                        pointer;
      typedef typename iterator::const_pointer                  const_pointer;
      typedef typename iterator::reference                      reference;
      typedef typename iterator::const_reference                const_reference;


      NodeArray() 
      : STK::Array<_Content>(), mpLast(NULL)
      { }

      NodeArray(pointer pInit)
      : STK::Array<_Content>(),  mpLast(NULL)
      { }


    protected:
      _Content*     mpLast;    ///< self descriptive
        
    };
  // StorageList : public Storage<_Content>
  //****************************************************************************


  template <typename _Network, typename _LinkContainer>
    class NodeLinkContainer 
    : public _LinkContainer
    {
    public:
      typedef _Network                                          Network;


      Network*      mpNetwork;
    };


  /** **************************************************************************
   ** **************************************************************************
   * @brief BidirectionalNodeIterator represents a more sophisticated type of
   *        network iterator, which not only iterates linearily through the 
   *        network nodes, but also offers movements through their links
   *        in both directions
   */
  template <typename _Network>
    class BidirectionalNodeIterator 
    : public _Network::iterator
    {
    public:
      typedef BidirectionalNodeIterator<_Network>               Self;
      typedef _Network                                          Network;
      typedef typename Network::LinkContainer                   LinkContainer;
      typedef typename Network::iterator                        iterator;


      // Constructors .........................................................
      BidirectionalNodeIterator() 
      : _Network::iterator()
      { }

      BidirectionalNodeIterator(iterator iter, Network* pNetwork) 
      : _Network::iterator(iter), mpNetwork(pNetwork)
      { }


      // Accessors ............................................................
      LinkContainer&
      OutLinks()
      { return mpNetwork->LinksFrom(static_cast<iterator>(*this)) ; }

      LinkContainer&
      InLinks()
      { return mpNetwork->LinksTo(static_cast<iterator>(*this)) ; }


    protected:
      // we need to have a refference to the parrent network
      Network*    mpNetwork;

    };



  /** **************************************************************************
   ** **************************************************************************
   * @brief BidirectionalNodeIterator represents a more sophisticated type of
   *        network iterator, which not only iterates linearily through the 
   *        network nodes, but also offers movements through their links
   *        in both directions
   */
  template <typename _Network>
    class BidirectionalLinkIterator 
    : public _Network::link_terator
    {
    public:
      typedef BidirectionalLinkIterator<_Network>               Self;
      typedef BidirectionalNodeIterator<_Network>               NodeIterator;
      typedef _Network                                          Network;
      typedef typename Network::iterator                        node_iterator;
      typedef typename Network::link_iterator                   link_iterator;


      // Constructors .........................................................
      BidirectionalLinkIterator() 
      : _Network::iterator()
      { }

      BidirectionalLinkIterator(link_iterator iter, Network* pNetwork) 
      : _Network::link_iterator(iter), mpNetwork(pNetwork)
      { }


      // Accessors ............................................................
      NodeIterator
      OutNode()
      { return mpNetwork->NodeFrom(static_cast<link_iterator>(*this)) ; }

      NodeIterator
      InNode()
      { return mpNetwork->NodeTo(static_cast<link_iterator>(*this)) ; }


    protected:
      // we need to have a refference to the parrent network
      Network*    mpNetwork;

    };



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
      typedef Network<_NodeContent, _LinkContent, _NodeContainer, _LinkContainer>
                                                                Self;
      typedef Node<_NodeContent, _LinkContent, _LinkContainer>  NodeType;
      typedef Link<NodeType, _LinkContent, _LinkContainer>          LinkType;


      typedef typename NodeType::NodeContent                        NodeContent;
      typedef _LinkContent                                      LinkContent;

      typedef _NodeContainer<NodeType>                              NodeContainer;
      typedef _LinkContainer<LinkType>                          LinkContainer;

      typedef typename NodeContainer::iterator                  iterator;
      typedef const iterator                                    const_iterator;

      typedef typename LinkContainer::iterator                  link_iterator;
      typedef const link_iterator                               const_link_iterator;

      typedef BidirectionalNodeIterator<Self>                   Iterator;


    public:

      // Construcotrs ..........................................................
      Network() 
      : NodeContainer(), /*mCompactRepresentation(false),*/ mIsExternal(false)
      { }

      /** 
       * @brief Takes control over an existing network specified by @c pNode.
       * 
       * In this case, no deallocation will take place when destructor is called
       */
      Network(NodeType* pNode) 
      : NodeContainer(pNode), /*mCompactRepresentation(false),*/
        mIsExternal(true) 
      { }

      // Destructor ............................................................
      ~Network()
      { 
        if (!mIsExternal)
          Clear(); 
      }


      Iterator
      Begin()
      { return Iterator(NodeContainer::begin(), this); }

      Iterator
      End()
      { return Iterator(NodeContainer::end(), this); }


      LinkContainer&
      LinksFrom(iterator iNode)
      { return iNode->Links(); }

      const LinkContainer&
      LinksFrom(iterator iNode) const
      { return iNode->Links(); }


      LinkContainer&
      LinksTo(iterator iNode)
      { return iNode->BackLinks(); }

      const LinkContainer&
      LinksTo(iterator iNode) const
      { return iNode->BackLinks(); }


      Iterator
      NodeFrom(link_iterator iLink)
      { return Iterator(iterator(iLink->pNode()), this); }

      Iterator
      NodeTo(link_iterator iLink)
      { // TODO: Implement this
        assert(false); 
      }

/*      const bool
      IsCompact() const
      { return mCompactRepresentation; }*/


      void
      Clear();

      // Various ............................................................... 
      Network&
      TopologicalSort();


      // accessors ............................................................. 
      NodeType*
      pFirst()
      { 
        return &(NodeContainer::front()); 
      }

      NodeType*
      pLast()
      { 
        //return NodeContainer::mpLast; 
        return &(NodeContainer::back()); 
      }


      const NodeType*
      pFirst() const
      { 
        return &(NodeContainer::front()); 
      }


      void
      SetFirst(NodeType* pFirst)
      { 
        // we don't want any memory leaks
        assert(IsEmpty());

        NodeContainer::mpPtr     = pFirst; 
        mIsExternal = true;
      }

      void
      SetLast(NodeType* pLast)
      { 
        NodeContainer::mpLast     = pLast; 
      }


      // helping functions
      /// Returns true if network is empty
      bool
      IsEmpty() const
      { return this->empty(); }

//      void
//      PruneNode(iterator pNode);

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
//      iterator
//      RemoveNode(NodeType* pNode);


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


//      bool          mCompactRepresentation;

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
