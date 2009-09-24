#include <cstdlib>
#include "Net.h"

#define SIGNIFICANT_PROB_DIFFERENCE (0.01)

namespace STK
{
  //***************************************************************************
  //***************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    void 
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    Clear() 
    {
      if (IsEmpty()) 
        return;

      NodeContainer::clear();

      return;
        
      /*
      if (!mCompactRepresentation)
      {
        NodeContainer::clear();
      }
      else
      {
        NodeBasic<NodeBasicContent, _LinkContent, LinkArray>*  p_node = reinterpret_cast<NodeBasic<NodeBasicContent, _LinkContent, LinkArray>* >(pFirst()) ;
        NodeBasic<NodeBasicContent, _LinkContent, LinkArray>*  p_tmp_node;

        for(p_tmp_node =  reinterpret_cast<NodeBasic<NodeBasicContent, _LinkContent, LinkArray>* >(p_node); 
            p_tmp_node->rNLinks() != 0; 
            p_tmp_node++) 
        {
          free(p_tmp_node->rpLinks());
        }
        
        free(p_node);
      }
      */
    }



  //***************************************************************************
  //***************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>&
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    Reverse()
    {
      NodeType*  node;
      NodeType*  p_first = pFirst();
      NodeType*  p_last = pLast();
      
      for (node = p_first; node != NULL; node = node->mpBackNext) 
      {
        LinkType*  links   = node->rpLinks();
        int         nlinks  = node->NLinks();

        node->rpLinks()       = node->rpBackLinks();
        node->rNLinks()       = node->rNBackLinks();
        node->rpBackLinks()   = links;
        node->rNBackLinks()   = nlinks;
        NodeType*  next     = node->mpNext;
        node->mpNext        = node->mpBackNext;
        node->mpBackNext    = next;
        node->mAux          = -node->mAux;
        p_last              = node;
      }

      this->splice(p_last, p_first);

      return *this;
    }


  //***************************************************************************
  //***************************************************************************
  template <typename _NodeContent, 
            typename _LinkContent, 
            template<class> class _StorageType, 
            template<class> class _LinkContainer>
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>&
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    TopologicalSort()
    {
      int       i;
      int       j;
      NodeType* p_node;
      NodeType* p_lastnode;

      // Sort nodes in topological order
      // printf("Sorting nodes...\n");
      pFirst()->mAux = 1;
      for (p_lastnode = p_node = pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      {
        for (i = 0; i < p_node->NLinks(); i++) 
        {
          NodeType* p_link_node = p_node->rpLinks()[i].pNode();
          
          if (p_link_node->mAux == 0) 
          {
            for (j = 0; j<p_link_node->rNBackLinks() && 
                        p_link_node->rpBackLinks()[j].pNode()->mAux==1; j++)
            {}
            
            if (j == p_link_node->rNBackLinks()) 
            {
              p_lastnode->mpNext  = p_link_node;
              p_lastnode          = p_link_node;
              p_link_node->mAux   = 1;
              p_link_node->mpNext = NULL;
            }
          }
        }
      }

      return *this;
    }
  // TopSort()


  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    void
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    IsolateNode(iterator iNode)
    {
      if (iNode->rpLinks())
      {
        for (int i(0); i< iNode->NLinks(); i++)
        {
          LinkType* p_link = &( iNode->rpLinks()[i] );
          LinkType* p_back_link;
          
          if (! p_link->PointsNowhere() 
          && (NULL != (p_back_link = p_link->pNode()->pFindBackLink(iNode.mpPtr))))
          {
            p_back_link->Detach();
          }
        }

        // TODO: use delete
        free(iNode->rpLinks());
        iNode->rpLinks() = NULL;
        iNode->rNLinks() = 0;
      }

      if (iNode->rpBackLinks())
      {
        for (int i(0); i< iNode->NBackLinks(); i++)
        {
          LinkType* p_link = &( iNode->rpBackLinks()[i] );
          LinkType* p_back_link;
          
          if (! p_link->PointsNowhere() 
          && (NULL != (p_back_link = p_link->pNode()->pFindLink(iNode.mpPtr))))
          {
            p_back_link->Detach();
          }
        }

        // TODO: use delete
        free(iNode->rpBackLinks());
        iNode->rpBackLinks() = NULL;
        iNode->rNBackLinks() = 0;

      }
    }
  // IsolateNode
  //**************************************************************************


  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    void
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    IsolateNode(NodeType* pNode)
    {
      if (pNode->rpLinks())
      {
        for (int i(0); i< pNode->NLinks(); i++)
        {
          // we delete the fist link since DeleteLink shifts the links to the left
          pNode->DeleteLink(pNode->rpLinks());
          LinkType* p_link = &( pNode->rpLinks()[i] );
          LinkType* p_back_link;
          
          if (! p_link->PointsNowhere() 
          && (NULL != (p_back_link = p_link->pNode()->pFindBackLink(pNode))))
          {
            p_back_link->Detach();
          }
        }

        pNode->Links().clear();
      }

      if (pNode->rpBackLinks())
      {
        for (int i(0); i< pNode->NBackLinks(); i++)
        {
          LinkType* p_link = &( pNode->rpBackLinks()[i] );
          LinkType* p_back_link;
          
          if (! p_link->PointsNowhere() 
          && (NULL != (p_back_link = p_link->pNode()->pFindLink(pNode))))
          {
            p_back_link->Detach();
          }
        }

        pNode->BackLinks().clear();
      }
    }
  // IsolateNode
  //**************************************************************************

#if 0
  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    typename Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>:: iterator
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    RemoveNode(Node* pNode)
    {
      IsolateNode(pNode);
      return erase(iterator(pNode));
    }
  // RemoveNode(Node* pNode);
  //**************************************************************************
#endif

  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    typename Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>:: iterator
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    RemoveNode(iterator iNode)
    {
      IsolateNode(iNode);
      return erase(iNode);
    }
  // RemoveNode(Node* pNode);
  //**************************************************************************

#if 0
  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    void
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    PruneNode(iterator pNode)
    {
      Node* p_aux_node;
      LinkType* p_aux_link;

      // Check whether I am the only successor or predecessor. If yes, 
      // recursively prune forward and backward

      // forward prune
      for (size_t i = 0; i < pNode->NLinks(); i++)
      {
        p_aux_link = &( pNode->rpLinks()[i] );

        // I am the only predecessor of p_aux_node, so I prune it first
        if (! p_aux_link->PointsNowhere())
        {
          p_aux_node = p_aux_link->pNode();

          if (1 == p_aux_node->NPredecessors())
          { 
            // Detach the link first, so no recursive pruning goes back
            p_aux_link->Detach();

            PruneNode(iterator(p_aux_node));
          }
          else
          {
            pNode->DetachLink(p_aux_link);
          }
        }
      }

      // backward prune
      for (size_t i = 0; i < pNode->NBackLinks(); i++)
      {
        p_aux_link = &( pNode->rpBackLinks()[i] );

        // I am the only predecessor of p_aux_node, so I prune it first
        if (! p_aux_link->PointsNowhere())
        {
          p_aux_node = p_aux_link->pNode();

          if (1 == p_aux_node->NSuccessors())
          { 
            // Detach the link first, so no recursive pruning goes back
            p_aux_link->Detach();

            PruneNode(iterator(p_aux_node));
          }
          else
          {
            pNode->DetachLink(p_aux_link);
          }
        }
      }

      NodeContainer::erase(pNode);
    }
  // RemoveNode(Node* pNode);
  //**************************************************************************
#endif

  //**************************************************************************
  //**************************************************************************
  template<typename _Content>
    typename ListStorage<_Content>::iterator
    ListStorage<_Content>::
    erase(iterator pos)
    {
      _Content* p_prev = pos.mpPtr->mpBackNext;
      _Content* p_next = pos.mpPtr->mpNext;
      iterator  tmp    = iterator(pos.mpPtr->mpNext);

      if (NULL != p_prev)
        p_prev->mpNext = p_next;

      if (NULL != p_next)
        p_next->mpBackNext = p_prev;
      
      // TODO: use delete
      free(pos.mpPtr);
      return tmp;
    }
  //**************************************************************************

  //**************************************************************************
  //**************************************************************************
  template<typename _Content>
    void
    ListStorage<_Content>::
    clear()
    {
      _Content* p_tmp;

      while (mpPtr != NULL)
      {
        p_tmp = mpPtr->mpNext;
        mpPtr->~_Content();
        free(mpPtr);
        mpPtr = p_tmp;
      }
      mpLast = NULL;
    }
  //**************************************************************************


  //**************************************************************************
  //**************************************************************************
  template<typename _Content>
    void
    ListStorage<_Content>::
    splice(_Content* pFrom, _Content* pTo)
    {
      // TODO: don't forget to change when fixing List
      mpPtr = pFrom;
      mpLast  = pTo;
    }
  //**************************************************************************
  

  //**************************************************************************
  //**************************************************************************
  template<typename _Content>
    void
    ListStorage<_Content>::
    splice(iterator pos, _Content* pFrom, _Content* pTo)
    {
      _Content* p_aux = (pos == end()) 
                      ? mpLast 
                      : pos.mpPtr->mpBackNext;

      // TODO: don't forget to change when fixing List

      pTo->mpNext = pos.mpPtr;

      if (pos == end())
        pFrom->mpBackNext=mpLast;
      else
        pFrom->mpBackNext = pos.mpPtr->mpBackNext;
        

      if (pos != end())
        pos.mpPtr->mpBackNext = pTo;

      p_aux->mpNext = pTo;
    }
  //**************************************************************************

  //**************************************************************************
  //**************************************************************************
  template<class _LinkContent>
    void
    LinkArray<_LinkContent>::
    defragment()
    {
      size_type last_free = mSize;

      for (size_type i = 0; i < mSize; ++i)
      {
        if (NULL == mpPtr[i].pNode() && last_free > i)
        { // mark the initial inconsistency
          last_free = i;
          continue;
        }

        if (NULL != mpPtr[i].pNode() && i > last_free)
        {
          mpPtr[last_free] = mpPtr[i];
          ++last_free;
        }
      }

      // free the memory if necessary
      //if (mSize != last_free)
      //  mpPtr = static_cast<value_type *> (realloc(mpPtr, last_free * 
      //        sizeof(value_type)));

      // update number of links
      mSize = last_free;
    }


  template< typename                _NodeContent, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    void
    Node<_NodeContent, _LinkContent, _LinkContainer>::
    DeleteLink(LinkType* pLink)
    {
      Node* p_node = pLink->pNode();
      LinkType* p_back_link = p_node->pFindBackLink(this);

      assert (NULL != p_back_link);

      ptrdiff_t offset = p_back_link - p_node->rpBackLinks();

      memmove(p_back_link, p_back_link + 1, (p_node->NBackLinks() - offset - 1) * sizeof(LinkType)); 
      p_node->rNLinks()--;

      offset = pLink - rpLinks();

      memmove(pLink, pLink + 1, (NLinks() - offset) * sizeof(LinkType));
      --rNLinks();
    }
  
  template< typename                _NodeContent, 
            typename                _LinkContent, 
            template<class> class   _LinkContainer>
    void
    Node<_NodeContent, _LinkContent, _LinkContainer>::
    DeleteBackLink(LinkType* pLink)
    {
      Node* p_back_node = pLink->pNode();
      LinkType* p_back_link = p_back_node->pFindLink(this);

      assert (NULL != p_back_link);

      ptrdiff_t offset = p_back_link - p_back_node->rpLinks();

      memmove(p_back_link, p_back_link + 1, (p_back_node->NLinks() - offset - 1) * sizeof(LinkType)); 
      p_back_node->rNLinks()--;

      offset = pLink - rpBackLinks();

      memmove(pLink, pLink + 1, (NBackLinks() - offset) * sizeof(LinkType));
      rNBackLinks() --;
    }
} // namespace STK

