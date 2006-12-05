#include <cstdlib>

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

      if (!mCompactRepresentation)
      {
        NodeContainer::clear();
      }
      else
      {
        NodeBasic<NodeBasicContent, _LinkContent, LinkArray>*  p_node(pFirst());
        NodeBasic<NodeBasicContent, _LinkContent, LinkArray>*  p_tmp_node;

        for(p_tmp_node =  reinterpret_cast<NodeBasic<NodeBasicContent, _LinkContent, LinkArray>* >(p_node); 
            p_tmp_node->rNLinks() != 0; 
            p_tmp_node++) 
        {
          free(p_tmp_node->rpLinks());
        }
        
        free(p_node);
      }

      NodeContainer::mpFirst = NULL;
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
      NodeType*  p_first(NodeContainer::mpFirst);
      NodeType*  p_last(NodeContainer::mpLast);
      
      for (node = NodeContainer::mpFirst; node != NULL; node = node->mpBackNext) 
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

      NodeContainer::mpFirst = p_last;
      NodeContainer::mpLast  = p_first;

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
      size_t    i;
      size_t    j;
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
        for (size_t i(0); i< iNode->NLinks(); i++)
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
        for (size_t i(0); i< iNode->NBackLinks(); i++)
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
        for (size_t i(0); i< pNode->NLinks(); i++)
        {
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
        for (size_t i(0); i< pNode->NBackLinks(); i++)
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


  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    typename Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>:: iterator
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    RemoveNode(NodeType* pNode)
    {
      IsolateNode(pNode);
      erase(iterator(pNode));
    }
  // RemoveNode(NodeType* pNode);
  //**************************************************************************


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
  // RemoveNode(NodeType* pNode);
  //**************************************************************************


  //**************************************************************************
  //**************************************************************************
  template<typename _NodeContent, 
           typename _LinkContent, 
           template<class> class _StorageType, 
           template<class> class _LinkContainer>
    void
    Network<_NodeContent, _LinkContent, _StorageType, _LinkContainer>::
    PruneNode(NodeType* pNode)
    {
      NodeType* p_aux_node;
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

            PruneNode(p_aux_node);
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

            PruneNode(p_aux_node);
          }
          else
          {
            pNode->DetachLink(p_aux_link);
          }
        }
      }

      erase(iterator(pNode));
    }
  // RemoveNode(NodeType* pNode);
  //**************************************************************************


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

      while (mpFirst != NULL)
      {
        p_tmp = mpFirst->mpNext;
        delete mpFirst;
        mpFirst = p_tmp;
      }
      mpLast = NULL;
    }
  //**************************************************************************

} // namespace STK

