#include <cstdlib>

#define SIGNIFICANT_PROB_DIFFERENCE (0.01)

namespace STK
{
  /*
  template<LinkRepresentationType _LR>
    Link<NODE_REGULAR, _LR>*
    Node<NODE_REGULAR, _LR>::
    pFindLink(const Node* pNode)
    {
      LinkType* p_link;

      for (size_t i(0); i < NLinks(); i++)
      {
        p_link = &( mpLinks[i] );

        if (p_link->pNode() == pNode)
          return p_link;
      }

      return NULL;
    }
    */

  
  //***************************************************************************
  //***************************************************************************
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, 
           NetworkStorageType _NetworkType, template<class> class _StorageType>
    void 
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>::
    Clear() 
    {
      if (IsEmpty()) 
        return;

      if (!mCompactRepresentation)
      {
        NodeType*  p_tmp_node;
        NodeType*  p_node(pFirst());
        
        while (p_node) 
        {
          p_tmp_node = p_node->mpNext;
          free(p_node->mpLinks);
          free(p_node->mpBackLinks);
          free(p_node);
          p_node = p_tmp_node;
        }
      }
      else
      {
        NodeBasic<NodeBasicContent, _LinkContent, NODE_REGULAR, LINK_REGULAR>*  p_node(pFirst());
        NodeBasic<NodeBasicContent, _LinkContent, NODE_REGULAR, LINK_REGULAR>*  p_tmp_node;

        for(p_tmp_node =  reinterpret_cast<NodeBasic<NodeBasicContent, _LinkContent, NODE_REGULAR, LINK_REGULAR>* >(p_node); 
            p_tmp_node->NLinks() != 0; 
            p_tmp_node++) 
        {
          free(p_tmp_node->mpLinks);
        }
        
        free(p_node);
      }

      StorageType::mpFirst = NULL;
    }



  //***************************************************************************
  //***************************************************************************
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, NetworkStorageType _NetworkType, template<class> class _StorageType>
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>&
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>::
    Reverse()
    {
      NodeType*  node;
      NodeType*  p_first(StorageType::mpFirst);
      NodeType*  p_last(StorageType::mpLast);
      
      for (node = StorageType::mpFirst; node != NULL; node = node->mpBackNext) 
      {
        LinkType*  links   = node->mpLinks;
        int         nlinks  = node->NLinks();

        node->mpLinks       = node->mpBackLinks;
        node->mNLinks       = node->mNBackLinks;
        node->mpBackLinks   = links;
        node->mNBackLinks   = nlinks;
        NodeType*  next     = node->mpNext;
        node->mpNext        = node->mpBackNext;
        node->mpBackNext    = next;
        node->mAux          = -node->mAux;
        p_last              = node;
      }

      StorageType::mpFirst = p_last;
      StorageType::mpLast  = p_first;

      return *this;
    }

  //***************************************************************************
  //***************************************************************************
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, 
           NetworkStorageType _NetworkType, template<class> class _StorageType>
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>&
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>::
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
          NodeType* p_link_node = p_node->mpLinks[i].pNode();
          
          if (p_link_node->mAux == 0) 
          {
            for (j = 0; j<p_link_node->mNBackLinks && 
                        p_link_node->mpBackLinks[j].pNode()->mAux==1; j++)
            {}
            
            if (j == p_link_node->mNBackLinks) 
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
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, 
           NetworkStorageType _NetworkType, template<class> class _StorageType>
    void
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>::
    IsolateNode(NodeType* pNode)
    {
      if (pNode->pLinks())
      {
        for (size_t i(0); i< pNode->NLinks(); i++)
        {
          LinkType* p_link = &( pNode->pLinks()[i] );
          LinkType* p_back_link;
          
          if (! p_link->PointsNowhere() 
          && (NULL != (p_back_link = p_link->pNode()->pFindBackLink(pNode))))
          {
            p_back_link->Detach();
          }
        }

        // TODO: use delete
        free(pNode->pLinks());
        pNode->mpLinks = NULL;
        pNode->mNLinks = 0;
      }

      if (pNode->pBackLinks())
      {
        for (size_t i(0); i< pNode->NBackLinks(); i++)
        {
          LinkType* p_link = &( pNode->pBackLinks()[i] );
          LinkType* p_back_link;
          
          if (! p_link->PointsNowhere() 
          && (NULL != (p_back_link = p_link->pNode()->pFindLink(pNode))))
          {
            p_back_link->Detach();
          }
        }

        // TODO: use delete
        free(pNode->pBackLinks());
        pNode->mpBackLinks = NULL;
        pNode->mNBackLinks = 0;

      }
    }
  // IsolateNode
  //**************************************************************************


  //**************************************************************************
  //**************************************************************************
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, 
           NetworkStorageType _NetworkType, template<class> class _StorageType>
    void
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>::
    RemoveNode(NodeType* pNode)
    {
      LinkType* p_link;

      IsolateNode(pNode);

      erase(iterator(pNode));
    }
  // RemoveNode(NodeType* pNode);
  //**************************************************************************


  //**************************************************************************
  //**************************************************************************
  template <typename _NodeContent, typename _LinkContent, NodeRepresentationType _NodeType, LinkRepresentationType _LinkType, 
           NetworkStorageType _NetworkType, template<class> class _StorageType>
    void
    Network<_NodeContent, _LinkContent, _NodeType, _LinkType, _NetworkType, _StorageType>::
    PruneNode(NodeType* pNode)
    {
      NodeType* p_aux_node;
      LinkType* p_aux_link;

      // Check whether I am the only successor or predecessor. If yes, 
      // recursively prune forward and backward

      // forward prune
      for (size_t i = 0; i < pNode->NLinks(); i++)
      {
        p_aux_link = &( pNode->pLinks()[i] );

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
        p_aux_link = &( pNode->pBackLinks()[i] );

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

} // namespace STK

