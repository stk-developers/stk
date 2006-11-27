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

#include "Lattice.h"

namespace STK
{
    
  //***************************************************************************
  //***************************************************************************
  void 
  Lattice::
  AllocateNodesForWordLinkRecords(WordLinkRecord* pWlr, 
      NodeType*& rpNode)
  {
    // mAux had been initialized to zero. Now, it is used to count how many
    // successors has been already processed. 

    if(++pWlr->mAux < pWlr->mNReferences) 
      return;
      
    // All successors has been already processed, so make new lattice node
    // corresponding to pWlr and initialize it according to pWlr->pNode().
    Lattice::NodeType* pNode = 
      (Lattice::NodeType *) calloc(1, sizeof(Lattice::NodeType));

    if (pNode == NULL) 
      Error("Insufficient memory");
        
    switch (pWlr->pNode()->mType & (NT_WORD | NT_MODEL | NT_PHONE)) 
    {
      case NT_WORD:  pNode->mpPronun = pWlr->pNode()->mpPronun; break;
      case NT_MODEL: pNode->mpHmm    = pWlr->pNode()->mpHmm;    break;
      case NT_PHONE: pNode->mpName   = pWlr->pNode()->mpName;   break;
      default:       Error("Fatal: Invalid node type");
    }
      
    pNode->mType       = pWlr->pNode()->mType;
    pNode->SetStart(UNDEF_TIME);
    pNode->SetStop(pWlr->mTime);
    pNode->mpAlphaBeta == NULL;
      
    pNode->mNBackLinks = pWlr->mpNext    == NULL ? 0 :
                         pWlr->mpAltHyps == NULL ? 1 : pWlr->mpAltHyps->size() + 1;
                         
    // If rpNode == NULL, pWlr is the wery last WLR referenced by token in last state
    pNode->mNLinks     = rpNode == NULL ? 0 : pWlr->mNReferences;


    // Allocate space for links and backlinks      
    pNode->mpLinks     = (Lattice::LinkType *) malloc(pNode->mNLinks     * sizeof(Lattice::LinkType));
    pNode->mpBackLinks = (Lattice::LinkType *) malloc(pNode->mNBackLinks * sizeof(Lattice::LinkType));
    if (pNode->mpLinks == NULL || pNode->mpBackLinks == NULL) 
      Error("Insufficient memory");
        
    // Push new node to lattice node list
    pNode->mpNext = rpNode;
    pNode->mpBackNext = NULL;
    rpNode = pNode;
    
    if (pNode->mpNext != NULL) 
      pNode->mpNext->mpBackNext = pNode;
      
      
    // Set pWlr->pNode() to the new lattice node, so we can latter easily
    // establish  links between the nodes in the lattice nodes
    pWlr->mpNode = pNode;
      
    // All successors has been already processed, so continue recursively with Wlr's predecessors
    if(pWlr->mpNext != NULL)
    {
      AllocateNodesForWordLinkRecords(pWlr->mpNext, rpNode);
        
      if(pWlr->mpAltHyps != NULL)
      {
        for(AltHypList::iterator i = pWlr->mpAltHyps->begin(); i != pWlr->mpAltHyps->end(); i++)
          AllocateNodesForWordLinkRecords(i->mpWlr, rpNode);
      }
    }
  }
  // MakeLatticeNodesForWordLinkRecords(WordLinkRecord* pWlr, Node*& rpNode)
  //***************************************************************************

  

  //***************************************************************************
  //***************************************************************************
  void 
  Lattice::
  EstablishLinks(WordLinkRecord* pWlr)
  {
    // After calling MakeLatticeNodesForWordLinkRecords, mAux was set to pWlr->mNReferencesen.
    // Use it  to count down the successors that has been already processed. 
    if(--pWlr->mAux > 0) 
      return;
      
    // All successors have been already processed, so continue recursively with 
    // Wlr's predecessors
    if(pWlr->mpNext != NULL)
    {
      int j = 0;
      WordLinkRecord*           p                   = pWlr->mpNext;
      WordLinkRecord::LikeType  aux_total_like      = pWlr->mLike - p->mLike;
      WordLinkRecord::LikeType  aux_acoustic_like   = pWlr->mAcousticLike - p->mAcousticLike;

      p->pNode()->mpLinks[p->mNReferences - p->mAux].SetNode(         pWlr->pNode());
      p->pNode()->mpLinks[p->mNReferences - p->mAux].SetLmLike(       float_safe_substract(aux_total_like, aux_acoustic_like, 3));
      p->pNode()->mpLinks[p->mNReferences - p->mAux].SetAcousticLike( aux_acoustic_like);

      pWlr->pNode()->mpBackLinks[j].SetNode(          p->pNode());
      pWlr->pNode()->mpBackLinks[j].SetLmLike(        float_safe_substract(aux_total_like, aux_acoustic_like, 3));
      pWlr->pNode()->mpBackLinks[j].SetAcousticLike(  aux_acoustic_like);

      //printf("1 i->mLike         = %f   p->mLike         = %f   diff = %e\n",  pWlr->mLike, p->mLike, static_cast<float>(pWlr->mLike - p->mLike)); 
      //printf("1 i->mAcousticLike = %f   p->mAcousticLike = %f   diff = %e\n",  pWlr->mAcousticLike, p->mAcousticLike, static_cast<float>(pWlr->mAcousticLike - p->mAcousticLike)); 
      //printf("1 i->mLike - i->mAcousticLike = %e\n",  static_cast<float>(pWlr->mLike - pWlr->mAcousticLike)); 
      //printf("1 p->mLike - p->mAcousticLike = %e\n",  static_cast<float>(p->mLike - p->mAcousticLike)); 

      // recursively do the rest of the WLR's
      EstablishLinks(pWlr->mpNext);
        
      if(pWlr->mpAltHyps != NULL)
      {
        AltHypList::iterator i;
        
        for(j = 1, i = pWlr->mpAltHyps->begin(); i != pWlr->mpAltHyps->end(); j++, i++)
        {
          WordLinkRecord*           p                   = i->mpWlr;
          WordLinkRecord::LikeType  aux_total_like      = i->mLike - p->mLike;
          WordLinkRecord::LikeType  aux_acoustic_like   = i->mAcousticLike - p->mAcousticLike;
          
          p->pNode()->mpLinks[p->mNReferences - p->mAux].SetNode(pWlr->pNode());
          p->pNode()->mpLinks[p->mNReferences - p->mAux].SetLmLike(float_safe_substract(aux_total_like, aux_acoustic_like, 3));
          p->pNode()->mpLinks[p->mNReferences - p->mAux].SetAcousticLike(aux_acoustic_like);

          pWlr->pNode()->mpBackLinks[j].SetNode(p->pNode());
          pWlr->pNode()->mpBackLinks[j].SetLmLike(float_safe_substract(aux_total_like, aux_acoustic_like, 3));
          pWlr->pNode()->mpBackLinks[j].SetAcousticLike(aux_acoustic_like);

          //printf("n i->mLike         = %f   p->mLike         = %f   diff = %e\n",  i->mLike, p->mLike, static_cast<float>(i->mLike - p->mLike)); 
          //printf("n i->mAcousticLike = %f   p->mAcousticLike = %f   diff = %e\n",  i->mAcousticLike, p->mAcousticLike, static_cast<float>(i->mAcousticLike - p->mAcousticLike)); 
          //printf("n i->mLike - i->mAcousticLike = %e\n",  static_cast<float>(i->mLike - i->mAcousticLike)); 
          //printf("n p->mLike - p->mAcousticLike = %e\n",  static_cast<float>(p->mLike - p->mAcousticLike)); 
        
          // recursively do the rest of the WLR's
          EstablishLinks(i->mpWlr);
        }
      }
    }
  }
  // EstablishLinksBetweenLatticeNodes(WordLinkRecord* pWlr)
  //***************************************************************************

  // ************************************************************************
  // ************************************************************************
  void
  Lattice::
  BuildFromWlr(WordLinkRecord* pWlr)
  {
    AllocateNodesForWordLinkRecords(pWlr, mpFirst);
    EstablishLinks(pWlr);

    NodeType* p_node;

    for (p_node = mpFirst; NULL != p_node && NULL != p_node->mpNext;
         p_node = p_node->mpNext)
    { }

    mpLast = p_node;
  }
  // BuildFromWlr(const WordLinkRecord* pWlr);
  // ************************************************************************

  
  // ************************************************************************
  // ************************************************************************
  FLOAT
  Lattice::
  ForwardBackward()
  {
    NodeType*  p_start_node;
    bool       viterbi = true;
    FLOAT      score(0.0);

    // We assume, that the list is topologically sorted
    
    // Clear all records
    for (p_start_node = pFirst(); NULL != p_start_node; 
         p_start_node = p_start_node->mpNext) 
    {
      if (NULL != p_start_node->mpAlphaBeta)
      {
        
      }
      p_start_node->mpAlphaBeta = new AlphaBeta();
      p_start_node->mpAlphaBeta->mAlpha = LOG_0;
      p_start_node->mpAlphaBeta->mBeta = LOG_0;
    }

    pFirst()->mpAlphaBeta->mAlpha = 0.0;

    // forward direction
    for (p_start_node = pFirst(); NULL != p_start_node; 
         p_start_node = p_start_node->mpNext) 
    {
      for (size_t i = 0; i < p_start_node->mNLinks; i++)
      {
        LinkType* p_link     (&(p_start_node->mpLinks[i]));
        NodeType* p_end_node (p_link->pNode());

        score = p_start_node->mpAlphaBeta->mAlpha + p_link->Like();
        
        if (viterbi)
        {
          if (score > p_end_node->mpAlphaBeta->mAlpha)
            p_end_node->mpAlphaBeta->mAlpha = score;
        }
        else
        {
          p_end_node->mpAlphaBeta->mAlpha = LogAdd(
              p_end_node->mpAlphaBeta->mAlpha, score);
        }
      }
    }

    pLast()->mpAlphaBeta->mBeta  = 0.0;

    // backward direction
    for (p_start_node = pLast(); NULL != p_start_node; 
         p_start_node = p_start_node->mpBackNext) 
    {
      for (size_t i = 0; i < p_start_node->mNBackLinks; i++)
      {
        LinkType* p_link     (&(p_start_node->mpBackLinks[i]));
        NodeType* p_end_node (p_link->pNode());

        score = p_start_node->mpAlphaBeta->mBeta + p_link->Like();
        
        if (viterbi)
        {
          if (score > p_end_node->mpAlphaBeta->mBeta)
            p_end_node->mpAlphaBeta->mBeta = score;
        }
        else
        {
          p_end_node->mpAlphaBeta->mBeta = LogAdd(
              p_end_node->mpAlphaBeta->mBeta, score);
        }
      }
    }

    return score;
  }
  // ForwardBackward()
  // ************************************************************************



  // ************************************************************************
  // ************************************************************************
  void
  Lattice::
  PosteriorPrune(const FLOAT& thresh)
  {
    assert(NULL != pLast()->mpAlphaBeta);

    FLOAT      end_alpha = pLast()->mpAlphaBeta->mAlpha;
    NodeType*  p_start_node;

    for (p_start_node = pLast(); NULL != p_start_node;
         p_start_node = p_start_node->mpBackNext) 
    {
      // prune the node if not the end node and has no successors
      if (p_start_node != pLast() && 0 == p_start_node->NSuccessors())
      {
        IsolateNode(p_start_node);
        continue;
      }
      
      // we will count the number of pruned links, so that if all links
      // are deleted, we can drop the node
      size_t links_changed = 0;

      for (size_t i = 0; i < p_start_node->mNBackLinks; i++)
      {
        LinkType* p_link     (&(p_start_node->mpBackLinks[i]));
        NodeType* p_end_node (p_link->pNode());

        // we don't physically delete anything... just mark for deletion...
        if (p_end_node->mpAlphaBeta->mAlpha + p_start_node->mpAlphaBeta->mBeta 
            + p_link->Like() + thresh < end_alpha)
        {
          LinkType* p_reverse_link = p_end_node->pFindLink(p_start_node);

          assert(NULL != p_reverse_link);

          p_link->Detach();
          p_reverse_link->Detach();

          links_changed++;
        }

        if (links_changed == p_start_node->mNBackLinks)
        {
          // we are going to remove the node if no backlinks point to it
          //PruneNode(p_start_node);
        }
      }
    }
  };
  // PosteriorPruning()
  // ************************************************************************

};
// namespace STK
//*****************************************************************************


//#############################################################################
//# EOF                                                                      ##
//#############################################################################

