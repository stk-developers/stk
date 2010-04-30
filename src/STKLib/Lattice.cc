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
#include "mymath.h"
#include <algorithm>
#include <vector>
#include <list>

namespace STK
{
    
  //***************************************************************************
  //***************************************************************************
#if 0
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
      
    pNode->rNBackLinks() = pWlr->mpNext    == NULL ? 0 :
                         pWlr->mpAltHyps == NULL ? 1 : pWlr->mpAltHyps->size() + 1;
                         
    // If rpNode == NULL, pWlr is the wery last WLR referenced by token in last state
    pNode->rNLinks()     = rpNode == NULL ? 0 : pWlr->mNReferences;


    // Allocate space for links and backlinks      
    pNode->rpLinks()     = (Lattice::LinkType *) malloc(pNode->rNLinks()     * sizeof(Lattice::LinkType));
    pNode->rpBackLinks() = (Lattice::LinkType *) malloc(pNode->rNBackLinks() * sizeof(Lattice::LinkType));
    if (pNode->rpLinks() == NULL || pNode->rpBackLinks() == NULL) 
      Error("Insufficient memory");
        

    for (size_t ii = 0; ii < pNode->NLinks(); ++ii)
    {
      pNode->rpLinks()[ii].mAcousticLike = LOG_0;
      pNode->rpLinks()[ii].mLmLike = LOG_0;
    }

    for (size_t ii = 0; ii < pNode->NBackLinks(); ++ii)
    {
      pNode->rpBackLinks()[ii].mAcousticLike = LOG_0;
      pNode->rpBackLinks()[ii].mLmLike = LOG_0;
    }

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
  // MakeLatticeNodesForWordLinkRecords(WordLinkRecord* pWlr, NodeType*& rpNode)
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

      p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetNode(         pWlr->pNode());
      p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetLmLike(       float_safe_substract(aux_total_like, aux_acoustic_like, 3));
      p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetAcousticLike( aux_acoustic_like);

      pWlr->pNode()->rpBackLinks()[j].SetNode(          p->pNode());
      pWlr->pNode()->rpBackLinks()[j].SetLmLike(        float_safe_substract(aux_total_like, aux_acoustic_like, 3));
      pWlr->pNode()->rpBackLinks()[j].SetAcousticLike(  aux_acoustic_like);

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

          
          p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetNode(pWlr->pNode());
          p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetLmLike(float_safe_substract(aux_total_like, aux_acoustic_like, 3));
          p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetAcousticLike(aux_acoustic_like);

          pWlr->pNode()->rpBackLinks()[j].SetNode(p->pNode());
          pWlr->pNode()->rpBackLinks()[j].SetLmLike(float_safe_substract(aux_total_like, aux_acoustic_like, 3));
          pWlr->pNode()->rpBackLinks()[j].SetAcousticLike(aux_acoustic_like);

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
  // EstablishLinks(WordLinkRecord* pWlr)
  //***************************************************************************

  // ************************************************************************
  // ************************************************************************
  void
  Lattice::
  BuildFromWlr(WordLinkRecord* pWlr)
  {
    AllocateNodesForWordLinkRecords(pWlr, mpPtr);
    EstablishLinks(pWlr);

    NodeType* p_node;

    for (p_node = mpPtr; NULL != p_node && NULL != p_node->mpNext;
         p_node = p_node->mpNext)
    { }

    mpLast = p_node;
  }
#endif
  // BuildFromWlr(const WordLinkRecord* pWlr);
  // ************************************************************************

  
  // ************************************************************************
  // ************************************************************************
  FLOAT
  Lattice::
  ForwardBackward(
    FLOAT wordPenalty,
    FLOAT modelPenalty,
    FLOAT lmScale,
    FLOAT posteriorScale,
    bool viterbi)
  {
    //NodeType*  p_start_node;
    double      score(0.0);

    iterator   i_node(begin());
    iterator   i_rbegin(pLast());
    iterator   i_rend(begin()); --i_rend;

    // We assume, that the list is topologically sorted
    // Clear all records
    for (i_node = begin(); i_node != end(); ++i_node)
    {
      i_node->mC.mpAlphaBeta = new AlphaBeta();
      i_node->mC.mpAlphaBeta->mAlpha = LOG_0;
      i_node->mC.mpAlphaBeta->mBeta = LOG_0;
      i_node->mAux  = 0;
    }

    begin()->mC.mpAlphaBeta->mAlpha = 0.0;

    // forward direction
    for (i_node = begin(); i_node != end(); ++i_node)
    {
      for (int i = 0; i < i_node->NLinks(); i++)
      {
        LinkType* p_link     (&(i_node->rpLinks()[i]));
        NodeType* p_end_node (p_link->pNode());

        score = i_node->mC.mpAlphaBeta->mAlpha 
	      + (p_link->LmLike() * lmScale + p_link->AcousticLike()) * posteriorScale;
        
        if (p_end_node->mC.mType & NT_MODEL || p_end_node->mC.mType & NT_PHONE)
          score += modelPenalty * posteriorScale; 
        else if (p_end_node->mC.mType & NT_WORD &&  p_end_node->mC.mpPronun != NULL)
          score += wordPenalty * posteriorScale;
        
        if (viterbi)
        {
          if (score > p_end_node->mC.mpAlphaBeta->mAlpha)
            p_end_node->mC.mpAlphaBeta->mAlpha = score;
        }
        else
        {
          p_end_node->mC.mpAlphaBeta->mAlpha = LogAdd(
              p_end_node->mC.mpAlphaBeta->mAlpha, score);
        }
      }
    }

    pLast()->mC.mpAlphaBeta->mBeta  = 0.0;

    // backward direction
    for (i_node = i_rbegin; i_node != i_rend; --i_node)
    {
      for (int i = 0; i < i_node->NBackLinks(); i++)
      {
        LinkType* p_link     (&(i_node->rpBackLinks()[i]));
        NodeType* p_end_node (p_link->pNode());

        score = i_node->mC.mpAlphaBeta->mBeta 
	      + (p_link->LmLike() * lmScale + p_link->AcousticLike()) * posteriorScale;
        
        if (i_node->mC.mType & NT_MODEL  || i_node->mC.mType & NT_PHONE)
          score += modelPenalty * posteriorScale; 
        else if (i_node->mC.mType & NT_WORD &&  i_node->mC.mpPronun != NULL)
          score += wordPenalty * posteriorScale;


        if (viterbi)
        {
          if (score > p_end_node->mC.mpAlphaBeta->mBeta)
            p_end_node->mC.mpAlphaBeta->mBeta = score;
        }
        else
        {
          p_end_node->mC.mpAlphaBeta->mBeta = LogAdd(
              p_end_node->mC.mpAlphaBeta->mBeta, score);
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
  PosteriorExpectedCounts(
    std::map<char *,float> &countMap)
  {
    FLOAT tot_log_like = HIGHER_OF(pLast()->mC.mpAlphaBeta->mAlpha, 
                                   begin()->mC.mpAlphaBeta->mBeta);

    for (iterator i_node = begin(); i_node != end(); ++i_node) {
      char *node_name;
      
      if (i_node->mC.mType & NT_MODEL) {
        node_name = i_node->mC.mpHmm->mpMacro->mpName;
      } else if (i_node->mC.mType & NT_WORD &&  i_node->mC.mpPronun != NULL) {
        node_name = i_node->mC.mpPronun->mpWord->mpName;
      } else {
        continue;
      }

//      countMap.insert(std::make_pair(node_name, 0.0)).first->second
//        += exp(i_node->mC.mpAlphaBeta->mAlpha + i_node->mC.mpAlphaBeta->mBeta - tot_log_like);
	
      countMap[node_name] += my_exp(i_node->mC.mpAlphaBeta->mAlpha + i_node->mC.mpAlphaBeta->mBeta - tot_log_like);
	
    }
  }
  // PosteriorExpectedCounts()
  // ************************************************************************
  
  
    
  // ************************************************************************
  // ************************************************************************
  void
  Lattice::
  FreePosteriors()
  {
    for (iterator i_node = begin(); i_node != end(); ++i_node)
    {
      delete i_node->mC.mpAlphaBeta;
      i_node->mC.mpAlphaBeta = NULL;
    }
  }
  // FreePosteriors()
  // ************************************************************************



  // ************************************************************************
  // ************************************************************************
  void
  Lattice::
  PosteriorPrune(
    const FLOAT& thresh,
    FLOAT wordPenalty,
    FLOAT modelPenalty,
    FLOAT lmScale,
    FLOAT posteriorScale)
  {
    assert(NULL != pLast()->mC.mpAlphaBeta);

    FLOAT      end_alpha = pLast()->mC.mpAlphaBeta->mAlpha;
    FLOAT      start_beta = begin()->mC.mpAlphaBeta->mBeta;


    iterator   i_node(begin());
    iterator   i_rbegin(pLast());
    iterator   i_rend(begin()); --i_rend;

    while (i_node != end())
    {
      // Prune the node if it has no predecessors
      if (i_node != begin() && i_node != i_rbegin 
      &&  0 == i_node->NPredecessors())
      {
        delete i_node->mC.mpAlphaBeta;
        i_node = RemoveNode(i_node);
        continue;
      }

      // ... else go through every link and prune the links
      for (int i = 0; i < i_node->NLinks(); ++i)
      {
        LinkType* p_link = &(i_node->rpLinks()[i]);
        iterator  p_end_node(p_link->pNode());
        FLOAT     score = i_node->mC.mpAlphaBeta->mAlpha + p_end_node->mC.mpAlphaBeta->mBeta 
                  + (p_link->LmLike() * lmScale + p_link->AcousticLike() + thresh) * posteriorScale;

        if (p_end_node->mC.mType & NT_MODEL)
          score += modelPenalty * posteriorScale;
        else if (p_end_node->mC.mType & NT_WORD &&  p_end_node->mC.mpPronun != NULL)
          score += wordPenalty * posteriorScale;

        if (score < end_alpha && score < start_beta)
        {
          i_node->DetachLink(p_link);
        }
      }

      if (i_node != begin() && i_node != i_rbegin && i_node->NSuccessors() == 0)
      {
        delete i_node->mC.mpAlphaBeta;
        i_node = RemoveNode(i_node);
      }
      else
      {
        ++i_node;
      }
    }

    // In the forward pass, we might have left some nodes with no successors.
    // It is now time to remove those
    for (i_node = i_rbegin; i_node != i_rend; --i_node)
    {
      if (i_node != begin() && i_node != i_rbegin 
      && (i_node->NSuccessors() == 0))
      {
        delete i_node->mC.mpAlphaBeta;
        i_node = RemoveNode(i_node);
      }
      else
      {
        i_node->Links().defragment();
        i_node->BackLinks().defragment();
      }
    }
  };
  // PosteriorPruning()
  // ************************************************************************

  static 
  bool 
  cmpNodesByTimes(Lattice::value_type *a, Lattice::value_type *b) 
  {
     return a->mC.Start() < b->mC.Start();
  }  

    struct ActiveLinkRecord 
    {
      ActiveLinkRecord(Lattice::value_type *pn, Lattice::LinkType *pl, TimingType s, FLOAT l) : mpNode(pn), mpLink(pl), mStop(s), mLike(l) {}
      Lattice::value_type *mpNode;
      Lattice::LinkType   *mpLink;
      TimingType  mStop;
      FLOAT       mLike;
      bool operator < (const ActiveLinkRecord &b) const { return mLike > b.mLike; }
    };


  // ************************************************************************
  // ************************************************************************
  void
  Lattice::
  DensityPrune(
    const int maxDensity,
    FLOAT wordPenalty,
    FLOAT modelPenalty,
    FLOAT lmScale)
  {

    assert(NULL != pLast()->mC.mpAlphaBeta);
    
    TimingType current_time = -1;

    //FLOAT      end_alpha = pLast()->mC.mpAlphaBeta->mAlpha;
    //FLOAT      start_beta = begin()->mC.mpAlphaBeta->mBeta;

    std::vector<value_type *> nodes_by_time(this->size());
    
    
    std::vector<ActiveLinkRecord> active_links;
    
    
    for(iterator   i_node = begin(); i_node != end(); i_node++)
      nodes_by_time.push_back(& *i_node);
    
    std::sort(nodes_by_time.begin(), nodes_by_time.end(), cmpNodesByTimes);
    
    for(std::vector<value_type *>::iterator i_pnode = nodes_by_time.begin(); 
        i_pnode != nodes_by_time.end();
        i_pnode++)
    {
      if(current_time < (*i_pnode)->mC.Start())
      {
        std::vector<ActiveLinkRecord>::iterator i_alr;
        
        if(static_cast<int>(active_links.size()) > maxDensity)
        {
          nth_element(active_links.begin(), active_links.end(), active_links.begin() + maxDensity);

          for(i_alr = active_links.begin() + maxDensity; i_alr != active_links.end(); i_alr++)
            i_alr->mpNode->DetachLink(i_alr->mpLink);
          
          active_links.erase(active_links.begin() + maxDensity, active_links.end());
        }
        
        current_time = (*i_pnode)->mC.Start();
        
        for (i_alr = active_links.begin(); i_alr != active_links.end(); i_alr++)
          if (i_alr->mStop < current_time)
          {
            std::swap(*i_alr, *(active_links.end()-1));
            active_links.pop_back();
          }
      }
      
      for (int i = 0; i < (*i_pnode)->NLinks(); ++i)
      {    
        LinkType* p_link = &((*i_pnode)->rpLinks()[i]);
        if((*i_pnode)->mC.Start() != p_link->pNode()->mC.Start())
        {
          ActiveLinkRecord tmpRec(*i_pnode, p_link,
                                  p_link->pNode()->mC.Start(),
                                  (*i_pnode)->mC.mpAlphaBeta->mAlpha 
                                  + p_link->pNode()->mC.mpAlphaBeta->mBeta 
                                  + p_link->LmLike() * lmScale 
                                  + p_link->AcousticLike());
                                  
          active_links.push_back(tmpRec);
        }
      }
    }

    iterator i_node(begin());
    iterator i_rbegin(pLast());
    iterator i_rend(begin()); --i_rend;

    while (i_node != end())
    {
      // Prune the node if it has no predecessors
      if (i_node != begin() 
      && i_node != i_rbegin 
      &&  0 == i_node->NPredecessors())
      {
        delete i_node->mC.mpAlphaBeta;
        i_node = RemoveNode(i_node);
        continue;
      }

      if (i_node != begin() 
      && i_node != i_rbegin 
      && i_node->NSuccessors() == 0)
      {
        delete i_node->mC.mpAlphaBeta;
        i_node = RemoveNode(i_node);
        continue;
      }

      ++i_node;
    }

    // In the forward pass, we might have left some nodes with no successors.
    // It is now time to remove those
    for (i_node = i_rbegin; i_node != i_rend; --i_node)
    {
      if (i_node != begin() 
      && i_node != i_rbegin 
      && (i_node->NSuccessors() == 0))
      {
        delete i_node->mC.mpAlphaBeta;
        i_node = RemoveNode(i_node);
      }
      else
      {
        i_node->Links().defragment();
        i_node->BackLinks().defragment();
      }
    }
  };
  // DensityPrune()
  // ************************************************************************


};
// namespace STK
//*****************************************************************************


//#############################################################################
//# EOF                                                                      ##
//#############################################################################

