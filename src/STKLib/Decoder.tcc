#include "Lattice.h"
#include "mymath.h"
#include <map>
#include <algorithm>



namespace STK
{

  FLOAT 
  compute_diag_c_gaussian_density(
    const FLOAT*  pObs, 
    const FLOAT   gConst,
    const FLOAT*  pMean,
    const FLOAT*  pVar,
    const size_t  vSize);

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    typename Decoder<_NetworkType>::WordLinkRecord*
    Decoder<_NetworkType>::
    TimePruning(int frame_delay)
    {
      size_t              i;
      NetworkNode*  p_node;
      Token*              token = mpBestToken;
      WordLinkRecord *    twlr;
      WordLinkRecord *    rwlr = NULL;
    
      if (frame_delay > mTime-1 || !token) 
        return NULL;
    
      if (token->mpTWlr != NULL &&
        token->mpTWlr->mTime == mTime-1 - frame_delay) 
      {
        rwlr = token->mpTWlr;
      }
    
      // !!! Now, we pass only through the list of active model nodes instead of all nodes
      // !!! CODE NOT TESTED YET
      for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mC.mpAnr->mpNextActiveModel)
      {
        assert((p_node->mC.mType & NT_MODEL) && p_node->mC.mpAnr->mIsActiveModel);
    
        for (i = 0; i < p_node->mC.mpHmm->mNStates-1; i++) 
        {
          if (p_node->mC.mpAnr->mpTokens[i].IsActive()) 
          {
            Token* token = &p_node->mC.mpAnr->mpTokens[i];
    
            if (token->mpTWlr != NULL &&
              token->mpTWlr->mTime == mTime-1 - frame_delay) 
            {
              if (rwlr != token->mpTWlr) 
              {
                KillToken(token);
              } 
              else 
              {
                if (token->mpWlr == token->mpTWlr) 
                {
                  token->mpTWlr = NULL;
                } 
                else 
                {
                  for (twlr = token->mpWlr; 
                       twlr->mpNext != token->mpTWlr; 
                       twlr = twlr->mpNext)
                  {}
                  token->mpTWlr = twlr;
                }
              }
            } 
            else if (rwlr) 
            {
              KillToken(token);
            }
          }
        }
      }
      return rwlr;
    }
    
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    WriteAlpha(int time, NetworkNode* pNode, int state, Token *token)
    {
      if (pNode->mC.rpAlphaBetaListReverse() == NULL ||
        pNode->mC.rpAlphaBetaListReverse()->mTime != time) 
      {
        size_t    i;
        FWBWR *   newrec;
        
        newrec  = (FWBWR*) malloc(sizeof(FWBWR) +
                                  sizeof(newrec->mpState[0]) * (pNode->mC.mpHmm->mNStates-1));
        if (newrec == NULL) 
          Error("Insufficient memory");
          
        newrec->mpNext = pNode->mC.rpAlphaBetaListReverse();
        newrec->mTime = time;
        
        for (i=0; i<pNode->mC.mpHmm->mNStates; i++) 
        {
          newrec->mpState[i].mAlpha = newrec->mpState[i].mBeta = LOG_0;
        }
        
        pNode->mC.rpAlphaBetaListReverse() = newrec;
      }
      pNode->mC.rpAlphaBetaListReverse()->mpState[state].mAlpha = token->mLike;
      pNode->mC.rpAlphaBetaListReverse()->mpState[state].mAlphaAccuracy = token->mAccuracy;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    WriteBeta(int time, NetworkNode* pNode, int state, Token *token)
    {  
      // Record for current time must be already moved to
      // mpAlphaBetaList by function BackwardPruning
      assert(pNode->mC.rpAlphaBetaListReverse() == NULL ||
            pNode->mC.rpAlphaBetaListReverse()->mTime < time);
    
      if (pNode->mC.rpAlphaBetaList() != NULL && pNode->mC.rpAlphaBetaList()->mTime == time) {
        pNode->mC.rpAlphaBetaList()->mpState[state].mBeta = token->mLike;
        pNode->mC.rpAlphaBetaList()->mpState[state].mBetaAccuracy = token->mAccuracy;
      }
    }
  //***************************************************************************

  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    int 
    Decoder<_NetworkType>::
    BackwardPruning(int time, NetworkNode* pNode, int state)
    {
      while (pNode->mC.rpAlphaBetaListReverse() != NULL &&
            pNode->mC.rpAlphaBetaListReverse()->mTime > time) 
      {
        FWBWR *fwbwr = pNode->mC.rpAlphaBetaListReverse();
        pNode->mC.rpAlphaBetaListReverse() = fwbwr->mpNext;
        free(fwbwr);
      }
    
      if (pNode->mC.rpAlphaBetaListReverse() != NULL &&
        pNode->mC.rpAlphaBetaListReverse()->mTime == time) 
      {
        FWBWR* fwbwr = pNode->mC.rpAlphaBetaListReverse();
        pNode->mC.rpAlphaBetaListReverse() = fwbwr->mpNext;
        fwbwr->mpNext = pNode->mC.rpAlphaBetaList();
        pNode->mC.rpAlphaBetaList() = fwbwr;
      }
    
      return !( pNode->mC.rpAlphaBetaList() != NULL &&
                pNode->mC.rpAlphaBetaList()->mTime == time &&
                pNode->mC.rpAlphaBetaList()->mpState[state].mAlpha > LOG_MIN);
    }
  //***************************************************************************
  
  
  
#ifndef NDEBUG
  extern int test_for_cycle;
  extern int HasCycleCounter;
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    int 
    Decoder<_NetworkType>::
    HasCycle() 
    {
      typename NetworkType::NodeType *node,  *prev_node = NULL;
      
      HasCycleCounter++;
      
      if (!test_for_cycle) 
        return 0;
        
      for (node = mpActiveNodes; node != NULL; prev_node = node, node = node->mC.mpAnr->mpNextActiveNode) 
      {
        int     i;
        int     n_links = InForwardPass() ? node->rNLinks() : node->rNBackLinks();
        typename NetworkType::LinkType* links   = InForwardPass() ? node->rpLinks() : node->rpBackLinks();
        
        if(node->mC.mpAnr == NULL)
        {
          printf("List of active nodes contains node with mpAnr == NULL\n");
          return 1;
        }
        
        if(prev_node != node->mC.mpAnr->mpPrevActiveNode)
        {
          printf("Inconsistency in mpNextActiveNode and mpPrevActiveNode pointers\n");
          return 2;
        }
        
        if (node->mC.mAux2 == HasCycleCounter) 
        {
          printf("Cycle in list of active nodes\n");
          return 3;
        }
        
        node->mC.mAux2 = HasCycleCounter;
    
        for (i = 0; i < n_links; i++)
        {
          if (links[i].pNode()->mC.mAux2 == HasCycleCounter &&
              (!(links[i].pNode()->mC.mType & NT_MODEL)
                || links[i].pNode()->mC.mType & NT_TEE)) 
          {
            printf("Active node %d listed after his non-model succesor %d\n",
                  (int) (node - rNetwork().pFirst()), (int) (links[i].pNode() - rNetwork().pFirst()));
    
            return 4;
          }
        }
      }
      
      for (prev_node = NULL, node = mpActiveModels; node != NULL; prev_node = node, node = node->mC.mpAnr->mpNextActiveModel)
      {
        if(node->mC.mpAnr == NULL)
        {
          printf("List of active MODELS contains node with mpAnr == NULL\n");
          return 5;
        }
        
        if(prev_node != node->mC.mpAnr->mpPrevActiveModel)
        {
          printf("Inconsistency in mpNextActiveModel and mpPrevActiveModel pointers\n");
          return 6;
        }
        
        if (node->mC.mAux2 != HasCycleCounter) 
        {
          printf("Active Model is not in the list of active nodes\n");
          return 7;
        }
      } 

      return 0;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    bool
    Decoder<_NetworkType>:: 
    AllWordSuccessorsAreActive() 
    {
      typename NetworkType::NodeType* node;
      
      if (!test_for_cycle) 
        return true;
    
      for (node = mpActiveNodes; node != NULL; node = node->mC.mpAnr->mpNextActiveNode) 
      {
        int     i;
        int     n_links = InForwardPass() ? node->rNLinks() : node->rNBackLinks();
        typename NetworkType::LinkType *  links   = InForwardPass() ? node->rpLinks() : node->rpBackLinks();
    
        for (i=0; i <n_links; i++)
        {
          if (links[i].pNode()->mC.mAux2 != HasCycleCounter &&
            links[i].pNode() != (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()) &&
            (!(links[i].pNode()->mC.mType & NT_MODEL)
            || links[i].pNode()->mC.mType & NT_TEE)) 
          {
            printf("Active node %d has nonactive non-model succesor %d\n",
                   (int) (node - rNetwork().pFirst()), (int) (links[i].pNode() - rNetwork().pFirst()));
            return false;
          }
        }
      }
      return true;
    }
#endif
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    MarkWordNodesLeadingFrom(NetworkNode* node)
    {
      int       i;
      int       n_links = InForwardPass() ? node->rNLinks() : node->rNBackLinks();
      typename  NetworkType::LinkType*     p_links = InForwardPass() ? node->rpLinks() : node->rpBackLinks();
    
      for (i = 0; i < n_links; i++) 
      {
        typename NetworkType::NodeType* p_lnode = p_links[i].pNode();
        
        assert(p_lnode != node);

        if ((p_lnode->mC.mType & NT_MODEL && !(p_lnode->mC.mType & NT_TEE))
        ||  (p_lnode == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()))) 
        {
          continue;
        }
        
        if (p_lnode->mC.mpAnr == NULL) 
          p_lnode->mC.mpAnr = new ActiveNodeRecord(p_lnode);

#ifdef REPORT_TOKEN_ACTIVITY
        p_lnode->mC.mpAnr->mActivationTime = mTime;
        p_lnode->mC.mpAnr->mNodeNumber = p_lnode-pFirst();
#endif
          
        if (p_lnode->mC.mpAnr->mActiveNodeFlag > 0) 
          continue;
    
        if (p_lnode->mC.mpAnr->mIsActiveModel) 
        {
          assert(p_lnode->mC.mType & NT_TEE);
          continue;
        }
        
        // when this t mActiveNodeFlag is 
        if (p_lnode->mC.mpAnr->mActiveNodeFlag-- == 0) 
        {
          p_lnode->mC.mpAnr->mAux = 0;
          MarkWordNodesLeadingFrom(p_lnode);
        }
      }
    }
  //***************************************************************************
  
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    typename Decoder<_NetworkType>::NetworkType::NodeType *
    Decoder<_NetworkType>::
    pActivateWordNodesLeadingFrom(typename _NetworkType::NodeType* pNode)
    {
      int     i;
      int     n_links = InForwardPass() ? pNode->rNLinks() : pNode->rNBackLinks();
      typename NetworkType::LinkType*   p_links   = InForwardPass() ? pNode->rpLinks() : pNode->rpBackLinks();
    
      for (i = 0; i < n_links; i++) 
      {
        typename NetworkType::NodeType* p_lnode = p_links[i].pNode();
        
        if (((p_lnode->mC.mType & NT_MODEL) && !(p_lnode->mC.mType & NT_TEE))
        || (p_lnode == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()))) 
        {
          continue;
        }
        
        assert(p_lnode->mC.mpAnr != NULL) ;
        
        if (p_lnode->mC.mpAnr->mActiveNodeFlag++ > 0) 
          continue;
    
        if (p_lnode->mC.mpAnr->mIsActiveModel) 
        {
          assert(p_lnode->mC.mType & NT_TEE);
          continue;
        }
    
        p_lnode->mC.mpAnr->mAux++;
        
        if (p_lnode->mC.mpAnr->mActiveNodeFlag < 0) 
          continue;
    
        assert(p_lnode->mC.mpAnr->mActiveNodeFlag == 0);
        
        p_lnode->mC.mpAnr->mActiveNodeFlag  = p_lnode->mC.mpAnr->mAux;
        p_lnode->mC.mpAnr->mpNextActiveNode = pNode->mC.mpAnr->mpNextActiveNode;
        p_lnode->mC.mpAnr->mpPrevActiveNode = pNode;
        
        if (pNode->mC.mpAnr->mpNextActiveNode) 
          pNode->mC.mpAnr->mpNextActiveNode->mC.mpAnr->mpPrevActiveNode = p_lnode;
          
        pNode->mC.mpAnr->mpNextActiveNode   = p_lnode;
        pNode = pActivateWordNodesLeadingFrom(p_lnode);
      }
      
      assert(!HasCycle());
      return pNode;
    }
  // pActivateWordNodesLeadingFrom(NetworkType::NodeType* pNode)
  //***************************************************************************
  
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void
    Decoder<_NetworkType>:: 
    ActivateModel(typename NetworkType::NodeType* pNode)
    {
      if(pNode->mC.mpAnr == NULL)
      {
        pNode->mC.mpAnr = new ActiveNodeRecord(pNode);
#ifdef REPORT_TOKEN_ACTIVITY
        pNode->mC.mpAnr->mActivationTime = mTime;
        pNode->mC.mpAnr->mNodeNumber = pNode-pFirst();
#endif
      }
      else
      {
        assert(pNode->mC.mpAnr->mIsActiveModel || pNode->mC.mpAnr->mActiveNodeFlag);
      }
      
      if (pNode->mC.mpAnr->mIsActiveModel) return;
      
      pNode->mC.mpAnr->mIsActiveModel    = true;
      pNode->mC.mpAnr->mpPrevActiveModel = NULL;
      pNode->mC.mpAnr->mpNextActiveModel = mpActiveModels;
      
      if (mpActiveModels != NULL) 
      {
        mpActiveModels->mC.mpAnr->mpPrevActiveModel = pNode;
      }
      
      mpActiveModels = pNode;
    
      if (pNode->mC.mpAnr->mActiveNodeFlag) {
        assert(pNode->mC.mType & NT_TEE);
        return;
      }
    
      // probably not necessary; when removed assert on line 555 shoud be allowed
      pNode->mC.mpAnr->mActiveNodeFlag = 1; 

      pNode->mC.mpAnr->mpPrevActiveNode = NULL;
      pNode->mC.mpAnr->mpNextActiveNode = mpActiveNodes;
      
      if (mpActiveNodes != NULL) 
      {
        mpActiveNodes->mC.mpAnr->mpPrevActiveNode = pNode;
      }
      
      mpActiveNodes = pNode;
    
      assert(!HasCycle());
      MarkWordNodesLeadingFrom(pNode);
      pActivateWordNodesLeadingFrom(pNode);
      assert(AllWordSuccessorsAreActive());
    }
  // ActivateModel(NetworkType::NodeType* pNode)
  //***************************************************************************
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void
    Decoder<_NetworkType>:: 
    DeactivateWordNodesLeadingFrom(typename NetworkType::NodeType* pNode)
    {
      int       i;
      int       n_links = InForwardPass() ? pNode->rNLinks() : pNode->rNBackLinks();
      typename NetworkType::LinkType* links   
        = InForwardPass() ? pNode->rpLinks() : pNode->rpBackLinks();
      
      assert(pNode->mC.mpAnr != NULL);
    
      for (i = 0; i < n_links; i++) 
      {
        typename NetworkType::NodeType* p_lnode = links[i].pNode();
        
        if ((p_lnode->mC.mType & NT_MODEL && !(p_lnode->mC.mType & NT_TEE))
        ||  (p_lnode->mC.mType & NT_STICKY)
        ||  (p_lnode == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()))) 
        {
          continue;
        }
              
        assert(/*!(p_lnode->mC.mType & NT_TEE) ||*/ p_lnode->mC.mpAnr != NULL 
               && p_lnode->mC.mpAnr->mActiveNodeFlag > 0);
        
        if (--p_lnode->mC.mpAnr->mActiveNodeFlag) 
          continue;
    
        if (p_lnode->mC.mType & NT_TEE && p_lnode->mC.mpAnr->mIsActiveModel) 
          return;
    
        DeactivateWordNodesLeadingFrom(p_lnode);
        assert(p_lnode->mC.mpAnr->mpPrevActiveNode);
        
        p_lnode->mC.mpAnr->mpPrevActiveNode->mC.mpAnr->mpNextActiveNode = p_lnode->mC.mpAnr->mpNextActiveNode;

        if (p_lnode->mC.mpAnr->mpNextActiveNode)
          p_lnode->mC.mpAnr->mpNextActiveNode->mC.mpAnr->mpPrevActiveNode = p_lnode->mC.mpAnr->mpPrevActiveNode;
#ifdef REPORT_TOKEN_ACTIVITY
        p_lnode->mC.mpAnr->ReportActivity(mTime);
#endif
        delete p_lnode->mC.mpAnr;
        p_lnode->mC.mpAnr = NULL;
    
      }
      
      assert(!HasCycle());
    }
  // DeactivateWordNodesLeadingFrom(NetworkType::NodeType* pNode)
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void
    Decoder<_NetworkType>:: 
    DeactivateModel(typename NetworkType::NodeType* pNode)
    {
      assert(pNode->mC.mpAnr != NULL && pNode->mC.mpAnr->mIsActiveModel);
        
      pNode->mC.mpAnr->mIsActiveModel = false;
    
      if (pNode->mC.mpAnr->mpNextActiveModel != NULL) 
      {
        pNode->mC.mpAnr->mpNextActiveModel->mC.mpAnr->mpPrevActiveModel = pNode->mC.mpAnr->mpPrevActiveModel;
      }
    
      if (pNode->mC.mpAnr->mpPrevActiveModel != NULL) 
      {
        pNode->mC.mpAnr->mpPrevActiveModel->mC.mpAnr->mpNextActiveModel = pNode->mC.mpAnr->mpNextActiveModel;
      } else {
        assert(mpActiveModels == pNode);
        mpActiveModels = pNode->mC.mpAnr->mpNextActiveModel;
      }
    
      assert(!HasCycle());
      
      if (pNode->mC.mType & NT_TEE && pNode->mC.mpAnr->mActiveNodeFlag)
        return;
        
      //pNode->mpAnr->mActiveNodeFlag = 0;
  //    assert(pNode->mC.mpAnr->mActiveNodeFlag == 0);
      
      DeactivateWordNodesLeadingFrom(pNode);
    
      if (pNode->mC.mpAnr->mpNextActiveNode != NULL) {
        pNode->mC.mpAnr->mpNextActiveNode->mC.mpAnr->mpPrevActiveNode = pNode->mC.mpAnr->mpPrevActiveNode;
      }
    
      if (pNode->mC.mpAnr->mpPrevActiveNode != NULL) {
        pNode->mC.mpAnr->mpPrevActiveNode->mC.mpAnr->mpNextActiveNode = pNode->mC.mpAnr->mpNextActiveNode;
      } else {
        assert(mpActiveNodes == pNode);
        mpActiveNodes = pNode->mC.mpAnr->mpNextActiveNode;
      }      
#ifdef REPORT_TOKEN_ACTIVITY
      pNode->mC.mpAnr->ReportActivity(mTime);
#endif
      delete pNode->mC.mpAnr;
      pNode->mC.mpAnr = NULL;
    
      assert(!HasCycle());
      assert(AllWordSuccessorsAreActive());
    } // DeactivateModel
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  template<typename _Node>
    void
    Token<_Node>::
    Penalize(LikeType penalty)
    {
      mLike += penalty;
      mAccuracy.logvalue += penalty;

      if (NULL != mpAltHyps)
      {
        for (typename AltHypList::iterator it = mpAltHyps->begin(); it != mpAltHyps->end(); 
             ++it)
        {
          it->mLike += penalty;
        }
      }
    }
  //***************************************************************************

  extern int nbest_lattices;

  template<typename _Wlr> 
  bool operator < (class WlrReference<_Wlr>& a, class WlrReference<_Wlr>& b)
  {
    return a.mLike < b.mLike;
  }

  //***************************************************************************
  //***************************************************************************
  template<typename _NodeType>
    void
    Token<_NodeType>::
    AddAlternativeHypothesis(WordLinkRecord* pWlr, Token<_NodeType>::LikeType like, 
        Token<_NodeType>::LikeType acousticLike)
    { 
      WlrReference  tmp_wlr_ref(pWlr, like, acousticLike);

      if (mpAltHyps == NULL)
        mpAltHyps = new AltHypList;
      
      assert(pWlr != NULL);

      if(nbest_lattices > 1 && static_cast<int>(mpAltHyps->size()) == nbest_lattices-1) // We do not want to keep more than nbest_lattices-1 alternative hyphotesis
      {
        typename AltHypList::iterator i = min_element(mpAltHyps->begin(), mpAltHyps->end());
        if(*i < tmp_wlr_ref)
        {
          FreeWordLinkRecords(i->mpWlr);
          pWlr->mNReferences++;
          *i = tmp_wlr_ref;
        }
      }
      else
      {
        mpAltHyps->push_back(tmp_wlr_ref);
        pWlr->mNReferences++;
      }
    }
  // Token::AddAlternativeHypothesis(...)
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    int 
    Decoder<_NetworkType>::
    PassTokenMaxForLattices(Token* pFrom, Token* pTo, FLOAT totalLike, 
      FLOAT acousticLike)
    {
      int ret = 0;
      
#ifdef TRACE_TOKENS
      printf("(%.2f + %.2f -> %.2f = ", from->mLike, mLike, to->mLike);
#endif
      
      // Since we generate lattices, WLR is created for each node and therefore
      // WLR for token 'from' should have been just created
      assert(pFrom->mpAltHyps == NULL);
      
      //printf("PassTokenMaxForLattices: totalLike = %f  acousticLike = %f\n", 
      //  totalLike, acousticLike);
      
      if (!pTo->IsActive() || pFrom->mLike + totalLike > pTo->mLike) 
      {
        typename Token::AltHypList* p_alt_hyps = NULL;
        
        // move the original best WLR to the alternate hyp list (if big enough)
        if(pTo->IsActive()) 
        {
          //printf("pTo->mpWlr->mLike = %f    pTo->mLike = %f    (float)pTo->mLike = %f    diff = %f\n",
          //    pTo->mpWlr->mLike, 
          //    pTo->mLike , 
          //    (FLOAT)pTo->mLike, 
          //    pTo->mLike - pTo->mpWlr->mLike);

          pTo->AddAlternativeHypothesis(pTo->mpWlr, 
              pTo->mLike, 
              pTo->mAcousticLike);

          p_alt_hyps = pTo->mpAltHyps;
        }
        
        KillToken(pTo);
    
        ret                  = 1;
        *pTo                 = *pFrom;
        pTo->mLike          += totalLike;
        pTo->mAcousticLike  += acousticLike;
        pTo->mpAltHyps       = p_alt_hyps;
        
        if (pTo->mpWlr) 
          pTo->mpWlr->mNReferences++;
      } 
      else 
      {
        //printf("pFrom->mLike + totalLike = %f\n", pFrom->mLike + totalLike);

        pTo->AddAlternativeHypothesis(pFrom->mpWlr, 
            pFrom->mLike + totalLike, 
            pFrom->mAcousticLike + acousticLike);
      }

#ifdef TRACE_TOKENS
      printf("%.2f)\n", to->totalLike);
#endif

      return ret;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    int 
    Decoder<_NetworkType>::
    PassTokenMax(Token* pFrom, Token* pTo, FLOAT totalLike, 
      FLOAT acousticLike)
    {
      assert(pFrom->mpAltHyps == NULL);
      int ret = 0;
      
#ifdef TRACE_TOKENS
      printf("(%.2f + %.2f -> %.2f = ", pFrom->mLike, mLike, to->mLike);
#endif
      
      //printf("PassTokenMax: totalLike = %f  acousticLike = %f\n", totalLike, acousticLike);
      //  std::cout <<"from " << std::fixed << pFrom->mLike+totalLike << " to " 
      //    << std::fixed << pTo->mLike << " acoustic is " << std::fixed << acousticLike << std::endl;

      if (!pTo->IsActive() || pFrom->mLike + totalLike > pTo->mLike) 
      {
        KillToken(pTo);
    
        ret = 1;
        *pTo = *pFrom;
        pTo->mLike          += totalLike;
        pTo->mAcousticLike  += acousticLike;
        
        if (pTo->mpWlr) 
          pTo->mpWlr->mNReferences++;
      }
    #ifdef TRACE_TOKENS
      printf("%.2f)\n", pTo->mLike);
    #endif
      return ret;
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    int 
    Decoder<_NetworkType>::
    PassTokenSum(Token* pFrom, Token* pTo, FLOAT totalLike, 
      FLOAT acousticLike)
    {
      double        tl;
      double        acoustic_like;
      FloatInLog    fe;
      FloatInLog    fil_from_like = {totalLike, 0};
      int           ret = 0;
      
#ifdef TRACE_TOKENS
      printf("(%.2f + %.2f -> %.2f = ", pFrom->mLike, totalLike, pTo->mLike);
#endif
      
      if (pTo->IsActive()) 
      {
        tl = LogAdd(pTo->mLike, pFrom->mLike + totalLike);
        acoustic_like = LogAdd(pTo->mAcousticLike, pFrom->mAcousticLike + acousticLike);
        fe = FIL_Add(pTo->mAccuracy, FIL_Mul(pFrom->mAccuracy, fil_from_like));
      } 
      else 
      {
        tl = pFrom->mLike + totalLike;
        acoustic_like = pFrom->mAcousticLike + acousticLike;
        fe = FIL_Mul(pFrom->mAccuracy, fil_from_like);
      }
    
      if (!pTo->IsActive() || pFrom->mLike + totalLike > pTo->mBestLike) 
      {
        KillToken(pTo);
    
        ret = 1;
        *pTo = *pFrom;
        pTo->mBestLike = pFrom->mLike + totalLike;
        
        if (pTo->mpWlr) 
        {
          pTo->mpWlr->mNReferences++;
        }  
      }
    
      pTo->mLike          = tl;
      pTo->mAcousticLike  = acoustic_like;
      pTo->mAccuracy      = fe;
#ifdef TRACE_TOKENS
      printf("%.2f)\n", pTo->mLike);
#endif
      return ret;
    }

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void
    Decoder<_NetworkType>:: 
    Init(ModelSet* pHmms, ModelSet* pHmmsToUpdate/*, bool compactRepresentation*/) 
    {
      int                       maxStatesInModel;

      if(NULL != pHmms) {
//      mCompactRepresentation = compactRepresentation;
        PhoneNodesToModelNodes(pHmms, pHmmsToUpdate, maxStatesInModel);
    
        // Need to contain at least 2 states for propper token passing
        if (maxStatesInModel < 2) {
          Error("The network does not contain any model "
                "(possibly no dictionary suplied for word-based network)");
        }
        mpAuxTokens = new Token[maxStatesInModel-1];
        mpOutPCache = (Cache*) malloc(pHmms->mNStates      * sizeof(Cache));
        mpMixPCache = (Cache*) malloc(pHmms->mNMixtures    * sizeof(Cache));
      
        if (mpOutPCache == NULL || mpMixPCache == NULL) {
          Error("Insufficient memory");
        }
        OutputProbability =
          pHmms->mOutPdfKind == KID_DiagC     ? &DiagCGaussianMixtureDensity :
          pHmms->mOutPdfKind == KID_PDFObsVec ? &FromObservationAtStateId    : NULL;
      } else {
        mpOutPCache = NULL;
	mpMixPCache = NULL;
      }
    
      mWPenalty          = 0.0;
      mMPenalty          = 0.0;
      mPronScale         = 1.0;
      mTranScale         = 1.0;
      mOutpScale         = 1.0;
      mOcpScale          = 1.0;
      mLmScale           = 1.0;
      mLatticeGeneration = false;
      mTimePruning       = false;
      
      PassTokenInNetwork  = &PassTokenMax;
      PassTokenInModel    = &PassTokenMax;
      
      mPropagDir          = FORWARD;
      mAlignment          = WORD_ALIGNMENT;
      mpThreshState       = NULL;
      mPruningThresh      = -LOG_0;
      mMaxThreshold       = -LOG_0;
      mMaxActiveModels    = 0;
      mMinActiveModels    = 0;
      mpModelSet          = pHmms;
      mpModelSetToUpdate  = pHmmsToUpdate;
      mCollectAlphaBeta   = 0;
      mAccumType          = AT_ML;            
      mSearchPaths        = SP_ALL;

      mKeepExitToken      = false;
    }

  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    Clear()
    {
      rNetwork().Clear();
      
      delete [] mpAuxTokens;
      free(mpOutPCache);
      free(mpMixPCache);
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    KillToken(Decoder<_NetworkType>::Token* token)
    {
      token->mLike = LOG_0;
      token->mAcousticLike = LOG_0;

      FreeWordLinkRecords(token->mpWlr);
      token->mpWlr  = NULL;
      token->mpTWlr = NULL;
    }

  
  //***************************************************************************
  //*************************************************************************** 
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs, Decoder<_NetworkType>* net)  
    {
      if (!net || net->mpMixPCache[mix->mID].mTime != net->mTime)
      {
        FLOAT m_like;
        
        // we call the computation
        m_like = compute_diag_c_gaussian_density(
          pObs,
          mix->GConst(),
          mix->mpMean->mVector.cpData(),
          mix->mpVariance->mVector.cpData(),        
          mix->mpMean->VectorSize());
      
        if (net)
        {
          net->mpMixPCache[mix->mID].mTime  = net->mTime;
          net->mpMixPCache[mix->mID].mValue = m_like;
        }
        
        return m_like;
      }
      else
      {
        return net->mpMixPCache[mix->mID].mValue;
      }
    }

  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    DiagCGaussianMixtureDensity(State* state, FLOAT* pObs, Decoder* net)
    {
      size_t  i;
      FLOAT   m_like = LOG_0;
      FLOAT*  l_obs;
      
      assert(state->mOutPdfKind == KID_DiagC);
          
      //if (!net || net->mpOutPCache[state->mID].mTime != net->mTime) 
      if (net && net->mpOutPCache[state->mID].mTime == net->mTime) 
      {
        return net->mpOutPCache[state->mID].mValue;
      }
      else
      {
        for (i = 0; i < state->mNMixtures; i++) 
        {
          FLOAT    glike;
          Mixture* mix = state->mpMixture[i].mpEstimates;
      
          l_obs = XformPass(mix->mpInputXform, pObs,
                            net ? net->mTime : UNDEF_TIME,
                            net ? net->mPropagDir : FORWARD);
          
          assert(l_obs != NULL);
          glike = DiagCGaussianDensity(mix, l_obs, net);
          m_like  = LogAdd(m_like, glike + state->mpMixture[i].mWeight);
        }
      
        
#ifdef DEBUG_MSGS
        gaus_computaions++;
#endif
      
        if (net) 
        {
          net->mpOutPCache[state->mID].mTime = net->mTime;
          net->mpOutPCache[state->mID].mValue = m_like;
        }
        return m_like;
      }
    }



  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    FromObservationAtStateId(State *state, FLOAT *obs, Decoder *net)
    {

      //if (!net || net->mpOutPCache[state->mID].mTime != net->mTime) 
      if (net && net->mpOutPCache[state->mID].mTime == net->mTime) 
      {
        return net->mpOutPCache[state->mID].mValue;
      }
      else
      {
        obs = XformPass(net->mpModelSet->mpInputXform, obs,
                        net ? net->mTime : UNDEF_TIME,
                        net ? net->mPropagDir : FORWARD);
        assert(obs != NULL);

        if (net) 
        {
          net->mpOutPCache[state->mID].mTime = net->mTime;
          net->mpOutPCache[state->mID].mValue = obs[state->PDF_obs_coef]; 
        }
        return obs[state->PDF_obs_coef]; 
      }
    }

  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    ReestState(
      typename NetworkType::NodeType*   pNode,
      int           stateIndex, 
      FLOAT         logPriorProb, 
      FLOAT         updateDir,
      FLOAT*        pObs, 
      FLOAT*        pObs2) 
    {
      size_t        i;
      size_t        j;
      size_t        k;
      size_t        m;
      State*        state  = pNode->mC.mpHmm->          mpState[stateIndex];
      State*        state2 = pNode->mC.rpHmmToUpdate()->mpState[stateIndex];
      FLOAT         bjtO   = LOG_0;
      size_t        n_mixtures;
    
      if (!mpModelSetToUpdate->mGaussLvl2ModelReest 
      && ( mpModelSet != mpModelSetToUpdate))
      {
        // Occupation probabilities of mixtures are computed using target model
        state = state2; pObs = pObs2;
      } 
      else if (state->mNMixtures <= state2->mNMixtures) 
      {
        bjtO = mpOutPCache[mpModelSet->mNStates * mTime + state->mID].mValue;
      }
    
      n_mixtures = LOWER_OF(state->mNMixtures, state2->mNMixtures);
    
      if (bjtO < LOG_MIN && n_mixtures > 1)
      {
        // State likelihood was not available in cache because
        // - occupation probabilities of mixtures are computed using target model
        // - not all mixtures of mAlignment model are used for computation of
        //   occupation probabilities (state2->num_mix < state->num_mix)
        for (m = 0; m < n_mixtures; m++)
        {
          Mixture *mix  = state->mpMixture[m].mpEstimates;
          FLOAT  cjm    = state->mpMixture[m].mWeight;
          FLOAT* xobs   = XformPass(mix->mpInputXform, pObs, mTime, FORWARD);
          FLOAT  bjmtO  = DiagCGaussianDensity(mix, xobs, this);
          bjtO          = LogAdd(bjtO, cjm + bjmtO);
        }
      }
    
      //for every emitting state
      for (m = 0; m < n_mixtures; m++)
      { 
        FLOAT Lqjmt = logPriorProb;
        
        if (n_mixtures > 1)
        {
          Mixture*  mix   = state->mpMixture[m].mpEstimates;
          FLOAT     cjm   = state->mpMixture[m].mWeight;
          FLOAT*    xobs  = XformPass(mix->mpInputXform, pObs, mTime, FORWARD);
          FLOAT     bjmtO = DiagCGaussianDensity(mix, xobs, this);
          
          Lqjmt +=  -bjtO + cjm + bjmtO;
        }
    
        if (Lqjmt > MIN_LOG_WEGIHT) 
        {
          Mixture*  mix      = state2->mpMixture[m].mpEstimates;
          size_t    vec_size = mix->mpMean->VectorSize();
          FLOAT*    mnvec    = mix->mpMean->mVector.pData();
          FLOAT*    mnacc    = mix->mpMean->mpAccums;
          FLOAT*    vvacc    = mix->mpVariance->mpAccums;
          FLOAT*    vmacc    = mix->mpVariance->mpAccums + vec_size;

          XformInstance* ixf = mix->mpInputXform;
          FLOAT*    xobs     = XformPass(ixf, pObs2, mTime, FORWARD);
    
  /*      if (mmi_den_pass) {
            mnacc += vec_size + 1;
            vvacc += 2 * vec_size + 1;
            vmacc += 2 * vec_size + 1;
          }*/
    
          // the real state occupation probability
          Lqjmt = my_exp(Lqjmt) * updateDir;
    
          // Update
          for (i = 0; i < vec_size; i++) 
          {                 
            mnacc[i] += Lqjmt * xobs[i];                  // mean
            if (!(mpModelSetToUpdate->mUpdateMask & UM_OLDMEANVAR)) 
            {
              vvacc[i] += Lqjmt * SQR(xobs[i]);           // scatter
              vmacc[i] += Lqjmt * xobs[i];                // mean for var.
            } 
            else 
            {
              vvacc[i] += Lqjmt * SQR(xobs[i]-mnvec[i]);  // var
            }
          }
    
          mnacc[vec_size] += Lqjmt; //norms for mean
          vmacc[vec_size] += Lqjmt; //norms for variance
          
  //      if (mmi_den_pass) {
  //        state2->mpMixture[m].mWeightAccumDen += Lqjmt;
  //      } else {
          state2->mpMixture[m].mWeightAccum     += Lqjmt; //  Update weight accum
  //      }
    
          if (Lqjmt < 0)
            state2->mpMixture[m].mWeightAccumDen -= Lqjmt;

          // collect statistics for cluster weight vectors 
          if (mpModelSetToUpdate->mUpdateMask & UM_CWEIGHTS
          && (0 < mix->mpMean->mCwvAccum.Cols()))
          {
            for (size_t vi = 0; vi < mix->mpMean->mCwvAccum.Rows(); vi++)
            {
              mix->mpMean->mpOccProbAccums[vi] += Lqjmt;
              for (size_t vj = 0; vj < mix->mpMean->mCwvAccum.Cols(); vj++)
              {
                mix->mpMean->mCwvAccum[vi][vj] += xobs[vj] * Lqjmt;
              }
            }
          }
          
          // if Cluster Parameters are updated:
          if (0 < mix->mAccumG.Rows())
          {
            mix->mPartialAccumG += Lqjmt;
            mix->mPartialAccumK.AddCVMul(Lqjmt, xobs);
            mix->mAccumL.AddCVVDotMul(Lqjmt, xobs, vec_size, xobs, vec_size);
            
            // in case of discriminative training, we collect positive Lqjmt
            if (Lqjmt > 0.0)
            {
              mix->mPartialAccumGd += Lqjmt;
            } 
          }
      
          if (ixf == NULL || ixf->mNumberOfXformStatCaches == 0) 
            continue;
    
          UpdateXformInstanceStatCaches(ixf, pObs2, mTime);
          
          for (i = 0; i < ixf->mNumberOfXformStatCaches; i++) 
          {
            XformStatCache *xfsc = &ixf->mpXformStatCache[i];
            Variance *var  = state2->mpMixture[m].mpEstimates->mpVariance;
            Mean     *mean = state2->mpMixture[m].mpEstimates->mpMean;
    
            if (!mpModelSetToUpdate->mCmllrStats) { // compute full covariance statistics (for HLDA)
              for (j = 0; j < var->mNumberOfXformStatAccums; j++) {
                XformStatAccum *xfsa = &var->mpXformStatAccum[j];
                if (xfsa->mpXform == xfsc->mpXform) {
                  size_t size = xfsc->mpXform->mInSize;
                  for (k = 0; k < size+size*(size+1)/2; k++) {
                    xfsa->mpStats[k] += xfsc->mpStats[k] * Lqjmt;
                  }
                  xfsa->mNorm += xfsc->mNorm * Lqjmt;
                  break;
                }
              }
              for (j = 0; j < mean->mNumberOfXformStatAccums; j++) {
                XformStatAccum *xfsa = &mean->mpXformStatAccum[j];                
                if (xfsa->mpXform == xfsc->mpXform) {
                  size_t size = xfsc->mpXform->mInSize;
                  for (k = 0; k < size; k++) {
                    xfsa->mpStats[k] += xfsc->mpStats[k] * Lqjmt;
                  }
                  xfsa->mNorm += xfsc->mNorm * Lqjmt;
                  break;
                }
              }
            } else { // mpModelSet.cmllrStats
              size_t outOffset = 0;
              for (j = 0; j < mean->mNumberOfXformStatAccums; j++) {
                if (mean->mpXformStatAccum[j].mpXform == xfsc->mpXform)  {
                  if(xfsc->mpXform->mOutSize != xfsc->mpXform->mInSize - 1) {
                    Error("Invalid size of CMMLR matrix %s", xfsc->mpXform->mpMacro->mpName);
                  }
                  size_t size = xfsc->mpXform->mInSize;
                  size_t stat_size = size+size*(size+1)/2;
                  size_t dim;
                  
                  for(dim = 0; dim < size - 1; dim++) {                  
                    for (k = 0; k < size; k++) {
                      xfsc->mpXform->mpCmllrStats[stat_size * dim + k]
                        += xfsc->mpStats[k] * Lqjmt * mean->mVector[dim+outOffset] * var->mVector[dim+outOffset];
                    }
                    for (k = size; k < stat_size; k++) {
                      xfsc->mpXform->mpCmllrStats[stat_size * dim + k]
                        += xfsc->mpStats[k] * Lqjmt *                                var->mVector[dim+outOffset];
                    }
                  }
                  xfsc->mpXform->mpCmllrStats[stat_size * dim] += xfsc->mNorm * Lqjmt;
                  break;
                }
                outOffset += mean->mpXformStatAccum[j].mpXform->mOutSize;
                // !!! This is hack to deal with block diaginal CMLLR's
                // !!! Works only if all block transformatins are listed in 'xormlist' in the correct order
              }
            }
          }
        }
      }
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    TokenPropagationInit()
    {
      size_t                    i;
      typename NetworkType::NodeType*   p_node;
      typename NetworkType::iterator    i_node;
    
      InitLogMath();
   
      if (mpOutPCache != NULL) 
      {
        for (i = 0; i < mpModelSet->mNStates; i++) 
        {
          mpOutPCache[i].mTime = UNDEF_TIME;
        }
      }
    
      if (mpMixPCache != NULL) 
      {
        for (i = 0; i < mpModelSet->mNMixtures; i++) 
        {
          mpMixPCache[i].mTime = UNDEF_TIME;
        }
      }

      mpBestToken  = NULL;
      mWordThresh = LOG_MIN;
      mpActiveModels = NULL;
      mpActiveNodes = NULL;
      mNActiveTokensForUtterance = 0;
      mNActiveModelsForUtterance = 0;
      
      if (0 < mMaxActiveModels)
        mActiveModelsBestLikes.reserve(3*mMaxActiveModels/2);
    
      if (mCollectAlphaBeta && InForwardPass()) 
      {
//        if (mCompactRepresentation)
//          Error("Fatal: CSTK format used for forward-backward");
          
        for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
        {
          if (!(i_node->mC.mType & NT_MODEL)) 
            continue;
            
          i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaListReverse() = NULL;
        }
      }
      
      // Needed to load last FWBWRs to mpAlphaBetaList
      if (mCollectAlphaBeta && !InForwardPass()) 
      {
//        if (mCompactRepresentation)
//          Error("Fatal: CSTK format used for forward-backward");

        for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
        {
          if (!(i_node->mC.mType & NT_MODEL)) 
            continue;
            
          BackwardPruning(mTime, i_node.mpPtr, i_node->mC.mpHmm->mNStates-1);
        }
      }
      
      // Last p_node is not activated and deactivated in usual way; allocate ActiveNodeRecord here
      p_node = InForwardPass() ?  rNetwork().pLast() : rNetwork().pFirst();
      p_node->mC.mpAnr = new ActiveNodeRecord(p_node);

      // First p_node is also not activated and deactivated in usual way...
      p_node = InForwardPass() ? rNetwork().pFirst() : rNetwork().pLast();
      mpActiveNodes = p_node;
      p_node->mC.mpAnr = new ActiveNodeRecord(p_node);
      p_node->mC.mpAnr->mpTokens[0].mLike = 0.0;
      p_node->mC.mpAnr->mpTokens[0].mAcousticLike = 0.0;
      p_node->mC.mpAnr->mpTokens[0].mAccuracy.logvalue = LOG_0;
      p_node->mC.mpAnr->mpTokens[0].mAccuracy.negative = 0;
  //    p_node->mC.mpAnr->mpTokens[0].mpWlr = NULL; 
#ifdef bordel_staff
  //    p_node->mC.mpTokens[0].mpTWlr = NULL;
      p_node->mC.mpAnr->mpTokens[0].mBestLike = 0.0;
#endif
    
      p_node->mC.mpAnr->mpPrevActiveNode = p_node->mC.mpAnr->mpNextActiveNode = NULL;
      
      //Is this needed?
      p_node->mC.mpAnr->mActiveNodeFlag = 1;

      MarkWordNodesLeadingFrom(p_node);
      pActivateWordNodesLeadingFrom(p_node);
      TokenPropagationInNetwork();
      DeactivateWordNodesLeadingFrom(p_node);
      
      //p_node->mC.mpAnr->mActiveNodeFlag = 0;
      
      if (p_node->mC.mpAnr->mpPrevActiveNode) 
      {
        p_node->mC.mpAnr->mpPrevActiveNode->mC.mpAnr->mpNextActiveNode = p_node->mC.mpAnr->mpNextActiveNode;
      } else {
        assert(p_node == mpActiveNodes);
        mpActiveNodes = p_node->mC.mpAnr->mpNextActiveNode; 
      }

      if (p_node->mC.mpAnr->mpNextActiveNode) 
      {
        p_node->mC.mpAnr->mpNextActiveNode->mC.mpAnr->mpPrevActiveNode = p_node->mC.mpAnr->mpPrevActiveNode;
      }
    
      delete p_node->mC.mpAnr;
      p_node->mC.mpAnr = NULL;    
    } 
  // Decoder::TokenPropagationInit()
  //***************************************************************************
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    ViterbiReest(
      const Matrix<FLOAT>&  rObsMx, 
      const Matrix<FLOAT>&  rObsMx2, 
      int                   nFrames, 
      FLOAT                 weight,
      BasicVector<FLOAT>*   pWeightVector)
    {
      int                     t;
      WordLinkRecord *        wlr;
      typename NetworkType::NodeType* prevnode =        NULL;
      FLOAT                   P;
      Cache *                 p_out_p_cache;
      ModelSet *              p_hmms_alig =   mpModelSet;
      ModelSet *              p_hmms_upd =    mpModelSetToUpdate;
    
      p_out_p_cache = 
        (Cache *) malloc(nFrames * p_hmms_alig->mNStates * sizeof(Cache));
        
      if (p_out_p_cache == NULL) 
        Error("Insufficient memory");
    
      for (t = 0; t < static_cast<int>(nFrames * p_hmms_alig->mNStates); t++) 
      {
        p_out_p_cache[t].mTime  = UNDEF_TIME;
        p_out_p_cache[t].mValue = LOG_0;
      }
    
      free(mpOutPCache);
      
      mpOutPCache = NULL;
      mAlignment = STATE_ALIGNMENT;
      
      ViterbiInit();
      
      nFrames += p_hmms_alig->mTotalDelay;
    
      for (t = 0; t < nFrames; t++) 
      {
        if (t >= p_hmms_alig->mTotalDelay)
          mpOutPCache = p_out_p_cache + p_hmms_alig->mNStates*(t-p_hmms_alig->mTotalDelay);
        
        ViterbiStep(rObsMx[t]);
      }
    
      mpOutPCache = p_out_p_cache;
    
      if (rNetwork().pLast()->mC.mpAnr == NULL || !rNetwork().pLast()->mC.mpAnr->mpExitToken->IsActive()) 
      {
        ViterbiDone(NULL);
        return LOG_0;
      }
    
      p_hmms_alig->ResetXformInstances();
    
      if (p_hmms_alig != p_hmms_upd) {
        p_hmms_upd->ResetXformInstances();
      }
    
      // invert order of WRLs
      wlr = rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr;
      
      while (wlr->mpNext != NULL) 
      {
        WordLinkRecord *tmp = wlr->mpNext->mpNext;
        wlr->mpNext->mpNext = rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr;
        rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr = wlr->mpNext;
        wlr->mpNext = tmp;
      }
    
      for (mTime = -p_hmms_alig->mTotalDelay; mTime < 0; mTime++) 
      {
        FLOAT *obs          = rObsMx[mTime+p_hmms_alig->mTotalDelay];

        p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
      }
    
      if (p_hmms_alig != p_hmms_upd) 
      {
        FLOAT *obs2 = rObsMx2[mTime+p_hmms_upd->mTotalDelay];
        
        for (mTime = -p_hmms_upd->mTotalDelay; mTime < 0; mTime++) 
        {
          p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
        }
      }
    
    // Update accumulators
      for (wlr = rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr; wlr != NULL; wlr = wlr->mpNext) 
      {
        typename NetworkType::NodeType *node   = wlr->mpNode;
        int Nq       = node->mC.rpHmmToUpdate()->mNStates;
        FLOAT *aqacc = node->mC.rpHmmToUpdate()->mpTransition->mpMatrixO + SQR(Nq);
        int currstate = wlr->mStateIdx+1;
        int nextstate = (wlr->mpNext && node == wlr->mpNext->pNode())
                        ? wlr->mpNext->mStateIdx+1 : Nq-1;
        int duration  = wlr->mTime - mTime;
    
    
    
        if (prevnode != node) {
    //      if (!mmi_den_pass)
          LOG_INC(aqacc[currstate], 0 /*ln(1)*/);
          prevnode = node;
        }
    
    //    if (!mmi_den_pass) { // So far we dont do any MMI estimation of trasitions
        LOG_INC(aqacc[currstate * Nq + currstate], my_log(duration-1));
        LOG_INC(aqacc[currstate * Nq + nextstate], 0 /*ln(1)*/);
    //    }
    
        for (; mTime < wlr->mTime; mTime++) 
        {
          //for every frame of segment
          FLOAT *obs  = rObsMx [mTime+p_hmms_alig->mTotalDelay];
          FLOAT *obs2 = rObsMx2[mTime+p_hmms_upd->mTotalDelay];
          FLOAT frame_weight  = (NULL == pWeightVector) ? 1.0 :
            (*pWeightVector)[mTime + p_hmms_alig->mTotalDelay];
    
          p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
          
          if (p_hmms_alig != p_hmms_upd)
            p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
          
          ReestState(node, currstate-1, 0.0, 1.0*weight * frame_weight, obs, 
              obs2);
        }
      }
    
      P = rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr->mLike;
      //ViterbiDone(net, NULL);
      ViterbiDone(NULL);
      return P;
    }
      
      
      
    //***************************************************************************
    //***************************************************************************
    // Token section
    //***************************************************************************
    //***************************************************************************
    template <typename _NodeType>
      void 
      Token<_NodeType>::
      AddWordLinkRecord(_NodeType* pNode, int stateIndex, int time)
      {
        WordLinkRecord* wlr;
      
        // if ((wlr = (WordLinkRecord*) malloc(sizeof(WordLinkRecord))) == NULL)
        //  Error("Insufficient memory");
        wlr = new WordLinkRecord;
      
        wlr->mStateIdx      = stateIndex;
        wlr->mLike          = mLike;
        wlr->mAcousticLike  = mAcousticLike;
        wlr->mpNode         = pNode;
        wlr->mTime          = time;
        wlr->mpNext         = mpWlr;
        wlr->mNReferences   = 1;
        wlr->mAux           = 0;
        mpWlr               = wlr;

        // Lattice generation: transfer information about
        // alternative hypothesis form token to WLR
        wlr->mpAltHyps = mpAltHyps;
        mpAltHyps      = NULL;
        
#ifdef bordel_staff
        if (!mpTWlr) 
          mpTWlr = wlr;
#endif // bordel_staff
      
#ifdef DEBUG_MSGS
        wlr->mpTmpNext = firstWLR;
        firstWLR = wlr;
        wlr->mIsFreed = false;
#endif
      }
    //***************************************************************************


  //***************************************************************************
  //***************************************************************************
  template <typename _NodeType>
    Label *
    Token<_NodeType>::
    pGetLabels(bool getTimesFromNetwork, FLOAT tot_log_like)
    {
      WordLinkRecord *  wlr;
      Label *           tmp;
      Label *           level[3] = {NULL, NULL, NULL};
      int               li = 0;
        
      if (!this->IsActive()) 
        return NULL;
	
      bool poster_score = tot_log_like > LOG_MIN;
      for (wlr = this->mpWlr; wlr != NULL; wlr = wlr->mpNext) 
      {
        if ((tmp = (Label *) malloc(sizeof(Label))) == NULL)
          Error("Insufficient memory");
    
        tmp->mScore = wlr->mLike;
	
	if(poster_score) {
	  assert(wlr->pNode()->mC.mpAlphaBeta);
	  tmp->mScore = my_exp(wlr->pNode()->mC.mpAlphaBeta->mAlpha + 
                            wlr->pNode()->mC.mpAlphaBeta->mBeta - tot_log_like);
	}
        tmp->mId    = wlr->mStateIdx;

	if(getTimesFromNetwork) {
	  tmp->mStart = wlr->pNode()->mC.Start();
          tmp->mStop  = wlr->pNode()->mC.Stop();
	} else {
	  tmp->mStop  = wlr->mTime;
	} 
    
        if (wlr->pNode()->mC.mType & NT_MODEL) 
        {
          li = wlr->mStateIdx >= 0 ? 0 : 1;
          tmp->mpData = wlr->pNode()->mC.mpHmm;
          tmp->mpName = wlr->pNode()->mC.mpHmm->mpMacro->mpName;
          tmp->mpNextLevel = level[li+1] ? level[li+1] : level[2];
        } else if (wlr->pNode()->mC.mType & NT_PHONE) {
          li = 1;
          tmp->mpData = wlr->pNode()->mC.mpName;
          tmp->mpName = wlr->pNode()->mC.mpName;
          tmp->mpNextLevel = level[li+1] ? level[li+1] : level[2];
        } else //if (wlr->pNode()->mC.mpPronun->outSymbol) 
	{
          li = 2;
          tmp->mpData = wlr->pNode()->mC.mpPronun ? wlr->pNode()->mC.mpPronun->mpWord    : NULL;
          tmp->mpName = wlr->pNode()->mC.mpPronun ? wlr->pNode()->mC.mpPronun->outSymbol : NULL;
          tmp->mpNextLevel = NULL;
        } 
        //else 
        //{
        //  free(tmp);
        //  continue;
        //}
    
        if (level[li]) 
        {
    
      // if previous label mis its label on lower level, just make it
    /*      if (li > 0 && (!level[li-1] || level[li-1]->mpNextLevel != level[li])) {
            Label *tmp2;
            if ((tmp2 = (Label *) malloc(sizeof(Label))) == NULL) {
              Error("Insufficient memory");
            }
            tmp2->mpNextLevel = level[li];
            tmp2->mpName  = level[li]->mpName;
            tmp2->mScore = level[li]->mScore;
            tmp2->mStop  = level[li]->stop;
            tmp2->mId    = level[li]->mId;
    
            if (level[li-1]) {
              level[li-1]->mScore -= tmp2->mScore;
              level[li-1]->mStart  = tmp2->Stop();
            }
    
            tmp2->mpNext  = level[li-1];
            level[li-1] = tmp2;
          }*/
    
//	  if(!getTimesFromNetwork) 
	  {
            level[li]->mStart  = tmp->mStop;
	  }
	  if(!poster_score)
            level[li]->mScore -= tmp->mScore;
        }
    
        tmp->mpNext = level[li];
        level[li] = tmp;
      }
      
      for (li = 0; li < 3; li++)
      {
        if (level[li]) 
          level[li]->mStart = 0;
      }
      
      if (level[0]) level[0]->mpNextLevel = level[1];
      if (level[1]) level[1]->mpNextLevel = level[2];
      
      return level[0] ? level[0] : level[1] ? level[1] : level[2];
    }



  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    TokenPropagationInNetwork()
    {
      typename NetworkType::NodeType*     p_node;
      typename NetworkType::LinkType* links;
      int                             i;
      int                             n_links;
    
      // Beam pruning is not active in backward pass. First, it is not necessary
      // since only token that fit into the forward pass beam are allowed (backward
      // pruning after forward pass). Second, it could be dangerous. For example,
      // when training from lattices where LM weights are pushed to the begining,
      // the lowest weight (the biggest penalty) will be on the wery first link.
      // If this weight was lower than minus pruning treshold, it would result in
      // always killing token in the first node, when pruning during the backward pass.
      //                                     |
      mBeamThresh = InForwardPass() && // <--' mpBestToken &&
                    mpBestToken &&
                    mpBestToken->mLike - mPruningThresh > LOG_MIN
                    ? LOWER_OF(mpBestToken->mLike - mPruningThresh, mMaxThreshold)
                    : LOG_MIN;
      // std::cout << "Beam thresh: " << mBeamThresh << std::endl;
    
      p_node = InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst();

      assert(p_node->mC.mpAnr != NULL);
      KillToken(p_node->mC.mpAnr->mpExitToken);

  //  NetworkType::NodeType *Xnode = mpActiveNodes;
  /*  for (p_node = InForwardPass() ? rNetwork().pFirst() : rNetwork().pLast();
        p_node != NULL;
        p_node = InForwardPass() ? p_node->mC.mpNext : p_node->mC.mpBackNext)  // */
      for (p_node = mpActiveNodes; p_node != NULL; p_node = p_node->mC.mpAnr->mpNextActiveNode) 
      {
        assert(p_node->mC.mpAnr != NULL);

        // If tee model is entered
        if ((p_node->mC.mType & NT_TEE) && p_node->mC.mpAnr->mpTokens[0].IsActive()) 
        {
  //      assert(p_node->mC.mActiveNodeFlag || (p_node->mC.mType & NT_TEE && p_node->mC.mIsActive));
  //      for (Xnode = mpActiveNodes; Xnode && Xnode != p_node; Xnode = Xnode->mpNextActiveNode);
  //      assert(Xnode);
      
          if (!(mCollectAlphaBeta && !InForwardPass()
          &&   BackwardPruning(mTime, p_node, 0))) // backward pruning after forward pass
          {
            Hmm*    p_hmm       = p_node->mC.mpHmm;
            FLOAT   trans_prob  = p_hmm->mpTransition->mpMatrixO[p_hmm->mNStates - 1] *
              mTranScale;

#ifdef TRACE_TOKENS
            printf("Tee model State 0 -> Exit State ");
#endif
            if (mLatticeGeneration)
            //||  p_node->mC.mpAnr->mpTokens[0].mLike > mBeamThresh)
            {
              p_node->mC.mpAnr->mpTokens[0].AddWordLinkRecord(rNetwork().pFirst(), 
                  -1, mTime);
            }
              
            PassTokenInModel(&p_node->mC.mpAnr->mpTokens[0], 
                p_node->mC.mpAnr->mpExitToken, trans_prob, trans_prob);
          }
        }
    
        if (p_node->mC.mpAnr->mpExitToken->IsActive()) 
        {
          //FLOAT total_penalty;

  //      assert(p_node->mC.mActiveNodeFlag || (p_node->mC.mType & NT_TEE && p_node->mC.mIsActive));
  //      for (Xnode = mpActiveNodes; Xnode && Xnode != p_node; Xnode = Xnode->mpNextActiveNode);
  //      assert(Xnode);
          if (p_node->mC.mType & (NT_MODEL | NT_PHONE)) 
          {
            if (mCollectAlphaBeta) 
            {
              if (InForwardPass()) {
                WriteAlpha(mTime, p_node, p_node->mC.mpHmm->mNStates-1, p_node->mC.mpAnr->mpExitToken);
                p_node->mC.mpAnr->mpExitToken->Penalize(mMPenalty);
              } else {
                WriteBeta(mTime, p_node, 0, p_node->mC.mpAnr->mpExitToken);
              }
            }
          }
          else if (p_node->mC.mType & NT_WORD && p_node->mC.mpPronun != NULL) 
          {
            p_node->mC.mpAnr->mpExitToken->Penalize(mWPenalty + mPronScale * 
                p_node->mC.mpPronun->prob);
            
            /*if (p_node->mC.mpExitToken->mLike < mWordThresh) {
              p_node->mC.mpExitToken->mLike = LOG_0;
            }*/
          }
    
          if (p_node->mC.mpAnr->mpExitToken->mLike > mBeamThresh) 
          {
            if (p_node->mC.mType & NT_WORD 
            && (p_node->mC.mpPronun != NULL || mLatticeGeneration)
            && mAlignment & WORD_ALIGNMENT) 
            {
              p_node->mC.mpAnr->mpExitToken->AddWordLinkRecord(p_node, -1, mTime);
            } 
            else if (p_node->mC.mType & (NT_MODEL | NT_PHONE) && mAlignment & MODEL_ALIGNMENT) 
            {
              p_node->mC.mpAnr->mpExitToken->AddWordLinkRecord(p_node, -1, mTime);
            }
    
            n_links = InForwardPass() ? p_node->rNLinks() : p_node->rNBackLinks();
            links  = InForwardPass() ? p_node->rpLinks()  : p_node->rpBackLinks();
    
            for (i = 0; i < n_links; i++) 
            {
//              FLOAT acoustic_like = links[i].AcousticLike() * mTranScale;
//              FLOAT lm_like       = (links[i].LmLike() - acoustic_like) * mLmScale;

              FLOAT acoustic_like = links[i].AcousticLike(); // * mAcScale;
              FLOAT tot_like      = links[i].LmLike() * mLmScale + acoustic_like;

              if (p_node->mC.mpAnr->mpExitToken->mLike + tot_like > mBeamThresh 
                  && (((!mTimePruning
		       || ((/*links[i].pNode()->mC.Start() == UNDEF_TIME ||*/
                             links[i].pNode()->mC.Start() <= mTime) 
                       && (  links[i].pNode()->mC.Stop()  == UNDEF_TIME        
                       ||    links[i].pNode()->mC.Stop()  >= mTime)))
                  && (mSearchPaths != SP_TRUE_ONLY
                  || (p_node->mC.mType & NT_TRUE)   
                  || !(links[i].pNode()->mC.mType & (NT_MODEL | NT_PHONE))))))
              {
                if (links[i].pNode()->mC.mType & NT_MODEL) 
                {
                  ActivateModel(links[i].pNode());
                } 
                else 
                {
                  assert(links[i].pNode()->mC.mpAnr->mActiveNodeFlag ||
                         links[i].pNode() == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()));
                }

                assert(links[i].pNode()->mC.mpAnr != NULL);

#ifdef TRACE_TOKENS
                printf("Node %d -> Node %d ", p_node->mC.mpAnr->mAux, links[i].pNode()->mC.mpAnr->mAux);
#endif

                // Current lattice generation algorithm expect that word link record is created for all tokens leaving
                // any network node (including model and !NULL nodes).
  //              if (mLatticeGeneration)
  //                AddLinkToLattice(p_node, links[i].pNode(), tot_like);
                
                PassTokenInNetwork(p_node->mC.mpAnr->mpExitToken,
                    &links[i].pNode()->mC.mpAnr->mpTokens[0], tot_like, acoustic_like);
              }
            }
          } 
          else if(mLatticeGeneration)
          {
	    // Just to grab alternative hypothesis. All will be dealocated in turn.
            p_node->mC.mpAnr->mpExitToken->AddWordLinkRecord(p_node, -1, mTime);
          }
          
          if (!(p_node->mC.mType & NT_STICKY)) 
            KillToken(p_node->mC.mpAnr->mpExitToken);
        }  
      }
    
    
      if (mCollectAlphaBeta) 
      {
        if (InForwardPass()) 
        {
          for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mC.mpAnr->mpNextActiveModel) 
          {
            assert(p_node->mC.mpAnr);
            if (/*!(p_node->mC.mType & NT_MODEL) ||*/ p_node->mC.mpAnr->mpTokens[0].mLike < LOG_MIN) 
              continue;
              
            WriteAlpha(mTime, p_node, 0, &p_node->mC.mpAnr->mpTokens[0]);
          }
        } 
        else 
        {
          for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mC.mpAnr->mpNextActiveModel) 
          {
            assert(p_node->mC.mpAnr);
            if (/*!(p_node->mC.mType & NT_MODEL) ||*/ p_node->mC.mpAnr->mpTokens[0].mLike < LOG_MIN
              || BackwardPruning(mTime, p_node, p_node->mC.mpHmm->mNStates-1)) 
              continue;

            p_node->mC.mpAnr->mpTokens[0].Penalize(mMPenalty);              
            WriteBeta(mTime, p_node, p_node->mC.mpHmm->mNStates-1, &p_node->mC.mpAnr->mpTokens[0]);
          }
        }
      }
      
      if (mLatticeGeneration)
      {
        for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mC.mpAnr->mpNextActiveModel)
        {
          if (p_node->mC.mpAnr->mpTokens[0].IsActive())
          {
            p_node->mC.mpAnr->mpTokens[0].AddWordLinkRecord(rNetwork().pFirst(), -1, mTime);
	                                                 // rNetwork().pFirst() ... just to point to some !NULL node
          }
        }

        if (rNetwork().pLast()->mC.mpAnr != NULL 
        &&rNetwork().pLast()->mC.mpAnr->mpExitToken->IsActive())
        {
          rNetwork().pLast()->mC.mpAnr->mpExitToken->AddWordLinkRecord(rNetwork().pLast(), -1, mTime);
        }
      }
    
      //  Go through newly activeted models and costruct list of active nodes
      //  for (p_node=mpActiveModels; p_node&&!p_node->mC.mIsActive; p_node=p_node->mC.mpNextActiveModel) {
      //    ActivateNode(net, p_node);
      //  }
      assert(!HasCycle());
    }
  // TokenPropagationInNetwork()
  //***************************************************************************

  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void
    Decoder<_NetworkType>::
    TokenPropagationInModels(FLOAT* pObservation)
    {
      typename NetworkType::NodeType*    p_node;
      typename NetworkType::NodeType*    p_next_active;
      Hmm*                      p_hmm;
      size_t                    winingToken = 0;
      size_t                    i;
      size_t                    j;
      int                       from;
      int                       to;
      int                       state_idx;
      //  FLOAT                 threshOutProb = LOG_0;
    
      mpBestToken = NULL;
      mNActiveTokensForObservation = 0;
      mNActiveModelsForObservation = 0;
      /*  if (mpThreshState) {
        threshOutProb = OutputProbability(mpThreshState, pObservation, net);
      } */
    
      for (p_node = mpActiveModels; p_node != NULL; p_node = p_next_active) 
      {
        // Store pointer to mpNextActiveModel here, because p_node->mC.mpAnr can be dealocated in this block
        p_next_active = p_node->mC.mpAnr->mpNextActiveModel;
        
        //  for (p_node = rNetwork().pFirst(); p_node != NULL; p_node = p_node->mC.mpNext) 
        //    if (!(p_node->mC.mType & NT_MODEL)) continue;
    
        assert(p_node->mC.mpAnr);
        p_hmm = p_node->mC.mpHmm;

//printf("%s %d %d %d", p_node->mC.mpHmm->mpMacro->mpName, (int) p_node->mC.Start(), (int) p_node->mC.Stop(), (int) mTime);

        if ((mTimePruning 
	    && ((/*p_node->Start() != UNDEF_TIME &&*/p_node->mC.Start() >= mTime)
            ||  (  p_node->mC.Stop()  != UNDEF_TIME &&  p_node->mC.Stop()  <  mTime)))
            ||  (mSearchPaths == SP_TRUE_ONLY && !(p_node->mC.mType & NT_TRUE)))
        {
          for (i = 0; i < p_hmm->mNStates-1; i++) 
          {
            KillToken(&p_node->mC.mpAnr->mpTokens[i]);
          }
//printf(" Killed\n");
          DeactivateModel(p_node);
          continue;
        }
//printf("\n");

    
        if (mAccumType == AT_MPE && InForwardPass() && p_node->mC.mpAnr->mpTokens[0].IsActive()) 
        {
          FloatInLog fil_lmpa =
            {p_node->mC.mpAnr->mpTokens[0].mLike + my_log(fabs(p_node->mC.PhoneAccuracy())),
             p_node->mC.PhoneAccuracy() < 0};
    
          p_node->mC.mpAnr->mpTokens[0].mAccuracy = FIL_Add(p_node->mC.mpAnr->mpTokens[0].mAccuracy, fil_lmpa);
        }
    
        // make a copy of active tokens into mpAuxTokens
        for (i = 0; i < p_hmm->mNStates-1; i++) 
        {
          assert(!mpAuxTokens[i].IsActive());
          if (p_node->mC.mpAnr->mpTokens[i].IsActive()) 
          {
            mpAuxTokens[i] = p_node->mC.mpAnr->mpTokens[i];
            assert(mpAuxTokens[i].mpAltHyps == NULL);
            p_node->mC.mpAnr->mpTokens[i].mLike = LOG_0;
            p_node->mC.mpAnr->mpTokens[i].mAcousticLike = LOG_0;
            p_node->mC.mpAnr->mpTokens[i].mpWlr  = NULL;
          }
        }
    
        bool keep_model_active = false;
        p_node->mC.mpAnr->mpBestToken = NULL;
    
        assert(!p_node->mC.mpAnr->mpTokens[p_hmm->mNStates-1].IsActive());
    
        for (j = 1; j < p_hmm->mNStates-1; j++) 
        {
          state_idx = (InForwardPass() ? j : p_hmm->mNStates-1 - j);
    
          if (mCollectAlphaBeta 
          &&  !InForwardPass() 
          &&  BackwardPruning(mTime, p_node, state_idx)) 
          {
            continue; // backward pruning after forward pass
          }
    
          for (i = 0; i < p_hmm->mNStates-1; i++) 
          {
            from = InForwardPass() ? i : p_hmm->mNStates-1 - j;
            to   = InForwardPass() ? j : p_hmm->mNStates-1 - i;
    
            assert(!mpAuxTokens[i].IsActive() || p_node->mC.mpAnr->mIsActiveModel);
    
            if (p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] > LOG_MIN 
            &&  mpAuxTokens[i].IsActive()) 
            {
              FLOAT trans_prob = 
                p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] * mTranScale;
    
#ifdef TRACE_TOKENS
              printf("Model %d State %d -> State %d ",  p_node->mC.mpAnr->mAux, (int) i, (int) j);
#endif

              if (PassTokenInModel(&mpAuxTokens[i], &p_node->mC.mpAnr->mpTokens[j], 
                    trans_prob, trans_prob)) 
              {
                winingToken = i;
              }
            }
          }
    
          // if (IS_ACTIVE(p_node->mC.mpTokens[j])) 
          // std::cout << "thresh: " << mBeamThresh << std::endl;
          if (p_node->mC.mpAnr->mpTokens[j].mLike > mBeamThresh) 
          {
            FLOAT out_prob = OutputProbability(p_hmm->mpState[state_idx-1],
                                              pObservation, this);
            // std::cout << "out_prob: " << out_prob << std::endl;
            
            out_prob *= mOutpScale;
    
            /*if (out_prob < threshOutProb) {
              out_prob = threshOutProb;
            }*/
    
            if (mCollectAlphaBeta && !InForwardPass())
              WriteBeta(mTime, p_node, state_idx, &p_node->mC.mpAnr->mpTokens[j]);
            
            if (mAccumType == AT_MFE && p_node->mC.mType & NT_TRUE) 
            {
              FloatInLog fil_like = {p_node->mC.mpAnr->mpTokens[j].mLike, 0};
              p_node->mC.mpAnr->mpTokens[j].mAccuracy = FIL_Add(p_node->mC.mpAnr->mpTokens[j].mAccuracy, fil_like);
            }
    
            // std::cout << "old like: " << p_node->mC.mpAnr->mpTokens[j].mLike  << std::endl;
            p_node->mC.mpAnr->mpTokens[j].mAccuracy.logvalue += out_prob;
            p_node->mC.mpAnr->mpTokens[j].mLike              += out_prob;
            p_node->mC.mpAnr->mpTokens[j].mAcousticLike      += out_prob;
    
            if (mCollectAlphaBeta && InForwardPass()) 
              WriteAlpha(mTime, p_node, state_idx, &p_node->mC.mpAnr->mpTokens[j]);
    
            if (mAlignment & STATE_ALIGNMENT && winingToken > 0 &&
              (winingToken != j || mAlignment & FRAME_ALIGNMENT)) 
            {
              p_node->mC.mpAnr->mpTokens[j].AddWordLinkRecord(
                  p_node,
                  (InForwardPass() ? winingToken : p_hmm->mNStates-1 - winingToken)-1,
                  mTime-1);            
            }
    
            keep_model_active = true;
            mNActiveTokensForObservation++;
            assert(p_node->mC.mpAnr->mIsActiveModel);
            
            if (NULL == p_node->mC.mpAnr->mpBestToken
            ||  p_node->mC.mpAnr->mpBestToken->mLike < p_node->mC.mpAnr->mpTokens[j].mLike)
            {
              p_node->mC.mpAnr->mpBestToken = &p_node->mC.mpAnr->mpTokens[j];
            }
          } 
          else 
          {
            assert(p_node->mC.mpAnr->mIsActiveModel || !p_node->mC.mpAnr->mpTokens[j].IsActive());
            KillToken(&p_node->mC.mpAnr->mpTokens[j]);
          }
        }
    
        for (i = 0; i < p_hmm->mNStates-1; i++) 
        {
          KillToken(&mpAuxTokens[i]);
        }
    
        if (!keep_model_active)
        {
          DeactivateModel(p_node);
          continue;
        }
        
        mNActiveModelsForObservation++;
        if (NULL == mpBestToken || mpBestToken->mLike < p_node->mC.mpAnr->mpBestToken->mLike) 
        {
          mpBestToken = p_node->mC.mpAnr->mpBestToken;
          mpBestNode  = p_node;
        }
        
        // Store likelihoods of best tokens of all active models to allow for
        // pruning (see begining of ), where only n-best models are kept active
        if (0 < mMaxActiveModels || 0 < mMinActiveModels)
          mActiveModelsBestLikes.push_back(p_node->mC.mpAnr->mpBestToken->mLike);

    
        if (mCollectAlphaBeta && !InForwardPass() &&
            BackwardPruning(mTime-1, p_node, 0))
        {
          // according to backward prunning (after forward pass) we should not
          // propagate token to first model's state, so we stop dealing with
          // the curent model here. The code below, must deal only with the problems
          // related to propagating token to last/first non-emiting state (in the
          // case of forward/backward propagation).
          continue;
        }
    
        // Propagate tokens form emmiting states to last/first non-emiting state        
        for (i = 1; i < p_hmm->mNStates-1; i++) 
        {
          from = InForwardPass() ? i : 0;
          to   = InForwardPass() ? p_hmm->mNStates-1 : p_hmm->mNStates-1 - i;
    
          if (p_node->mC.mpAnr->mpTokens[i].IsActive()) 
          {
            if (p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] > LOG_MIN) 
            {
              FLOAT trans_prob = p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] *
                mTranScale;
#ifdef TRACE_TOKENS
              printf("Model %d State %d -> Exit State ",  p_node->mC.mpAnr->mAux, (int) i);
#endif
              if (PassTokenInModel(&p_node->mC.mpAnr->mpTokens[i],
                                   &p_node->mC.mpAnr->mpTokens[p_hmm->mNStates - 1],
                                   trans_prob, trans_prob)) 
              {
                winingToken = i;
              }
            }
          }
        }
    
        if (p_node->mC.mpAnr->mpTokens[p_hmm->mNStates - 1].IsActive()) 
        {
          if (mAccumType == AT_MPE && !InForwardPass()) 
          {
            FloatInLog fil_lmpa =
              {p_node->mC.mpAnr->mpTokens[p_hmm->mNStates - 1].mLike + my_log(fabs(p_node->mC.PhoneAccuracy())),
               p_node->mC.PhoneAccuracy() < 0};
    
            p_node->mC.mpAnr->mpTokens[p_hmm->mNStates - 1].mAccuracy =
              FIL_Add(p_node->mC.mpAnr->mpTokens[p_hmm->mNStates - 1].mAccuracy, fil_lmpa);
          }
    
          if (mAlignment & STATE_ALIGNMENT) 
          {
            p_node->mC.mpAnr->mpTokens[p_hmm->mNStates - 1].AddWordLinkRecord(
                p_node,
                (InForwardPass() ? winingToken : p_hmm->mNStates-1 - winingToken-1)-1,
                mTime);
          }
        }
      }

      assert(!HasCycle());
//printf("%d: %d\n", mTime, mNActiveModelsForObservation);
      mNActiveTokensForUtterance += mNActiveTokensForObservation;
      mNActiveModelsForUtterance += mNActiveModelsForObservation;
      
      // Do pruning, where only mMaxActiveModels models with the highest
      // likelihood tokens are left active.
      if (0 < mMaxActiveModels && mActiveModelsBestLikes.size() > static_cast<size_t>(mMaxActiveModels))
      {
        assert(mActiveModelsBestLikes.size() == static_cast<size_t>(mNActiveModelsForObservation));

        // Find n-th best likelihood in the vector
        std::nth_element(mActiveModelsBestLikes.begin(),
                    mActiveModelsBestLikes.end() - mMaxActiveModels,
                    mActiveModelsBestLikes.end());
                    
        FLOAT pruning_threshold = *(mActiveModelsBestLikes.end() - mMaxActiveModels);
                        
        for (p_node = mpActiveModels; p_node != NULL; p_node = p_next_active)
        {
          p_next_active = p_node->mC.mpAnr->mpNextActiveModel;
          
          if(p_node->mC.mpAnr->mpBestToken->mLike < pruning_threshold)
          {
            for (i = 1; i < p_node->mC.mpHmm->mNStates; i++)
              KillToken(&p_node->mC.mpAnr->mpTokens[i]);
          
            DeactivateModel(p_node);
          }
        }
      }
      
      mMaxThreshold = -LOG_0;
      
      if (0 < mMinActiveModels)
      {
        assert(mActiveModelsBestLikes.size() == static_cast<size_t>(mNActiveModelsForObservation));

        if(mActiveModelsBestLikes.size() < static_cast<size_t>(mMinActiveModels))
          mMaxThreshold = LOG_MIN;
        else
        {
          // Find n-th best likelihood in the vector
          nth_element(mActiveModelsBestLikes.begin(),
                      mActiveModelsBestLikes.end() - mMinActiveModels,
                      mActiveModelsBestLikes.end());
                    
          mMaxThreshold = *(mActiveModelsBestLikes.end() - mMinActiveModels);
        }
      }
      
      mActiveModelsBestLikes.clear();    
      // std::cout << "Best token like: " << mpBestToken->mLike << std::endl;
    }
  // TokenPropagationInModels(FLOAT* pObservation)
  //***************************************************************************
  
  

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    TokenPropagationDone()
    // {{{
    {
      size_t  j;
      typename NetworkType::NodeType*   node = InForwardPass() 
        ? rNetwork().pLast() 
        : rNetwork().pFirst();
    
      // Last node is not activated and deactivated in usual way; deallocate ActiveNodeRecord here
      assert(node->mC.mpAnr != NULL);    
      KillToken(node->mC.mpAnr->mpExitToken);
      delete node->mC.mpAnr;
      node->mC.mpAnr = NULL;
    
      while(mpActiveModels != NULL) 
      {
        for (j=0; j < mpActiveModels->mC.mpHmm->mNStates; j++) 
        {
          KillToken(&mpActiveModels->mC.mpAnr->mpTokens[j]);
        }
        
        DeactivateModel(mpActiveModels);
      }
      
//#ifndef NDEBUG
      for (typename NetworkType::iterator i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node)
      {
        delete i_node->mC.mpAnr;
        i_node->mC.mpAnr = NULL;
        assert(i_node->mC.mpAnr == NULL);
      }
//#endif      
  //    if(mLatticeGeneration)
  //      my_hdestroy_r(&mLatticeNodeHash, 1);
  //
    } // TokenPropagationDone() }}}                 
  //***************************************************************************

    

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    FreeFWBWRecords()
    {
      typename NetworkType::iterator i_node;
      
      // go through all nodes in the network and free the alpha/beta list
      // and alpha/beta reverse list
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); 
           ++i_node) 
      {
        if (!(i_node->mC.mType & NT_MODEL)) 
          continue;
    
        while (i_node->mC.rpAlphaBetaList()) 
        {
          FWBWR *fwbwr = i_node->mC.rpAlphaBetaList();
          i_node->mC.rpAlphaBetaList() = fwbwr->mpNext;
          free(fwbwr);
        }
    
        while (i_node->mC.rpAlphaBetaListReverse()) 
        {
          FWBWR *fwbwr = i_node->mC.rpAlphaBetaListReverse();
          i_node->mC.rpAlphaBetaListReverse() = fwbwr->mpNext;
          free(fwbwr);
        }
      }
    }



  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    typename Decoder<_NetworkType>::NetworkType::NodeType* 
    Decoder<_NetworkType>:: 
    PhoneNodesToModelNodes(ModelSet * pHmms, ModelSet *pHmmsToUpdate, int& maxStatesInModel)
    {
      typename NetworkType::iterator p_node;
      typename NetworkType::iterator p_last_node;
    
//      if (pHmmsToUpdate == NULL) 
//        pHmmsToUpdate = pHmms;

      mNumberOfNetStates = 0;
      maxStatesInModel   = 0;
    
      // for (p_node =  rNetwork().begin(); 
      // p_node != rNetwork().end(); 
      // p_last_node = p_node, p_node = mCompactRepresentation 
      // ? (p_node->mC.rNLinks() ? reinterpret_cast<NetworkType::NodeType*>(reinterpret_cast<NodeBasic<NODE_REGULAR, LINK_REGULAR>* >(p_node)+1) : NULL)
      // : p_node->mC.mpNext)
      
      for (p_node =  rNetwork().begin(); 
           p_node != rNetwork().end(); 
           p_last_node = p_node, ++p_node)
      {
        if (p_node->mC.mType & NT_PHONE)  
        {
          Macro *macro;
	  char* phone_name;
          p_node->mC.mType &= ~NT_PHONE;
          p_node->mC.mType |= NT_MODEL;
	  phone_name = p_node->mC.mpName; // mpName is in union -> will be rewritten by mpHmm
          macro = FindMacro(&pHmms->mHmmHash, phone_name);
      
          if (macro == NULL) 
          {
            Error("Model %s not defined in %sHMM set", phone_name,
                  pHmmsToUpdate != pHmms ? "alignment " : "");
          }
          p_node->mC.mpHmm =  (Hmm *) macro->mpData;
          
          if (NULL != pHmmsToUpdate) 
          {
            if (pHmmsToUpdate != pHmms) 
            {
              macro = FindMacro(&pHmmsToUpdate->mHmmHash, phone_name);
              if (macro == NULL) {
                Error("Model %s not defined in HMM set", phone_name,
                      pHmmsToUpdate != pHmms ? "" : "target ");
              }
              p_node->mC.rpHmmToUpdate() = (Hmm *) macro->mpData;
            } 
            else
            {
              p_node->mC.rpHmmToUpdate() = p_node->mC.mpHmm;
            }
            p_node->mC.rpHmmToUpdate()->mpMacro->mOccurances++;
          }
          
          int nstates = p_node->mC.mpHmm->mNStates;
          
          if (p_node->mC.mpHmm->mpTransition->mpMatrixO[nstates - 1] > LOG_MIN) 
          {
            p_node->mC.mType |= NT_TEE;
          }
#ifndef NDEBUG
          p_node->mC.mAux2 = 0;
          p_node->mC.mEmittingStateId = mNumberOfNetStates;
#endif      
          if (maxStatesInModel < nstates) 
            maxStatesInModel = nstates;
          
          assert(nstates >= 2); // two non-emiting states
          mNumberOfNetStates += nstates;
        }
      }

      return &(*p_last_node);
    }

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void 
    Decoder<_NetworkType>::
    ViterbiInit()
    {
      PassTokenInModel    = &PassTokenMax;
      PassTokenInNetwork  = mLatticeGeneration ? &PassTokenMaxForLattices : &PassTokenMax;
      mPropagDir          = FORWARD;

      if(NULL != mpModelSet)
        mpModelSet->ResetXformInstances();
    
      mTime = 0; // Must not be set to -mpModelSet->mTotalDelay yet
                 // otherwise token cannot enter first model node
                 // with start set to 0
                 
      TokenPropagationInit();
      
      if(NULL != mpModelSet)
        mTime = -mpModelSet->mTotalDelay;
    }  
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    void
    Decoder<_NetworkType>:: 
    ViterbiStep(FLOAT* pObservation)
    {
      mTime++;
      mpModelSet->UpdateStacks(pObservation, mTime, mPropagDir);
    
      if (mTime <= 0)
        return;
    
      TokenPropagationInModels(pObservation);
      TokenPropagationInNetwork();
    #ifdef DEBUG_MSGS
      printf("Frame: %ld Nuberm of WordLinkRecord Active: ", mTime); PrintNumOfWLR();
    #endif
    }
  
  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    ViterbiDone(Label** pLabels, Lattice* pLattice, bool getTimesFromNetwork)
    {
      FLOAT tot_like = LOG_0;
      
      if (pLabels  != NULL) *pLabels  = NULL;
      if (pLattice != NULL) pLattice->Clear();
      
      // Although this peace of code is also in TokenPropagationDone(),
      // we must make sure that no token other then the token in the last node
      // is active before calling pGetLattice(). This is because pGetLattice() rely
      // on mNReferences and expect that (except the last node) WLR is referenced
      // only by other WRLs and not by any token.
      
      // !!! It would be better to rewrite pGetLattice() so that it do not depend on mNReferences
      
      while(mpActiveModels != NULL)
      {
        for (size_t j=0; j < mpActiveModels->mC.mpHmm->mNStates; j++)
          KillToken(&mpActiveModels->mC.mpAnr->mpTokens[j]);
          ;
        DeactivateModel(mpActiveModels);
      }

      if (rNetwork().pLast()->mC.mpAnr &&
          rNetwork().pLast()->mC.mpAnr->mpExitToken && 
          rNetwork().pLast()->mC.mpAnr->mpExitToken->IsActive()) 
      {
        tot_like = rNetwork().pLast()->mC.mpAnr->mpExitToken->mLike;

        if (pLabels != NULL)
          *pLabels = rNetwork().pLast()->mC.mpAnr->mpExitToken->pGetLabels(getTimesFromNetwork, 
	                                             rNetwork().pFirst()->mC.mpAlphaBeta 
						     ? rNetwork().pFirst()->mC.mpAlphaBeta->mBeta
						     : LOG_0);
          
        if (/*mLatticeGeneration && */pLattice != NULL)
        {
          //pLattice->BuildFromWlr(rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr);
          Wlr2Lattice(rNetwork().pLast()->mC.mpAnr->mpExitToken->mpWlr, *pLattice);
        }
      }
      
      TokenPropagationDone();
    
    #ifdef DEBUG_MSGS
      printf("Number of WordLinkRecord Unreleased: "); PrintNumOfWLR();
      printf("Number of output prob. computations: %d\n", gaus_computaions);
    #endif
      return tot_like;
    }

  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    typename Decoder<_NetworkType>::FWBWRet
    Decoder<_NetworkType>::
    ForwardBackward(const Matrix<FLOAT>& rFeatureMatrix, size_t nFrames)
    {
      int         i;
      Cache*      p_out_p_cache;
      FWBWRet     ret;
      ModelSet*   hmms = mpModelSet;
      
      p_out_p_cache = (Cache*) malloc(nFrames * hmms->mNStates * sizeof(Cache));
    
      if (p_out_p_cache == NULL) 
      {
        Error("Insufficient memory (nFrames=%d, nStates=%d)", nFrames, hmms->mNStates);
      }
    
      // clear the cache
      for (i = 0; i < static_cast<int>(nFrames * hmms->mNStates); i++) 
      {
        p_out_p_cache[i].mTime  = UNDEF_TIME;
        p_out_p_cache[i].mValue = LOG_0;
      }
    
      free(mpOutPCache);
      mpOutPCache = NULL;
    
      PassTokenInModel    = &PassTokenSum;
      PassTokenInNetwork  = &PassTokenSum;
      mAlignment          = NO_ALIGNMENT;
    
      //Forward Pass
      mPropagDir          = FORWARD;
      mCollectAlphaBeta   = 1;
      mTime               = -hmms->mTotalDelay;
    
      mpModelSet->ResetXformInstances();
      
      for (i = 0; i < hmms->mTotalDelay; i++) 
      {
        mTime++;
        hmms->UpdateStacks(rFeatureMatrix[i], mTime, FORWARD);
      }
    
      TokenPropagationInit();
    
      for (i = hmms->mTotalDelay; 
           i < static_cast<int>(nFrames) + hmms->mTotalDelay; 
           i++) 
      {
        mTime++;
        mpOutPCache = p_out_p_cache + hmms->mNStates * (mTime-1);
    
        mpModelSet->UpdateStacks(rFeatureMatrix[i], mTime, mPropagDir);
    
        TokenPropagationInModels(rFeatureMatrix[i]);
        TokenPropagationInNetwork();
      }
    
      if (rNetwork().pLast()->mC.mpAnr == NULL 
      || !rNetwork().pLast()->mC.mpAnr->mpExitToken->IsActive())  
      { // No token survivered
        TokenPropagationDone();
        FreeFWBWRecords();
        mpOutPCache = p_out_p_cache;
        ret.totLike = LOG_0;
        return ret;
      }
    
      ret.totLike = rNetwork().pLast()->mC.mpAnr->mpExitToken->mLike; //  totalLikelihood;
      TokenPropagationDone();
    
      //Backward Pass
      mPropagDir = BACKWARD;
      mTime      = nFrames + hmms->mTotalDelay;
                            
      for (i = nFrames + hmms->mTotalDelay - 1; i >= static_cast<int>(nFrames); i--) 
      {
        //  We do not need any features, in backward prop. All output probab. are 
        //  cached.
        //  UpdateStacks(hmms, obsMx + hmms->mInputVectorSize * i,
        //               mTime, BACKWARD);
        mTime--;
      }
    
      mpOutPCache = NULL;  // TokenPropagationInit would reset valid 1st frm cache
      TokenPropagationInit();
    
      for (i = nFrames-1; i >= 0; i--) 
      {
        mpOutPCache = p_out_p_cache + hmms->mNStates * (mTime-1);
    
        //  We do not need any features, in backward prop. All output probab. are 
        //  cached.
        //  UpdateStacks(mpModelSet, obsMx + hmms->mInputVectorSize * i,     |
        //                             mTime, mPropagDir);                   |
        //                                                                   V
        TokenPropagationInModels(NULL);   //obsMx + hmms->mInputVectorSize * i);
        mTime--;
        TokenPropagationInNetwork();
      }
    
      mpOutPCache = p_out_p_cache;
    
      if (rNetwork().pFirst()->mC.mpAnr == NULL || !rNetwork().pFirst()->mC.mpAnr->mpExitToken->IsActive()) 
      { // No token survivered
        TokenPropagationDone();
        FreeFWBWRecords();
        ret.totLike = LOG_0;
        return ret;
      }
      ret.totLike = HIGHER_OF(ret.totLike, rNetwork().pFirst()->mC.mpAnr->mpExitToken->mLike); //  totalLikelihood;
      // Backward pass P can differ from forward pass P because of the precision
      // problems. Take the higher one to decrease the possibility of getting
      // an occupation probability (when normalizing by P) higher that one.
    
      FloatInLog fil_ret_totLike = {ret.totLike, 0};
      ret.avgAccuracy  = FIL_Div(rNetwork().pFirst()->mC.mpAnr->mpExitToken->mAccuracy, fil_ret_totLike);
      TokenPropagationDone();
    
      // There may be remaining records in mpAlphaBetaListReverse unused in
      // backward pass. Free them and set mpAlphaBetaListReverse's to NULLs;
      typename NetworkType::iterator i_node;
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (!(i_node->mC.mType & NT_MODEL)) 
          continue;
    
        while (i_node->mC.rpAlphaBetaListReverse()) 
        {
          FWBWR* fwbwr = i_node->mC.rpAlphaBetaListReverse();
          i_node->mC.rpAlphaBetaListReverse() = fwbwr->mpNext;
          free(fwbwr);
        }
      }
    
      return ret;
    }
  // ForwardBackward(const Matrix<FLOAT>& rFeatureMatrix, size_t nFrames)
  //***************************************************************************

  
  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    MCEReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, 
           int nFrames, FLOAT weight, FLOAT sigSlope)
    {
      struct FWBWRet    fwbw;
      FLOAT             TP;
      FLOAT             P;
      FLOAT             F;
      
      int               i;
      int               j;
      int               k;
      int               t;
      
      ModelSet*         p_hmms_alig = mpModelSet;
      ModelSet*         p_hmms_upd = mpModelSetToUpdate;
      //typename NetworkType::NodeType*             node;
      typename NetworkType::iterator              i_node;
    
      AccumType origAccumType = mAccumType;
      
      mAccumType          = AT_ML;
      mPropagDir          = FORWARD;
      mAlignment          = NO_ALIGNMENT;
      
      // pointers to functions
      PassTokenInModel    = &PassTokenSum;  
      PassTokenInNetwork  = &PassTokenSum;
    
      mSearchPaths        = SP_TRUE_ONLY;
      mpModelSet->ResetXformInstances();
    
      mTime = 0; // Must not be set to -mpModelSet->totalDelay yet
                 // otherwise token cannot enter first model node
                 // with start set to 0
                 
      TokenPropagationInit();
      
      mTime = -mpModelSet->mTotalDelay;
    
      for (t = 0; t < nFrames + p_hmms_alig->mTotalDelay; t++) 
      {
        ViterbiStep(rObsMx[t]);
      }
    
      TP = rNetwork().pLast()->mC.mpAnr == NULL ? LOG_0 : rNetwork().pLast()->mC.mpAnr->mpExitToken->mLike;
      
      ViterbiDone(NULL);
    
      if (TP <= LOG_MIN) 
        return LOG_0;
    
    
      ////////////////// Denominator accumulation //////////////////
      mSearchPaths = SP_ALL;
    
      fwbw = ForwardBackward(rObsMx, nFrames);
      P = fwbw.totLike;
    
      assert(P >= TP);
      if(sigSlope > 0.0) 
      {
        F = TP - LogSub(P, TP);
        printf("MCE distance: %g; ", F);
        F = my_exp(-sigSlope * F);
        F = (sigSlope*F) / SQR(1+F);
        printf("weight: %g\n", F);
        weight *= F;
      }

      if (P < LOG_MIN) return LOG_0;
    
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.mType & NT_MODEL && i_node->mC.rpAlphaBetaList() != NULL &&
          i_node->mC.rpAlphaBetaList()->mTime == 0) 
        {
          i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
          i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
        }
      }
    
      p_hmms_alig->ResetXformInstances();
    
      if (p_hmms_alig != p_hmms_upd) {
        p_hmms_upd->ResetXformInstances();
      }
    
      for (i = 0; i < p_hmms_alig->mTotalDelay; i++) 
      {
        p_hmms_alig->UpdateStacks(rObsMx[i],
                    i-p_hmms_alig->mTotalDelay, FORWARD);
      }
    
      if (p_hmms_alig != p_hmms_upd) 
      {
        for (i = 0; i < p_hmms_upd->mTotalDelay; i++) 
        {
          p_hmms_upd->UpdateStacks(rObsMx2[i],
                      i - p_hmms_upd->mTotalDelay, FORWARD);
        }
      }
    
    
      // mpMixPCache might be used to cache likelihoods of mixtures of target
      // models. Reallocate the cache to fit mixtures of both models and reset it.
      k = HIGHER_OF(p_hmms_upd->mNMixtures, p_hmms_alig->mNMixtures);
      mpMixPCache = (Cache*) realloc(mpMixPCache, k * sizeof(Cache));
      if (mpMixPCache == NULL) Error("Insufficient memory");
    
      for (i = 0; i < k; i++)
      {
        mpMixPCache[i].mTime = UNDEF_TIME;
      }
    
      // Update accumulators
      for (mTime = 0; mTime < nFrames; mTime++) 
      { // for every frame
        FLOAT* obs  =rObsMx[mTime+p_hmms_alig->mTotalDelay];
        FLOAT* obs2 =rObsMx2[mTime+p_hmms_upd->mTotalDelay];
        
        p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
        
        if (p_hmms_alig != p_hmms_upd) {
          p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
        }
    
        for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
        { //for every model
    //    for (k=0; k < nnodes; k++) 
    //      NetworkType::NodeType *node = &mpNodes[k];
          if (i_node->mC.mType & NT_MODEL &&
            i_node->mC.rpAlphaBetaList() != NULL &&
            i_node->mC.rpAlphaBetaList()->mTime == mTime+1) {
    
            struct AlphaBetaMPE *st;
            int Nq       = i_node->mC.mpHmm->mNStates;
            st = i_node->mC.rpAlphaBetaList()->mpState;
    
            for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
              if (st[j].mAlpha + st[j].mBeta - P > MIN_LOG_WEGIHT) {
                assert(i_node->mC.rpAlphaBetaListReverse()->mTime == mTime);
    
                ReestState(
                    &(*i_node),
                    j-1,
                    (st[j].mAlpha + st[j].mBeta - P)  * mOcpScale,
                    -weight,
                    obs,
                    obs2);
              }
            }
            
            if (i_node->mC.rpAlphaBetaListReverse()) 
              free(i_node->mC.rpAlphaBetaListReverse());
            
            i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
            i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
          }
        }
      }
    
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.rpAlphaBetaListReverse() != NULL)
          free(i_node->mC.rpAlphaBetaListReverse());
      }
    
    
      ////////////////// Numerator accumulation //////////////////
      mSearchPaths = SP_TRUE_ONLY;
    
    
      ForwardBackward(rObsMx, nFrames);
    
      //for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) {
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.mType & NT_MODEL && i_node->mC.rpAlphaBetaList() != NULL &&
          i_node->mC.rpAlphaBetaList()->mTime == 0) {
          i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
          i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
        }
      }
    
      p_hmms_alig->ResetXformInstances();
    
      if (p_hmms_alig != p_hmms_upd) 
      {
        p_hmms_upd->ResetXformInstances();
      }
    
      for (i = 0; i < p_hmms_alig->mTotalDelay; i++) 
      {
        p_hmms_alig->UpdateStacks(rObsMx[i],
                    i-p_hmms_alig->mTotalDelay, FORWARD);
      }
    
      if (p_hmms_alig != p_hmms_upd) 
      {
        for (i = 0; i < p_hmms_upd->mTotalDelay; i++) 
        {
          p_hmms_upd->UpdateStacks(rObsMx2[i], i-p_hmms_upd->mTotalDelay, FORWARD);
        }
      }  
    
      // mpMixPCache might be used to cache likelihoods of mixtures of target
      // models. Reallocate the cache to fit mixtures of both models and reset it.
      k = HIGHER_OF(p_hmms_upd->mNMixtures, p_hmms_alig->mNMixtures);
      mpMixPCache = (Cache *) realloc(mpMixPCache, k * sizeof(Cache));
      if (mpMixPCache == NULL) Error("Insufficient memory");
    
      for (i = 0; i < k; i++) mpMixPCache[i].mTime = UNDEF_TIME;
    
    // Update accumulators
      for (mTime = 0; mTime < nFrames; mTime++) 
      {//for every frame
        FLOAT* obs  = rObsMx [mTime+p_hmms_alig->mTotalDelay];
        FLOAT* obs2 = rObsMx2[mTime+p_hmms_upd->mTotalDelay];
        
        p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
        
        if (p_hmms_alig != p_hmms_upd) 
        {
          p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
        }
    
        for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
        { //for every model
          if (i_node->mC.mType & NT_MODEL 
          && i_node->mC.rpAlphaBetaList() != NULL 
          && i_node->mC.rpAlphaBetaList()->mTime == mTime+1) 
          {
    
            AlphaBetaMPE *st;
            int Nq       = i_node->mC.mpHmm->mNStates;
            FLOAT *aq    = i_node->mC.mpHmm->        mpTransition->mpMatrixO;
            FLOAT *aqacc = i_node->mC.rpHmmToUpdate()->mpTransition->mpMatrixO + SQR(Nq);
    //        int qt_1 = (mnumberofnetstates * mtime) + i_node->mC.memittingstateid;
    //        int qt = qt_1 + mnumberofnetstates;
    
            st = i_node->mC.rpAlphaBetaList()->mpState;
    
            if (//!mmi_den_pass &&
              st[Nq-1].mAlpha + st[Nq-1].mBeta - TP > MIN_LOG_WEGIHT) 
            {
              for (i = 0; i < Nq - 1; i++) 
              {
                LOG_INC(aqacc[i * Nq + Nq-1], aq[i * Nq + Nq-1]  * mTranScale +
                   (st[i].mAlpha + st[Nq-1].mBeta - TP) * mOcpScale);
              }
            }
    
            for (j = 1; j < Nq - 1; j++) 
            {                   //for every emitting state
              if (st[j].mAlpha + st[j].mBeta - TP > MIN_LOG_WEGIHT) 
              {
                FLOAT bjto =mpOutPCache[p_hmms_alig->mNStates * mTime +
                                          i_node->mC.mpHmm->mpState[j-1]->mID].mValue;
                // forwardbackward() set mpoutpcache to contain out prob. for all frames
    
                assert(i_node->mC.rpAlphaBetaListReverse()->mTime == mTime);
    
    //            if (!mmi_den_pass) {
                for (i = 0; i < Nq - 1; i++) 
                {
                  LOG_INC(aqacc[i * Nq + j],
                          aq[i * Nq + j]    * mTranScale +
                          (i_node->mC.rpAlphaBetaListReverse()->mpState[i].mAlpha +
                          bjto              * mOutpScale +
                          st[j].mBeta - TP)   * mOcpScale);
                }
    //            }
    
                if (origAccumType == AT_MMI)
                {
                  ReestState(
                      &(*i_node),
                      j-1,
                      (st[j].mAlpha + st[j].mBeta - TP)  * mOcpScale,
                      weight, 
                      obs, 
                      obs2);
                }
                else // origaccumtype == at_mce
                {
                  ReestState(&(*i_node), j-1,
                      (st[j].mAlpha + st[j].mBeta - TP + LogAdd(TP,P) - P)  * mOcpScale,
                      weight, 
                      obs, 
                      obs2);
                }
              }
            }
            
            if (i_node->mC.rpAlphaBetaListReverse()) 
              free(i_node->mC.rpAlphaBetaListReverse());
              
            i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
            i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
          }
        }
      }
    
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.rpAlphaBetaListReverse() != NULL)
          free(i_node->mC.rpAlphaBetaListReverse());
      }
    
      mAccumType = origAccumType;
      return TP;
    } 
  // MCEReest
  //***************************************************************************
  
extern std::map<std::string,FLOAT> state_posteriors;

  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    BaumWelchReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, 
        int nFrames, FLOAT weight, BasicVector<FLOAT>* pWeightVector)
    {
      struct FWBWRet          fwbw;
      FLOAT                   P;
      FLOAT                   update_dir;
      int                     i;
      int                     j;
      int                     k;
      ModelSet*               p_hmms_alig = mpModelSet;
      ModelSet*               p_hmms_upd = mpModelSetToUpdate;
      //typename NetworkType::NodeType*                   node;
      typename NetworkType::iterator                    i_node;
    
      fwbw = ForwardBackward(rObsMx, nFrames);
      P    = fwbw.totLike;
      
      if (P < LOG_MIN) 
        return LOG_0;
    
#ifdef MOTIF
      FLOAT* ocprob = (FLOAT*) malloc(mNumberOfNetStates * (nFrames+1) *
          sizeof(FLOAT));
      for (i=0; i<mNumberOfNetStates * (nFrames+1); i++) ocprob[i] = 0;
#endif
    
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.mType & NT_MODEL && i_node->mC.rpAlphaBetaList() != NULL &&
          i_node->mC.rpAlphaBetaList()->mTime == 0) 
        {
          i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
          i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
        }
      }
    
      p_hmms_alig->ResetXformInstances();
    
      if (p_hmms_alig != p_hmms_upd)
        p_hmms_upd->ResetXformInstances();
    
      for (i = 0; i < p_hmms_alig->mTotalDelay; i++)
      {
        p_hmms_alig->UpdateStacks(rObsMx[i],
                                  i-p_hmms_alig->mTotalDelay, FORWARD);
      }
    
      if (p_hmms_alig != p_hmms_upd)
      {
        for (i = 0; i < p_hmms_upd->mTotalDelay; i++)
        {
          p_hmms_upd->UpdateStacks(rObsMx2[i],
                      i-p_hmms_upd->mTotalDelay, FORWARD);
        }
      }
    
      // mpMixPCache might be used to cache likelihoods of mixtures of target
      // models. Reallocate the cache to fit mixtures of both models and reset it.
      k = HIGHER_OF(p_hmms_upd->mNMixtures, p_hmms_alig->mNMixtures);
      mpMixPCache = (Cache *) realloc(mpMixPCache, k * sizeof(Cache));
      
      if (mpMixPCache == NULL) { 
        Error("Insufficient memory");
      }
      
      for (i = 0; i < k; i++) { 
        mpMixPCache[i].mTime = UNDEF_TIME;
      }
      
      // Update accumulators
      for (mTime = 0; mTime < nFrames; mTime++) 
      { //for every frame
        FLOAT* obs  = rObsMx [mTime + p_hmms_alig->mTotalDelay];
        FLOAT* obs2 = rObsMx2[mTime + p_hmms_upd->mTotalDelay];
        FLOAT  frame_weight = (NULL == pWeightVector) ? 1.0 :
             (*pWeightVector)[mTime];
        
        p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
        
        if (p_hmms_alig != p_hmms_upd) {
          p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
        }
    
        for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
        { //for every model
    //    for (k=0; k < nnodes; k++) {\\}
    //      NetworkType::NodeType *i_node = &mpNodes[k];
          if (i_node->mC.mType & NT_MODEL &&
            i_node->mC.rpAlphaBetaList() != NULL &&
            i_node->mC.rpAlphaBetaList()->mTime == mTime+1) 
          {
            struct AlphaBetaMPE *st;
            int Nq       = i_node->mC.mpHmm->mNStates;
            FLOAT* aq    = i_node->mC.mpHmm->        mpTransition->mpMatrixO;
            FLOAT* aqacc = i_node->mC.rpHmmToUpdate()->mpTransition->mpMatrixO + SQR(Nq);
  //        int qt_1 = (mNumberOfNetStates * mTime) + i_node->mC.mEmittingStateId;
  //        int qt = qt_1 + mNumberOfNetStates;
    
            st = i_node->mC.rpAlphaBetaList()->mpState;
    
            if (//!mmi_den_pass &&
              st[Nq-1].mAlpha + st[Nq-1].mBeta - P > MIN_LOG_WEGIHT) 
            {
              for (i = 0; i < Nq - 1; i++) 
              {
                LOG_INC(aqacc[i * Nq + Nq-1], aq[i * Nq + Nq-1]  * mTranScale +
                    (st[i].mAlpha + st[Nq-1].mBeta - P) * mOcpScale);
              }
            }
    
            for (j = 1; j < Nq - 1; j++) 
            { //for every emitting state
              if (st[j].mAlpha + st[j].mBeta - P > MIN_LOG_WEGIHT) 
              {
#ifdef MOTIF
                ocprob[mNumberOfNetStates * (mTime+1) + 
                  i_node->mC.mEmittingStateId + j]
  //            = i_node->mC.mPhoneAccuracy;
                = my_exp(st[j].mAlpha+st[j].mBeta-P) *
                  ((1-2*static_cast<int>(st[j].mAlphaAccuracy.negative)) * my_exp(st[j].mAlphaAccuracy.logvalue - st[j].mAlpha) +
                  (1-2*static_cast<int>(st[j].mBetaAccuracy.negative))  * my_exp(st[j].mBetaAccuracy.logvalue  - st[j].mBeta)
                  - (1-2*static_cast<int>(fwbw.avgAccuracy.negative)) * my_exp(fwbw.avgAccuracy.logvalue)
                  );
    
#endif
  //            int qt_1   = qt - mNumberOfNetStates;
                FLOAT bjtO =mpOutPCache[p_hmms_alig->mNStates * mTime +
                                          i_node->mC.mpHmm->mpState[j-1]->mID].mValue;
                // ForwardBackward() set mpOutPCache to contain out prob. for all frames
    
                assert(i_node->mC.rpAlphaBetaListReverse()->mTime == mTime);
    
/*printf("AB: time: %d node: %s state: %d alpha: %9g beta: %9g P: %9g occup: %9g\n",
       mTime,
       i_node->mC.mpHmm->mpMacro->mpName, j, 
       i_node->mC.rpAlphaBetaListReverse()->mpState[j].mAlpha,
       i_node->mC.rpAlphaBetaListReverse()->mpState[j].mBeta, P,
       st[j].mAlpha + st[j].mBeta - P);*/

  //            if (!mmi_den_pass) {
                for (i = 0; i < Nq - 1; i++) 
                {
                  LOG_INC(aqacc[i * Nq + j],
                          aq[i * Nq + j] * mTranScale +
                          (i_node->mC.rpAlphaBetaListReverse()->mpState[i].mAlpha +
                          bjtO * mOutpScale + st[j].mBeta - P) * mOcpScale);
                }
  //            }
    
                if (mAccumType == AT_MFE || mAccumType == AT_MPE) 
                {
                  update_dir = (1-2*static_cast<int>(st[j].mAlphaAccuracy.negative)) * my_exp(st[j].mAlphaAccuracy.logvalue - st[j].mAlpha) +
                               (1-2*static_cast<int>(st[j].mBetaAccuracy.negative))  * my_exp(st[j].mBetaAccuracy.logvalue  - st[j].mBeta)  -
                               (1-2*static_cast<int>(fwbw.avgAccuracy.negative))     * my_exp(fwbw.avgAccuracy.logvalue);
                } 
                else 
                {
                  update_dir = 1.0;
                }
//string state_name = i_node->mC.mpHmm->mpState[j - 1]->mpMacro->mpName;
//state_posteriors[state_name.substr(0, state_name.find('_'))] += my_exp(st[j].mAlpha + st[j].mBeta - P);

                ReestState(
                  &(*i_node), 
                  j - 1,
                  (st[j].mAlpha + st[j].mBeta - P)  * mOcpScale,
                  update_dir*weight*frame_weight, 
                  obs, 
                  obs2);
              }
            }
            
             if (i_node->mC.rpAlphaBetaListReverse()) 
              free(i_node->mC.rpAlphaBetaListReverse());
              
            i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
            i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
          }
        }
//for(std::map<string,FLOAT>::iterator i_sp = state_posteriors.begin(); i_sp != state_posteriors.end(); ++i_sp) {
//  printf("%d %s %f\n", mTime, i_sp->first.c_str(), i_sp->second);
//}
//state_posteriors.clear();
      }
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.rpAlphaBetaListReverse() != NULL)
          free(i_node->mC.rpAlphaBetaListReverse());
      }
    
#ifdef MOTIF
      FLOAT max = LOG_0;
      printf("mTranScale: %f\noutpScale: %f\n",mTranScale, mOutpScale);
      for (i = 0; i < mNumberOfNetStates * (nFrames+1); i++) 
        max = HIGHER_OF(max, ocprob[i]);
    
      imagesc(ocprob, mNumberOfNetStates, (nFrames+1),
              sizeof(FLOAT) == 4 ? "float" : "double",
              NULL,  cm_color, "OccProb");
    
      for (j=0; j<nFrames+1; j++) 
      {
        for (i=0; i<mNumberOfNetStates; i++) 
        {
          printf("%6.3g ", (ocprob[mNumberOfNetStates * j + i]));
        }
        printf("\n");
      }
  //  for (i = 0; i < mNumberOfNetStates * (nFrames+1); i++) {
  //    if (ocprob[i] < LOG_MIN) ocprob[i] = max;
  //  }
  //  imagesc(ocprob, mNumberOfNetStates, (nFrames+1), STR(FLOAT), NULL, cm_color, "Log OccProb");
      free(ocprob);
#endif
    
      return P;
    }
    //***************************************************************************




  //***************************************************************************
  //***************************************************************************
  template<typename _NetworkType>
    FLOAT 
    Decoder<_NetworkType>::
    GetMpeGamma(const Matrix<FLOAT>& rObsMx, Matrix<FLOAT>& rGamma, FLOAT& avgAcc, 
        int nFrames, FLOAT weight, BasicVector<FLOAT>* pWeightVector)
    {
      struct FWBWRet          fwbw;
      FLOAT                   P;
      FLOAT                   update_dir;
      int                     i;
      int                     j;
      int                     k;
      ModelSet*               p_hmms_alig = mpModelSet;
      ModelSet*               p_hmms_upd = mpModelSetToUpdate;
      //typename NetworkType::NodeType*                   node;
      typename NetworkType::iterator                    i_node;
    
      fwbw = ForwardBackward(rObsMx, nFrames);
      P    = fwbw.totLike;
      avgAcc = (1-2*static_cast<int>(fwbw.avgAccuracy.negative))
                 * my_exp(fwbw.avgAccuracy.logvalue);

      rGamma.Init(nFrames,mpModelSet->mNStates);
      
      if (P < LOG_MIN) 
        return LOG_0;
    
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.mType & NT_MODEL && i_node->mC.rpAlphaBetaList() != NULL &&
          i_node->mC.rpAlphaBetaList()->mTime == 0) 
        {
          //store top of the stack of AlphaBeta
          i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
          //pop the stack AlphaBeta
          i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
        }
      }
    
           
      // Compute gammas:
      // We assume the AlphaBeta stack has earliest time on top,
      // for each frame we pop a AlphaBeta record of each model,
      // if the model was active at the time
      for (mTime = 0; mTime < nFrames; mTime++) 
      { //for every frame
        FLOAT  frame_weight = (NULL == pWeightVector) ? 1.0 :
             (*pWeightVector)[mTime];
        
        for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
        { //for every model
          if (i_node->mC.mType & NT_MODEL &&
            i_node->mC.rpAlphaBetaList() != NULL &&
            i_node->mC.rpAlphaBetaList()->mTime == mTime+1) 
          {
            struct AlphaBetaMPE *st;
            int Nq       = i_node->mC.mpHmm->mNStates;
            st = i_node->mC.rpAlphaBetaList()->mpState;
    
            for (j = 1; j < Nq - 1; j++) 
            { //for every emitting state
              if (st[j].mAlpha + st[j].mBeta - P > MIN_LOG_WEGIHT) 
              {
#ifdef MOTIF
                ocprob[mNumberOfNetStates * (mTime+1) + 
                  i_node->mC.mEmittingStateId + j]
  //            = i_node->mC.mPhoneAccuracy;
                = my_exp(st[j].mAlpha+st[j].mBeta-P) *
                  ((1-2*static_cast<int>(st[j].mAlphaAccuracy.negative)) * my_exp(st[j].mAlphaAccuracy.logvalue - st[j].mAlpha) +
                  (1-2*static_cast<int>(st[j].mBetaAccuracy.negative))  * my_exp(st[j].mBetaAccuracy.logvalue  - st[j].mBeta)
                  - (1-2*static_cast<int>(fwbw.avgAccuracy.negative)) * my_exp(fwbw.avgAccuracy.logvalue)
                  );
    
#endif
                    
                assert(i_node->mC.rpAlphaBetaListReverse()->mTime == mTime);
    
                update_dir = my_exp(st[j].mAlpha + st[j].mBeta - P); //update_dir = gama_q
    
                if (mAccumType == AT_MFE || mAccumType == AT_MPE) 
                {
                  //update_dir = gama_q^MPE = gama_q * (aplha_q'+beata_q' -c_avg)
                  update_dir *= (1-2*static_cast<int>(st[j].mAlphaAccuracy.negative)) * my_exp(st[j].mAlphaAccuracy.logvalue - st[j].mAlpha) +
                                (1-2*static_cast<int>(st[j].mBetaAccuracy.negative))  * my_exp(st[j].mBetaAccuracy.logvalue  - st[j].mBeta)  -
                                (1-2*static_cast<int>(fwbw.avgAccuracy.negative))     * my_exp(fwbw.avgAccuracy.logvalue);
                } 
                //store the gamma
                size_t PDF_obs_coef = i_node->mC.mpHmm->mpState[j-1]->PDF_obs_coef;
                rGamma(mTime, PDF_obs_coef) += weight*frame_weight*update_dir;
              }
            }
            
             if (i_node->mC.rpAlphaBetaListReverse()) 
              free(i_node->mC.rpAlphaBetaListReverse());
              
            i_node->mC.rpAlphaBetaListReverse() = i_node->mC.rpAlphaBetaList();
            i_node->mC.rpAlphaBetaList() = i_node->mC.rpAlphaBetaList()->mpNext;
          }
        }
      }
      for (i_node = rNetwork().begin(); i_node != rNetwork().end(); ++i_node) 
      {
        if (i_node->mC.rpAlphaBetaListReverse() != NULL)
          free(i_node->mC.rpAlphaBetaListReverse());
      }
    
#ifdef MOTIF
      FLOAT max = LOG_0;
      printf("mTranScale: %f\noutpScale: %f\n",mTranScale, mOutpScale);
      for (i = 0; i < mNumberOfNetStates * (nFrames+1); i++) 
        max = HIGHER_OF(max, ocprob[i]);
    
      imagesc(ocprob, mNumberOfNetStates, (nFrames+1),
              sizeof(FLOAT) == 4 ? "float" : "double",
              NULL,  cm_color, "OccProb");
    
      for (j=0; j<nFrames+1; j++) 
      {
        for (i=0; i<mNumberOfNetStates; i++) 
        {
          printf("%6.3g ", (ocprob[mNumberOfNetStates * j + i]));
        }
        printf("\n");
      }
  //  for (i = 0; i < mNumberOfNetStates * (nFrames+1); i++) {
  //    if (ocprob[i] < LOG_MIN) ocprob[i] = max;
  //  }
  //  imagesc(ocprob, mNumberOfNetStates, (nFrames+1), STR(FLOAT), NULL, cm_color, "Log OccProb");
      free(ocprob);
#endif
    
      return P;
    }
    //***************************************************************************




  //***************************************************************************
  //***************************************************************************
  template<typename _NodeType>
    void 
    FreeWordLinkRecords(WordLinkRecord<_NodeType>* wlr)
    {
      typedef typename WordLinkRecord<_NodeType>::AltHypList AltHypList;

      if (wlr != NULL) 
      {
        --wlr->mNReferences;
        assert(wlr->mNReferences >= 0);
    
        if (wlr->mNReferences == 0) 
        {
          FreeWordLinkRecords(wlr->mpNext);
          
          if (wlr->mpAltHyps)
          {
            for(typename AltHypList::iterator i = wlr->mpAltHyps->begin(); i != wlr->mpAltHyps->end(); i++)
              FreeWordLinkRecords(i->mpWlr);
          
            delete wlr->mpAltHyps;
          }
          
#ifdef DEBUG_MSGS
          assert(wlr->mIsFreed == false);
          wlr->mIsFreed = true;
#else
          //free(wlr);
          delete wlr;
#endif
        }
      }
    }
    //***************************************************************************

  
  //***************************************************************************
  //***************************************************************************
  template<typename _NodeType>
    void 
    Decoder<_NodeType>::
    Wlr2Lattice(WordLinkRecord* pWlr, Lattice& rLattice)
    {
      
      Lattice::NodeType* p_first = NULL;
      
      Wlr2Lattice_AllocateNodes(pWlr, p_first);
      Wlr2Lattice_EstablishLinks(pWlr);
      
      Lattice::NodeType* p_node;
      
      for (p_node = p_first; NULL != p_node && NULL != p_node->mpNext;
           p_node = p_node->mpNext)
      { }

      rLattice.splice(p_first, p_node);
      
    }
  // Wlr2Lattice(WordLinkRecord* pWlr, Lattice& rLattice)
  //***************************************************************************

  //***************************************************************************
  //***************************************************************************
  template<typename _NodeType>
    void 
    Decoder<_NodeType>::
    Wlr2Lattice_AllocateNodes(WordLinkRecord* pWlr, 
        typename STK::Lattice::NodeType*& rpNode)
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
          
      switch (pWlr->pNode()->mC.mType & (NT_WORD | NT_MODEL | NT_PHONE)) 
      {
        case NT_WORD:  pNode->mC.mpPronun = pWlr->pNode()->mC.mpPronun; break;
        case NT_MODEL: pNode->mC.mpHmm    = pWlr->pNode()->mC.mpHmm;    break;
        case NT_PHONE: pNode->mC.mpName   = pWlr->pNode()->mC.mpName;   break;
        default:       Error("Fatal: Invalid node type");
      }
        
      pNode->mC.mType       = pWlr->pNode()->mC.mType;
      pNode->mC.SetStart(UNDEF_TIME);
      pNode->mC.SetStop(pWlr->mTime);
      pNode->mC.mpAlphaBeta = NULL;
        
      pNode->rNBackLinks() = pWlr->mpNext    == NULL ? 0 :
                           pWlr->mpAltHyps == NULL ? 1 : pWlr->mpAltHyps->size() + 1;
                           
      // If rpNode == NULL, pWlr is the wery last WLR referenced by token in last state
      pNode->rNLinks()     = rpNode == NULL ? 0 : pWlr->mNReferences;


      // Allocate space for links and backlinks      
      pNode->rpLinks()     = (Lattice::LinkType *) malloc(pNode->rNLinks()     * sizeof(Lattice::LinkType));
      pNode->rpBackLinks() = (Lattice::LinkType *) malloc(pNode->rNBackLinks() * sizeof(Lattice::LinkType));
      if (pNode->rpLinks() == NULL || pNode->rpBackLinks() == NULL) 
        Error("Insufficient memory");
          

      for (int ii = 0; ii < pNode->NLinks(); ++ii)
      {
        pNode->rpLinks()[ii].mAcousticLike = LOG_0;
        pNode->rpLinks()[ii].mLmLike = LOG_0;
      }

      for (int ii = 0; ii < pNode->NBackLinks(); ++ii)
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
        
        
      // Set pWlr->mpNode to the new lattice node, so we can latter easily
      // establish  links between the nodes in the lattice nodes
      pWlr->mpNode = reinterpret_cast<typename WordLinkRecord::NodeType*>(pNode);
        
      // All successors have been already processed, so continue recursively with Wlr's predecessors
      if(pWlr->mpNext != NULL)
      {
        Wlr2Lattice_AllocateNodes(pWlr->mpNext, rpNode);
          
        if(pWlr->mpAltHyps != NULL)
        {
          for(typename WordLinkRecord::AltHypList::iterator i = pWlr->mpAltHyps->begin(); i != pWlr->mpAltHyps->end(); i++)
            Wlr2Lattice_AllocateNodes(i->mpWlr, rpNode);
        }
      }
    }
  // MakeLatticeNodesForWordLinkRecords(WordLinkRecord* pWlr, NodeType*& rpNode)
  //***************************************************************************

  

  //***************************************************************************
  //***************************************************************************
  template<typename _NodeType>
    void 
    Decoder<_NodeType>::
    Wlr2Lattice_EstablishLinks(WordLinkRecord* pWlr)
    {
      typedef typename WordLinkRecord::AltHypList AltHypList;

      // After calling MakeLatticeNodesForWordLinkRecords, mAux was set to 
      // pWlr->mNReferencesen.
      // Use it  to count down the successors that has been already processed. 
      if(--pWlr->mAux > 0) 
        return;
        
      // All successors have been already processed, so continue recursively with 
      // Wlr's predecessors
      if(pWlr->mpNext != NULL)
      {
        int j = 0;
        WordLinkRecord* p = pWlr->mpNext;
        typename WordLinkRecord::LikeType aux_total_like = 
          pWlr->mLike - p->mLike;
        typename WordLinkRecord::LikeType aux_acoustic_like = 
          pWlr->mAcousticLike - p->mAcousticLike;
        typename WordLinkRecord::LikeType aux_lm_like = 
          float_safe_substract(aux_total_like, aux_acoustic_like, 3);

        if (pWlr->pNode()->mC.mType & NT_MODEL) {
          aux_lm_like = (aux_lm_like - mMPenalty) /mLmScale;
        }
        else if (pWlr->pNode()->mC.mType & NT_WORD
        &&       pWlr->pNode()->mC.mpPronun != NULL) {
          aux_lm_like = (aux_lm_like - mWPenalty) / mLmScale;
        }
        else {
          aux_lm_like = aux_lm_like / mLmScale;
        }

        p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetNode(pWlr->pNode());
        p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetLmLike(aux_lm_like);
        p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetAcousticLike(aux_acoustic_like);

        pWlr->pNode()->rpBackLinks()[j].SetNode(p->pNode());
        pWlr->pNode()->rpBackLinks()[j].SetLmLike(aux_lm_like);
        pWlr->pNode()->rpBackLinks()[j].SetAcousticLike(aux_acoustic_like);

        //printf("1 i->mLike         = %f   p->mLike         = %f   diff = %e\n",  pWlr->mLike, p->mLike, static_cast<float>(pWlr->mLike - p->mLike)); 
        //printf("1 i->mAcousticLike = %f   p->mAcousticLike = %f   diff = %e\n",  pWlr->mAcousticLike, p->mAcousticLike, static_cast<float>(pWlr->mAcousticLike - p->mAcousticLike)); 
        //printf("1 i->mLike - i->mAcousticLike = %e\n",  static_cast<float>(pWlr->mLike - pWlr->mAcousticLike)); 
        //printf("1 p->mLike - p->mAcousticLike = %e\n",  static_cast<float>(p->mLike - p->mAcousticLike)); 

        // recursively do the rest of the WLR's
        Wlr2Lattice_EstablishLinks(pWlr->mpNext);
          
        if(pWlr->mpAltHyps != NULL)
        {
          typename AltHypList::iterator i;
          
          for(j = 1, i = pWlr->mpAltHyps->begin(); i != pWlr->mpAltHyps->end(); j++, i++)
          {
            WordLinkRecord* p = i->mpWlr;
            typename WordLinkRecord::LikeType aux_total_like  = i->mLike - p->mLike;
            typename WordLinkRecord::LikeType aux_acoustic_like = i->mAcousticLike - p->mAcousticLike;
            typename WordLinkRecord::LikeType aux_lm_like = float_safe_substract(aux_total_like, aux_acoustic_like, 3);

            if (pWlr->pNode()->mC.mType & NT_MODEL) {
              aux_lm_like = (aux_lm_like - mMPenalty) / mLmScale;
            }
            else if (pWlr->pNode()->mC.mType & NT_WORD
            &&       pWlr->pNode()->mC.mpPronun != NULL) {
              aux_lm_like = (aux_lm_like - mWPenalty) / mLmScale;
            }
            else {
              aux_lm_like = aux_lm_like / mLmScale;
            }
            
            p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetNode(pWlr->pNode());
            p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetLmLike(aux_lm_like);
            p->pNode()->rpLinks()[p->mNReferences - p->mAux].SetAcousticLike(aux_acoustic_like);

            pWlr->pNode()->rpBackLinks()[j].SetNode(p->pNode());
            pWlr->pNode()->rpBackLinks()[j].SetLmLike(aux_lm_like);
            pWlr->pNode()->rpBackLinks()[j].SetAcousticLike(aux_acoustic_like);

            //printf("n i->mLike         = %f   p->mLike         = %f   diff = %e\n",  i->mLike, p->mLike, static_cast<float>(i->mLike - p->mLike)); 
            //printf("n i->mAcousticLike = %f   p->mAcousticLike = %f   diff = %e\n",  i->mAcousticLike, p->mAcousticLike, static_cast<float>(i->mAcousticLike - p->mAcousticLike)); 
            //printf("n i->mLike - i->mAcousticLike = %e\n",  static_cast<float>(i->mLike - i->mAcousticLike)); 
            //printf("n p->mLike - p->mAcousticLike = %e\n",  static_cast<float>(p->mLike - p->mAcousticLike)); 
          
            // recursively do the rest of the WLR's
            Wlr2Lattice_EstablishLinks(i->mpWlr);
          }
        }
      }
    }
  // EstablishLinks(WordLinkRecord* pWlr)
  //***************************************************************************
  

#if ENABLE_SSE
# define OPTIMIZE_GAUSSIAN_COMPUTATION      
#endif

  //***************************************************************************
  //***************************************************************************
  // The function does the gaussian dist exponent operations.
  // It is declared inline as it is called from one place only
  inline FLOAT 
  compute_diag_c_gaussian_density(
    const FLOAT*  pObs, 
    const FLOAT   gConst,
    const FLOAT*  pMean,
    const FLOAT*  pVar,
    const size_t  vSize)
  {
#ifdef __GNUC__
    FLOAT m_like __attribute__ ((aligned (16))) = 0.0;
#else
    FLOAT m_like = 0.0;
#endif
    size_t j;

#if defined(OPTIMIZE_GAUSSIAN_COMPUTATION) && defined(__GNUC__ )
// THE PIECE OF CODE IN THIS "IF" BRANCH IS THE SSE OPTIMIZED
// WAY OF COMPUTING MULTI GAUSSIAN DENSITY. IT ASSUMES THAT 
// ALL POINTERS ARE 16-BYTES ALIGNED (SSE REQUIREMENT)

#  if !DOUBLEPRECISION
    // this is a stric optimization
    f4vector        l = {{0.0F}}; 
    const f4vector* o = reinterpret_cast<const f4vector*>(pObs); 
    const f4vector* m = reinterpret_cast<const f4vector*>(pMean);
    const f4vector* v = reinterpret_cast<const f4vector*>(pVar);
    
    l.f[0] = l.f[1] = l.f[2] = l.f[3] = 0.0F;

    for (j = 0; j < vSize; j += 4) 
    {
      l.v += SQR(o->v - m->v) * v->v;
      // v4sf tmp;
      // tmp = __builtin_ia32_subps(o->v, m->v);
      // tmp = __builtin_ia32_mulps(tmp,  tmp);
      // tmp = __builtin_ia32_mulps(tmp,  v->v);
      // l.v = __builtin_ia32_addps(tmp,  l.v);
      o++;
      m++;
      v++;
    }
    
    m_like = l.f[0] + l.f[1] + l.f[2] + l.f[3];
#  else
    // this is a stric optimization
    d2vector        l = {{0.0F}}; 
    const d2vector* o = reinterpret_cast<const d2vector*>(pObs); 
    const d2vector* m = reinterpret_cast<const d2vector*>(pMean);
    const d2vector* v = reinterpret_cast<const d2vector*>(pVar);
    
    l.f[0] = l.f[1] = 0.0;

    for (j = 0; j < vSize; j += 2) 
    {
      l.v += SQR(o->v - m->v) * v->v;
      o++;
      m++;
      v++;
    }
    
    m_like = l.f[0] + l.f[1];
#  endif // DOUBLEPRECISION
    
#else
// NO OPTIMIZATION IS PERFORMED IN THIS PIECE OF CODE
    // the original loop
    for (j = 0; j < vSize; j++) 
    {                   
      m_like += SQR(pObs[j] - pMean[j]) * pVar[j];
    }
#endif
    
    m_like = -0.5 * (gConst + m_like);
    return m_like;
  } // compute_diag_c_gaussian_density(...)
  //***************************************************************************

}
