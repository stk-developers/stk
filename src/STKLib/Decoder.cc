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


//#############################################################################
//#############################################################################
// PROJECT INCLUDES
//#############################################################################
//#############################################################################
#include "Decoder.h"
#include "Lattice.h"
#include "labels.h"
#include "common.h"


//#############################################################################
//#############################################################################
// SYSTEM INCLUDES
//#############################################################################
//#############################################################################
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <iostream>
#include <sstream>

#include <boost/tokenizer.hpp>

#ifdef MOTIF
#include "imagesc.h"
#endif

//#define TRACE_TOKENS
#define SQR(x) ((x) * (x))


//#############################################################################
//#############################################################################
// CODE
//#############################################################################
//#############################################################################

namespace STK
{

  #ifdef DEBUG_MSGS
  WordLinkRecord * firstWLR;
  #endif
  
  
  //***************************************************************************
  //***************************************************************************
  WordLinkRecord*
  Decoder::
  TimePruning(int frame_delay)
  {
    size_t              i;
    NetworkType::NodeType*  p_node;
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
    for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mpAnr->mpNextActiveModel)
    {
      assert((p_node->mType & NT_MODEL) && p_node->mpAnr->mIsActiveModel);
  
      for (i = 0; i < p_node->mpHmm->mNStates-1; i++) 
      {
        if (p_node->mpAnr->mpTokens[i].IsActive()) 
        {
          Token* token = &p_node->mpAnr->mpTokens[i];
  
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
  void 
  Decoder::
  WriteAlpha(int time, NetworkType::NodeType* pNode, int state, Token *token)
  {
    if (pNode->mpAlphaBetaListReverse == NULL ||
      pNode->mpAlphaBetaListReverse->mTime != time) 
    {
      size_t    i;
      FWBWR *   newrec;
      
      newrec  = (FWBWR*) malloc(sizeof(FWBWR) +
                                sizeof(newrec->mpState[0]) * (pNode->mpHmm->mNStates-1));
      if (newrec == NULL) 
        Error("Insufficient memory");
        
      newrec->mpNext = pNode->mpAlphaBetaListReverse;
      newrec->mTime = time;
      
      for (i=0; i<pNode->mpHmm->mNStates; i++) 
      {
        newrec->mpState[i].mAlpha = newrec->mpState[i].mBeta = LOG_0;
      }
      
      pNode->mpAlphaBetaListReverse = newrec;
    }
    pNode->mpAlphaBetaListReverse->mpState[state].mAlpha = token->mLike;
    pNode->mpAlphaBetaListReverse->mpState[state].mAlphaAccuracy = token->mAccuracy;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  WriteBeta(int time, NetworkType::NodeType* pNode, int state, Token *token)
  {  
    // Record for current time must be already moved to
    // mpAlphaBetaList by function BackwardPruning
    assert(pNode->mpAlphaBetaListReverse == NULL ||
          pNode->mpAlphaBetaListReverse->mTime < time);
  
    if (pNode->mpAlphaBetaList != NULL && pNode->mpAlphaBetaList->mTime == time) {
      pNode->mpAlphaBetaList->mpState[state].mBeta = token->mLike;
      pNode->mpAlphaBetaList->mpState[state].mBetaAccuracy = token->mAccuracy;
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  Decoder::
  BackwardPruning(int time, NetworkType::NodeType* pNode, int state)
  {
    while (pNode->mpAlphaBetaListReverse != NULL &&
          pNode->mpAlphaBetaListReverse->mTime > time) 
    {
      FWBWR *fwbwr = pNode->mpAlphaBetaListReverse;
      pNode->mpAlphaBetaListReverse = fwbwr->mpNext;
      free(fwbwr);
    }
  
    if (pNode->mpAlphaBetaListReverse != NULL &&
      pNode->mpAlphaBetaListReverse->mTime == time) 
    {
      FWBWR* fwbwr = pNode->mpAlphaBetaListReverse;
      pNode->mpAlphaBetaListReverse = fwbwr->mpNext;
      fwbwr->mpNext = pNode->mpAlphaBetaList;
      pNode->mpAlphaBetaList = fwbwr;
    }
  
    return !(pNode->mpAlphaBetaList != NULL &&
            pNode->mpAlphaBetaList->mTime == time &&
            pNode->mpAlphaBetaList->mpState[state].mAlpha > LOG_MIN);
  }
  
  
  
  
  #ifndef NDEBUG
  int test_for_cycle = 0;
  int HasCycleCounter = 100000;
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  Decoder::
  HasCycle() 
  {
    NetworkType::NodeType* node;
    
    HasCycleCounter++;
    
    if (!test_for_cycle) 
      return 0;
      
    for (node = mpActiveNodes; node != NULL; node = node->mpAnr->mpNextActiveNode) 
    {
      int     i;
      int     n_links = InForwardPass() ? node->mNLinks : node->mNBackLinks;
      NetworkType::LinkType*   links   = InForwardPass() ? node->mpLinks : node->mpBackLinks;
      
      if (node->mAux2 == HasCycleCounter) 
      {
        printf("Cycle in list of active nodes\n");
        return 1;
      }
      
      node->mAux2 = HasCycleCounter;
  
      for (i = 0; i < n_links; i++)
      {
        if (links[i].pNode()->mAux2 == HasCycleCounter &&
            (!(links[i].pNode()->mType & NT_MODEL)
              || links[i].pNode()->mType & NT_TEE)) 
        {
          printf("Active node %d listed after his non-model succesor %d\n",
                (int) (node - rNetwork().pFirst()), (int) (links[i].pNode() - rNetwork().pFirst()));
  
          return 2;
        }
      }
    }
    return 0;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  bool
  Decoder:: 
  AllWordSuccessorsAreActive() 
  {
    NetworkType::NodeType* node;
    
    if (!test_for_cycle) 
      return true;
  
    for (node = mpActiveNodes; node != NULL; node = node->mpAnr->mpNextActiveNode) 
    {
      int     i;
      int     n_links = InForwardPass() ? node->mNLinks : node->mNBackLinks;
      NetworkType::LinkType *  links   = InForwardPass() ? node->mpLinks : node->mpBackLinks;
  
      for (i=0; i <n_links; i++)
      {
        if (links[i].pNode()->mAux2 != HasCycleCounter &&
          links[i].pNode() != (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()) &&
          (!(links[i].pNode()->mType & NT_MODEL)
          || links[i].pNode()->mType & NT_TEE)) 
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
  void 
  Decoder::
  MarkWordNodesLeadingFrom(NetworkType::NodeType* node)
  {
    int       i;
    int       n_links = InForwardPass() ? node->mNLinks : node->mNBackLinks;
    NetworkType::LinkType*     p_links = InForwardPass() ? node->mpLinks : node->mpBackLinks;
  
    for (i = 0; i < n_links; i++) 
    {
      NetworkType::NodeType* p_lnode = p_links[i].pNode();
      
      if ((p_lnode->mType & NT_MODEL && !(p_lnode->mType & NT_TEE))
      ||  (p_lnode == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()))) 
      {
        continue;
      }
      
      if (p_lnode->mpAnr == NULL) 
        p_lnode->mpAnr = new ActiveNodeRecord(p_lnode);

        
      if (p_lnode->mpAnr->mActiveNodeFlag > 0) 
        continue;
  
      if (p_lnode->mpAnr->mIsActiveModel) 
      {
        assert(p_lnode->mType & NT_TEE);
        continue;
      }
      
      // when this t mActiveNodeFlag is 
      if (p_lnode->mpAnr->mActiveNodeFlag-- == 0) 
      {
        p_lnode->mpAnr->mAux = 0;
        MarkWordNodesLeadingFrom(p_lnode);
      }
    }
  }
  //***************************************************************************
  
  
  
  //***************************************************************************
  //***************************************************************************
  Decoder::NetworkType::NodeType *
  Decoder::
  pActivateWordNodesLeadingFrom(NetworkType::NodeType* pNode)
  {
    int     i;
    int     n_links = InForwardPass() ? pNode->mNLinks : pNode->mNBackLinks;
    NetworkType::LinkType*   p_links   = InForwardPass() ? pNode->mpLinks : pNode->mpBackLinks;
  
    for (i = 0; i < n_links; i++) 
    {
      NetworkType::NodeType* p_lnode = p_links[i].pNode();
      
      if (((p_lnode->mType & NT_MODEL) && !(p_lnode->mType & NT_TEE))
      || (p_lnode == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()))) 
      {
        continue;
      }
      
      assert(p_lnode->mpAnr != NULL) ;
      
      if (p_lnode->mpAnr->mActiveNodeFlag++ > 0) 
        continue;
  
      if (p_lnode->mpAnr->mIsActiveModel) 
      {
        assert(p_lnode->mType & NT_TEE);
        continue;
      }
  
      p_lnode->mpAnr->mAux++;
      
      if (p_lnode->mpAnr->mActiveNodeFlag < 0) 
        continue;
  
      assert(p_lnode->mpAnr->mActiveNodeFlag == 0);
      
      p_lnode->mpAnr->mActiveNodeFlag    = p_lnode->mpAnr->mAux;
      p_lnode->mpAnr->mpNextActiveNode = pNode->mpAnr->mpNextActiveNode;
      p_lnode->mpAnr->mpPrevActiveNode = pNode;
      
      if (pNode->mpAnr->mpNextActiveNode) 
        pNode->mpAnr->mpNextActiveNode->mpAnr->mpPrevActiveNode = p_lnode;
        
      pNode->mpAnr->mpNextActiveNode   = p_lnode;
      pNode = pActivateWordNodesLeadingFrom(p_lnode);
    }
    
    assert(!HasCycle());
    return pNode;
  }
  // pActivateWordNodesLeadingFrom(NetworkType::NodeType* pNode)
  //***************************************************************************
  
  
  
  //***************************************************************************
  //***************************************************************************
  void
  Decoder:: 
  ActivateModel(NetworkType::NodeType* pNode)
  {
    if(pNode->mpAnr == NULL)
    {
      pNode->mpAnr = new ActiveNodeRecord(pNode);
    }
    else
    {
      assert(pNode->mpAnr->mIsActiveModel || pNode->mpAnr->mActiveNodeFlag);
    }
    
    if (pNode->mpAnr->mIsActiveModel) return;
    
    pNode->mpAnr->mIsActiveModel    = true;
    pNode->mpAnr->mpPrevActiveModel = NULL;
    pNode->mpAnr->mpNextActiveModel = mpActiveModels;
    
    if (mpActiveModels != NULL) 
    {
      mpActiveModels->mpAnr->mpPrevActiveModel = pNode;
    }
    
    mpActiveModels = pNode;
  
    if (pNode->mpAnr->mActiveNodeFlag) {
      assert(pNode->mType & NT_TEE);
      return;
    }
  
    // probably not necessary; when removed assert on line 555 shoud be allowed
    pNode->mpAnr->mActiveNodeFlag = 1; 

    pNode->mpAnr->mpPrevActiveNode = NULL;
    pNode->mpAnr->mpNextActiveNode = mpActiveNodes;
    
    if (mpActiveNodes != NULL) 
    {
      mpActiveNodes->mpAnr->mpPrevActiveNode = pNode;
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
  void
  Decoder:: 
  DeactivateWordNodesLeadingFrom(NetworkType::NodeType* pNode)
  {
    int       i;
    int       n_links = InForwardPass() ? pNode->mNLinks : pNode->mNBackLinks;
    NetworkType::LinkType *    links   = InForwardPass() ? pNode->mpLinks : pNode->mpBackLinks;
    
    assert(pNode->mpAnr != NULL);
  
    for (i = 0; i < n_links; i++) 
    {
      NetworkType::NodeType* p_lnode = links[i].pNode();
      
      if ((p_lnode->mType & NT_MODEL && !(p_lnode->mType & NT_TEE))
      ||  (p_lnode->mType & NT_STICKY)
      ||  (p_lnode == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()))) 
      {
        continue;
      }
            
      assert(/*!(p_lnode->mType & NT_TEE) ||*/ p_lnode->mpAnr != NULL && p_lnode->mpAnr->mActiveNodeFlag > 0);
      
      if (--p_lnode->mpAnr->mActiveNodeFlag) 
        continue;
  
      if (p_lnode->mType & NT_TEE && p_lnode->mpAnr->mIsActiveModel) 
        return;
  
      DeactivateWordNodesLeadingFrom(p_lnode);
      assert(p_lnode->mpAnr->mpPrevActiveNode);
      
      p_lnode->mpAnr->mpPrevActiveNode->mpAnr->mpNextActiveNode = p_lnode->mpAnr->mpNextActiveNode;

      if (p_lnode->mpAnr->mpNextActiveNode)
        p_lnode->mpAnr->mpNextActiveNode->mpAnr->mpPrevActiveNode = p_lnode->mpAnr->mpPrevActiveNode;
    
      delete p_lnode->mpAnr;
      p_lnode->mpAnr = NULL;
  
    }
    
    assert(!HasCycle());
    assert(AllWordSuccessorsAreActive());
  }
  // DeactivateWordNodesLeadingFrom(NetworkType::NodeType* pNode)
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  void
  Decoder:: 
  DeactivateModel(NetworkType::NodeType* pNode)
  {
    assert(pNode->mpAnr != NULL && pNode->mpAnr->mIsActiveModel);
      
    pNode->mpAnr->mIsActiveModel = false;
  
    if (pNode->mpAnr->mpNextActiveModel != NULL) 
    {
      pNode->mpAnr->mpNextActiveModel->mpAnr->mpPrevActiveModel = pNode->mpAnr->mpPrevActiveModel;
    }
  
    if (pNode->mpAnr->mpPrevActiveModel != NULL) 
    {
      pNode->mpAnr->mpPrevActiveModel->mpAnr->mpNextActiveModel = pNode->mpAnr->mpNextActiveModel;
    } else {
      assert(mpActiveModels == pNode);
      mpActiveModels = pNode->mpAnr->mpNextActiveModel;
    }
  
    assert(!HasCycle());
    
    if (pNode->mType & NT_TEE && pNode->mpAnr->mActiveNodeFlag)
      return;
      
//    assert(pNode->mpAnr->mActiveNodeFlag == 0);
    
    DeactivateWordNodesLeadingFrom(pNode);
  
    if (pNode->mpAnr->mpNextActiveNode != NULL) {
      pNode->mpAnr->mpNextActiveNode->mpAnr->mpPrevActiveNode = pNode->mpAnr->mpPrevActiveNode;
    }
  
    if (pNode->mpAnr->mpPrevActiveNode != NULL) {
      pNode->mpAnr->mpPrevActiveNode->mpAnr->mpNextActiveNode = pNode->mpAnr->mpNextActiveNode;
    } else {
      assert(mpActiveNodes == pNode);
      mpActiveNodes = pNode->mpAnr->mpNextActiveNode;
    }
    
    delete pNode->mpAnr;
    pNode->mpAnr = NULL;
  
    assert(!HasCycle());
  }
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  void
  Token::AddAlternativeHypothesis(WordLinkRecord* pWlr, Token::LikeType like, 
      Token::LikeType acousticLike)
  { 
    WlrReference  tmp_wlr_ref(pWlr, like, acousticLike);

    if (mpAltHyps == NULL)
      mpAltHyps = new AltHypList;
    
    assert(pWlr != NULL);

    mpAltHyps->push_back(tmp_wlr_ref);
    pWlr->mNReferences++;
  }
  // Token::AddAlternativeHypothesis(...)
  //***************************************************************************

  
  //***************************************************************************
  //***************************************************************************
  void
  Token::AddAlternativeHypothesis(WordLinkRecord* pWlr)
  { 
    WlrReference  tmp_wlr_ref(pWlr, 0.0, 0.0);

    if (mpAltHyps == NULL)
      mpAltHyps = new AltHypList;
    
    assert(pWlr != NULL);

    mpAltHyps->push_back(tmp_wlr_ref);
    pWlr->mNReferences++;
  }
  //***************************************************************************
  
  
#ifdef USE_OLD_TOKEN_PASSING
  //***************************************************************************
  //***************************************************************************
  int 
  PassTokenMaxForLattices(Token* from, Token* to, FLOAT mLike)
  {
    int ret = 0;
    
#ifdef TRACE_TOKENS
    printf("(%.2f + %.2f -> %.2f = ", from->mLike, mLike, to->mLike);
#endif
    
    // Since we generate lattices, WLR is created for each node and therefor 
    // WLR for token 'from' should have been just created
    assert(from->mpAltHyps == NULL);
    
    if (!to->IsActive() || from->mLike + mLike > to->mLike) 
    {
      Token::AltHypList * pAltHyps = NULL;
      
      if(to->IsActive()) 
      {
        to->AddAlternativeHypothesis(to->mpWlr);
        pAltHyps = to->mpAltHyps;
      }
      
      KillToken(to);
  
      ret = 1;
      *to = *from;
      to->mLike += mLike;
      to->mpAltHyps = pAltHyps;
      
      if (to->mpWlr) 
        to->mpWlr->mNReferences++;
    } else {
      to->AddAlternativeHypothesis(from->mpWlr);
    }

#ifdef TRACE_TOKENS
    printf("%.2f)\n", to->mLike);
#endif

    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  PassTokenMax(Token* from, Token* to, FLOAT mLike)
  {
    assert(from->mpAltHyps == NULL);
    int ret = 0;
    
  #ifdef TRACE_TOKENS
    printf("(%.2f + %.2f -> %.2f = ", from->mLike, mLike, to->mLike);
  #endif
    
    if (!to->IsActive() || from->mLike + mLike > to->mLike) 
    {

      KillToken(to);
  
      ret = 1;
      *to = *from;
      to->mLike += mLike;
      
      if (to->mpWlr) 
        to->mpWlr->mNReferences++;
    }
  #ifdef TRACE_TOKENS
    printf("%.2f)\n", to->mLike);
  #endif
    return ret;
  }
  
  /*int PassTokenSum(Token *from, Token *to, FLOAT mLike)
  {
    double tl;
    int ret = 0;
  #ifdef TRACE_TOKENS
    printf("(%.2f + %.2f -> %.2f = ", from->mLike, mLike, to->mLike);
  #endif
    if (IS_ACTIVE(*to)) {
      tl = LogAdd(to->mLike, from->mLike + mLike);
    } else {
      tl = from->mLike + mLike;
    }
  
    if (!IS_ACTIVE(*to) || from->mLike + mLike > to->mBestLike) {
      KillToken(to);
  
      ret = 1;
      *to = *from;
      to->mBestLike = from->mLike + mLike;
      if (to->mpWlr) {
        to->mpWlr->mNReferences++;
      }
    }
  
    to->mLike = tl;
  #ifdef TRACE_TOKENS
    printf("%.2f)\n", to->mLike);
  #endif
    return ret;
  }*/
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  PassTokenSum(Token* from, Token* to, FLOAT like)
  {
    double        tl;
    FloatInLog    fe;
    FloatInLog    fil_from_like = {like, 0};
    int           ret = 0;
    
#ifdef TRACE_TOKENS
    printf("(%.2f + %.2f -> %.2f = ", from->mLike, like, to->mLike);
#endif
    
    if (to->IsActive()) 
    {
      tl = LogAdd(to->mLike, from->mLike + like);
      fe = FIL_Add(to->mAccuracy, FIL_Mul(from->mAccuracy, fil_from_like));
    } 
    else 
    {
      tl = from->mLike + like;
      fe = FIL_Mul(from->mAccuracy, fil_from_like);
    }
  
    if (!to->IsActive() || from->mLike + like > to->mBestLike) 
    {
      KillToken(to);
  
      ret = 1;
      *to = *from;
      to->mBestLike = from->mLike + like;
      
      if (to->mpWlr) 
      {
        to->mpWlr->mNReferences++;
      }
    }
  
    to->mLike     = tl;
    to->mAccuracy = fe;
#ifdef TRACE_TOKENS
    printf("%.2f)\n", to->mLike);
#endif
    return ret;
  }

#else


  //***************************************************************************
  //***************************************************************************
  int 
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
    
    //printf("PassTokenMaxForLattices: totalLike = %f  acousticLike = %f\n", totalLike, acousticLike);
    
    if (!pTo->IsActive() || pFrom->mLike + totalLike > pTo->mLike) 
    {
      AltHypList* p_alt_hyps = NULL;
      
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
  int 
  PassTokenMax(Token* pFrom, Token* pTo, FLOAT totalLike, 
      FLOAT acousticLike)
  {
    assert(pFrom->mpAltHyps == NULL);
    int ret = 0;
    
  #ifdef TRACE_TOKENS
    printf("(%.2f + %.2f -> %.2f = ", pFrom->mLike, mLike, to->mLike);
  #endif
    
    //printf("PassTokenMax: totalLike = %f  acousticLike = %f\n", totalLike, acousticLike);

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
  int 
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
  
#endif // USE_OLD_TOKEN_PASSIGN  



  //***************************************************************************
  //***************************************************************************
  void 
  FreeWordLinkRecords(WordLinkRecord* wlr)
  {
    if (wlr != NULL) 
    {
      --wlr->mNReferences;
      assert(wlr->mNReferences >= 0);
  
      if (wlr->mNReferences == 0) 
      {
        FreeWordLinkRecords(wlr->mpNext);
        
        if (wlr->mpAltHyps)
        {
          for(AltHypList::iterator i = wlr->mpAltHyps->begin(); i != wlr->mpAltHyps->end(); i++)
            FreeWordLinkRecords(i->mpWlr);
        
          delete wlr->mpAltHyps;
        }
        
#ifdef DEBUG_MSGS
        assert(wlr->mIsFreed == false);
        wlr->mIsFreed = true;
#else
        free(wlr);
#endif
      }
    }
  }

  
  
  //***************************************************************************
  //***************************************************************************
  void 
  KillToken(Token* token)
  {
    token->mLike = LOG_0;
    token->mAcousticLike = LOG_0;

    FreeWordLinkRecords(token->mpWlr);
    token->mpWlr  = NULL;
    token->mpTWlr = NULL;
  }

  

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
  //***************************************************************************
  FLOAT 
  DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs, Decoder* net)  
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
  FLOAT 
  DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs, Decoder* net, 
      size_t nOverhead)  
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
  FLOAT 
  Decoder::
  DiagCGaussianDensity(const Mixture* mix, const FLOAT* pObs)
  {
    if (mpMixPCache[mix->mID].mTime != mTime)
    {
      FLOAT m_like;
      
      // we call the computation
      m_like = compute_diag_c_gaussian_density(
        pObs,
        mix->GConst(),
        mix->mpMean->mVector.cpData(),
        mix->mpVariance->mVector.cpData(),        
        mix->mpMean->VectorSize());
    
      mpMixPCache[mix->mID].mTime  = mTime;
      mpMixPCache[mix->mID].mValue = m_like;
      
      return m_like;
    }
    else
    {
      return mpMixPCache[mix->mID].mValue;
    }
  }

      
#ifdef DEBUG_MSGS
  int gaus_computaions = 0;
#endif
  
  
  //***************************************************************************
  //***************************************************************************
  void
  FindNBestMixtures(State* pState, FLOAT* pObs, NBestRecord* pNBest, 
      size_t nBest, int time)
  {
    // go through all mixtures and find N-best
    for (size_t i = 0; i < pState->mNMixtures; i++)
    {
      size_t    j;
      Mixture*  mix   = pState->mpMixture[i].mpEstimates;
      FLOAT*    l_obs = XformPass(mix->mpInputXform, pObs, time, FORWARD);
      assert(l_obs != NULL);
      FLOAT     glike = DiagCGaussianDensity(mix, l_obs, NULL) + pState->mpMixture[i].mWeight;

      // the new result is better than some of the already chosen ones
      for (j = 0; j < nBest && pNBest[j].like < glike; j++)
      {}

      if (j > 0)
      {
        // rearange
        for (size_t k = 1; k < j; k++)
          pNBest[k - 1] = pNBest[k];

        pNBest[j - 1].like  = glike;
        pNBest[j - 1].index = i;
      }
    }
  }


  //***************************************************************************
  //***************************************************************************
  FLOAT 
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

  
/*  
  // ***************************************************************************
  // ***************************************************************************
  FLOAT 
  DiagCGaussianMixtureDensity(State* state, FLOAT* pObs, Decoder* net)
  {
    size_t  i;
    size_t  c; // cache counter
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
      const size_t cache_offset  = net->mTime % OUT_P_CACHES;
      const size_t states        = net->mpModelSet->mNStates;
      const long   orig_time     = net->mTime;
      
      // precount probs
      for (c = 0; c < OUT_P_CACHES - cache_offset; c++)
      {
        for (i = 0; i < state->mNMixtures; i++) 
        {
          FLOAT    glike;
          Mixture* mix = state->mpMixture[i].mpEstimates;
      
          l_obs   = XformPass(mix->mpInputXform, pObs,
                          net ? net->mTime : UNDEF_TIME,
                          net ? net->mPropagDir : FORWARD);
          
          assert(l_obs != NULL);
          
          glike   = DiagCGaussianDensity(mix, l_obs, net);
          m_like  = LogAdd(m_like, glike + state->mpMixture[i].mWeight);
        }
        
#ifdef DEBUG_MSGS
        gaus_computaions++;
#endif
      
        if (net)
        {
          net->mpOutPCache[state->mID + c * states].mTime = net->mTime;
          net->mpOutPCache[state->mID + c * states].mValue = m_like;
        }
        else
        {        
          return m_like;
        }
        
        net->mTime++;
      }
      
      net->mTime = orig_time;
      return  net->mpOutPCache[state->mID].mValue = m_like;
    }
  }
*/
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  Decoder::
  DiagCGaussianMixtureDensity(State* state, FLOAT* pObs)
  {
    size_t  i;
    FLOAT   m_like = LOG_0;
    FLOAT*  l_obs;
    assert(state->mOutPdfKind == KID_DiagC);
  
        
    //if (!net || net->mpOutPCache[state->mID].mTime != net->mTime) 
    if (mpOutPCache[state->mID].mTime == mTime) 
    {
      return mpOutPCache[state->mID].mValue;
    }
    else
    {
      for (i = 0; i < state->mNMixtures; i++) 
      {
        FLOAT    g_like;
        Mixture* mix = state->mpMixture[i].mpEstimates;
    
        l_obs = XformPass(mix->mpInputXform, pObs, mTime, mPropagDir);
        
        assert(l_obs != NULL);
        g_like = DiagCGaussianDensity(mix, l_obs);
        m_like  = LogAdd(m_like, g_like + state->mpMixture[i].mWeight);
      }
    
  #ifdef DEBUG_MSGS
      gaus_computaions++;
  #endif
    
      mpOutPCache[state->mID].mTime = mTime;
      mpOutPCache[state->mID].mValue = m_like;
      
      return m_like;
    }
    
    /*
    if (net)
    {
      // compute which offset to use
      const size_t cache_segment = state->mID * OUT_P_CACHES;
      const size_t cache_offset  = net->mTime & (OUT_P_CACHES -1);
      const size_t cache_index   = cache_segment + cache_offset;
      const size_t stride        = net->mpModelSet->mInputVectorStride;      
      
      if (net->mpOutPCache[cache_index].mTime == net->mTime)
      {
        return net->mpOutPCache[cache_index].mValue;
      }
      else
      {
        long orig_time = net->mTime;
        
        for (size_t ci = cache_offset; ci < OUT_P_CACHES; ci++)
        {
          m_like = LOG_0;
          for (i = 0; i < state->mNMixtures; i++) 
          {
            FLOAT    glike;
            Mixture* mix = state->mpMixture[i].mpEstimates;
        
            l_obs = XformPass(mix->mpInputXform, pObs,
                              net->mTime, net->mPropagDir);
        
            assert(l_obs != NULL);
            glike = DiagCGaussianDensity(mix, l_obs, net);
            m_like  = LogAdd(m_like, glike + state->mpMixture[i].mWeight);
          }
        
#ifdef DEBUG_MSGS
          gaus_computaions++;
#endif
        
          // store in cache
          net->mpOutPCache[cache_index].mTime = net->mTime;
          net->mpOutPCache[cache_index].mValue = m_like;
          
          // move one frame ahead
          pObs += stride;
          net->mTime ++;
        }
        
        // restore the original network time
        net->mTime =  orig_time;
      }
      
      return net->mpOutPCache[cache_index].mValue;
    }
    else
    {
      for (i = 0; i < state->mNMixtures; i++) 
      {
        FLOAT    glike;
        Mixture* mix = state->mpMixture[i].mpEstimates;
    
        l_obs = XformPass(mix->mpInputXform, pObs,
                        UNDEF_TIME, FORWARD);
    
        assert(l_obs != NULL);
        glike = DiagCGaussianDensity(mix, l_obs, NULL);
        m_like  = LogAdd(m_like, glike + state->mpMixture[i].mWeight);
      }
      return m_like;
    }
    */    
  }
  
  
    
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  FromObservationAtStateId(State *state, FLOAT *obs, Decoder *net)
  {
    obs = XformPass(net->mpModelSet->mpInputXform, obs,
                    net ? net->mTime : UNDEF_TIME,
                    net ? net->mPropagDir : FORWARD);
    assert(obs != NULL);
    return obs[state->PDF_obs_coef];
  }

  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  Decoder::
  FromObservationAtStateId(State* pState, FLOAT* pObs)
  {
    pObs = XformPass(mpModelSet->mpInputXform, pObs, mTime, mPropagDir);
    assert(pObs != NULL);
    return pObs[pState->PDF_obs_coef];
  }
    
  
  
  //***************************************************************************
  //***************************************************************************
  #ifdef DEBUG_MSGS
  void 
  PrintNumOfWLR() 
  {
    WordLinkRecord *  wlr = firstWLR;
    int               nwlrs = 0;
    int               nwlrsnf = 0;
    
    while (wlr) 
    {
      nwlrs++;
      
      if (!wlr->mIsFreed)
        nwlrsnf++;
      
      wlr = wlr->mpTmpNext;
    }
    printf("%d Released: %d\n", nwlrsnf, nwlrs - nwlrsnf);
  }
  #endif
  
  
  
  //***************************************************************************
  //***************************************************************************
  /*void BaumWelchInit(ModelSet *pHmms) 
  {
    //ResetAccumsForHMMSet(pHmms);
    pHmms.ResetAccums();
  }*/
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  UpdateXformStatCache(
    XformStatCache*  xfsc,
    Xform*           topXform, //to locate positions in input vector
    FLOAT*           input) 
  {
    size_t          i;
    size_t          j;
    const size_t    size(xfsc->mpXform->mInSize);
  
    if (topXform == NULL)
      return;
  
    if (topXform == xfsc->mpXform) 
    {
      if (xfsc->mNorm > 0) 
      {
        for (i=0; i < size; i++) 
        {
          xfsc->mpStats[i] += input[i];
          for (j=0; j <= i; j++) 
          {
            xfsc->mpStats[size + i*(i+1)/2 + j] += input[i] * input[j];
          }
        }
      } 
      else 
      {
        for (i=0; i < size; i++) 
        {
          xfsc->mpStats[i] = input[i];
          for (j=0; j <= i; j++) 
          {
            xfsc->mpStats[size + i*(i+1)/2 + j] = input[i] * input[j];
          }
        }
      }
      xfsc->mNorm++;
    } 
    
    else if (topXform->mXformType == XT_COMPOSITE) 
    {
      CompositeXform* cxf = (CompositeXform *) topXform;
      for (i=0; i<cxf->mpLayer[0].mNBlocks; i++) 
      {
        UpdateXformStatCache(xfsc, cxf->mpLayer[0].mpBlock[i], input);
        input += cxf->mpLayer[0].mpBlock[i]->mInSize;
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  UpdateXformInstanceStatCaches(XformInstance* xformInstance,
                                FLOAT* pObservation, int time)
  {
    size_t  i;
    size_t  j;
    FLOAT*  obs;
  
    if (xformInstance == NULL || xformInstance->mStatCacheTime == time) return;
  
    xformInstance->mStatCacheTime = time;
  
    if (xformInstance->mNumberOfXformStatCaches == 0) return;
  
    if (xformInstance->mpInput) 
    {
      UpdateXformInstanceStatCaches(xformInstance->mpInput, pObservation, time);
    }
  
    for (i = 0; i < xformInstance->mNumberOfXformStatCaches; i++) 
    {
      XformStatCache* xfsc = &xformInstance->mpXformStatCache[i];
  
      //just link to upper level?
      if (xfsc->mpUpperLevelStats != NULL 
      &&  xfsc->mpUpperLevelStats->mpStats == xfsc->mpStats) 
      { 
        xfsc->mNorm = xfsc->mpUpperLevelStats->mNorm;
        continue;
      }
  
      obs = XformPass(xformInstance->mpInput, pObservation, time, FORWARD);
      xfsc->mNorm = 0;
      UpdateXformStatCache(xfsc, xformInstance->mpXform, obs);
      
      if (xfsc->mpUpperLevelStats != NULL) 
      {
        size_t size = xfsc->mpXform->mInSize;
        for (j = 0; j < size + size*(size+1)/2; j++) 
        {
          xfsc->mpStats[j] += xfsc->mpUpperLevelStats->mpStats[j];
        }
        xfsc->mNorm += xfsc->mpUpperLevelStats->mNorm;
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  ReestState(
    NetworkType::NodeType*         pNode,
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
    State*        state  = pNode->mpHmm->        mpState[stateIndex];
    State*        state2 = pNode->mpHmmToUpdate->mpState[stateIndex];
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
        FLOAT  bjmtO  = ::DiagCGaussianDensity(mix, xobs, this);
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
        FLOAT     bjmtO = ::DiagCGaussianDensity(mix, xobs, this);
        
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
        Lqjmt = exp(Lqjmt) * updateDir;
  
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
        
//      if (mmi_den_pass) {
//        state2->mpMixture[m].mWeightAccumDen += Lqjmt;
//      } else {
        state2->mpMixture[m].mWeightAccum     += Lqjmt; //  Update weight accum
//      }
  
        if (Lqjmt < 0)
          state2->mpMixture[m].mWeightAccumDen -= Lqjmt;
  
        if (ixf == NULL || ixf->mNumberOfXformStatCaches == 0) 
          continue;
  
        UpdateXformInstanceStatCaches(ixf, pObs2, mTime);
  
        for (i = 0; i < ixf->mNumberOfXformStatCaches; i++) 
        {
          XformStatCache *xfsc = &ixf->mpXformStatCache[i];
          Variance *var  = state2->mpMixture[m].mpEstimates->mpVariance;
          Mean     *mean = state2->mpMixture[m].mpEstimates->mpMean;
  
          for (j = 0; j < var->mNumberOfXformStatAccums; j++) 
          {
            XformStatAccum *xfsa = &var->mpXformStatAccum[j];
            
            if (xfsa->mpXform == xfsc->mpXform) 
            {
              size_t size = xfsc->mpXform->mInSize;
              for (k = 0; k < size+size*(size+1)/2; k++) 
              {
                xfsa->mpStats[k] += xfsc->mpStats[k] * Lqjmt;
              }
              xfsa->mNorm += xfsc->mNorm * Lqjmt;
              break;
            }
          }
  
          for (j = 0; j < mean->mNumberOfXformStatAccums; j++) 
          {
            XformStatAccum *xfsa = &mean->mpXformStatAccum[j];
            
            if (xfsa->mpXform == xfsc->mpXform) 
            {
              size_t size = xfsc->mpXform->mInSize;
              for (k = 0; k < size; k++) 
              {
                xfsa->mpStats[k] += xfsc->mpStats[k] * Lqjmt;
              }
              xfsa->mNorm += xfsc->mNorm * Lqjmt;
              break;
            }
          }
        }
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  /*FLOAT *
  StateOccupationProbability(Decoder *net, FLOAT *obsMx, ModelSet *hmms,
                             int nFrames, FLOAT **outProbOrMahDist, int getMahalDist)
  {
    int i, j, k;
    int nNetModels;
    int nEmitingStates;
    FLOAT totalLike;
    FLOAT *occupProb;
  
    FLOAT *beta = (FLOAT *) malloc(net->mNumberOfNetStates * (nFrames+1) * sizeof(FLOAT));
    FLOAT *alfa = (FLOAT *) malloc(net->mNumberOfNetStates * (nFrames+1) * sizeof(FLOAT));
    Cache *outPCache = (Cache *) malloc(nFrames * hmms->mNStates * sizeof(Cache));
  
    if (outPCache == NULL || beta == NULL || alfa == NULL) {
      Error("Insufficient memory");
    }
  
    for (i = 0; i < net->mNumberOfNetStates * (nFrames+1); i++) {
      beta[i] = alfa[i] = LOG_0;
    }
  
    nNetModels = 0;
    for (i=0; i < net->nnodes; i++) {
      if (net->mpNodes[i].type & NT_MODEL) nNetModels++;
    }
  
    nEmitingStates = net->mNumberOfNetStates - 2 * nNetModels;
  
    for (i = 0; i < nFrames * hmms->mNStates; i++) {
      outPCache[i].mTime = UNDEF_TIME;
      outPCache[i].mValue = LOG_0;
    }
  
    free(net->mpOutPCache);
  
  
    occupProb = (FLOAT *) malloc(nEmitingStates * nFrames * sizeof(FLOAT));
    if (occupProb == NULL) {
      Error("Insufficient memory");
    }
  
    if (outProbOrMahDist != NULL) {
      *outProbOrMahDist = (FLOAT *) malloc(nEmitingStates * nFrames * sizeof(FLOAT));
      if (*outProbOrMahDist == NULL) {
        Error("Insufficient memory");
      }
    }
  
    net->PassTokenInModel   = &PassTokenSum;
    net->PassTokenInNetwork = &PassTokenSum;
    net->mAlignment          = NO_ALIGNMENT;
  
    //Backward Pass
    net->mPropagDir = BACKWARD;
  
    net->mTime = nFrames;
    TokenPropagationInit(net,
                        beta + net->mNumberOfNetStates * net->mTime,
                        NULL);
  
    for (i = nFrames-1; i >= 0; i--) {
  
      net->mpOutPCache = outPCache + hmms->mNStates * i;
      TokenPropagationInModels(net,  obsMx + hmms->mInputVectorSize * i,
                              beta + net->mNumberOfNetStates * net->mTime,
                              NULL);
  
      if (outProbOrMahDist != NULL && getMahalDist !=4) {
        int state_counter = 0;
  
        for (k=0; k < net->nnodes; k++) {
          NetworkType::NodeType *node = &net->mpNodes[k];
          if (node->mType & NT_MODEL) {
            for (j = 0; j < node->mpHmm->mNStates - 2; j++, state_counter++) {
              FLOAT tmpf = net->OutputProbability(node->mpHmm->mpState[j],
                                        obsMx + hmms->mInputVectorSize * i, net);
              switch (getMahalDist) {
                case 0:
                  break;
                case 1:
                tmpf = log((tmpf / -0.5) - node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->GConst()) * -0.5;
                break;
                case 2:
                tmpf = log((tmpf / -0.5) - node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->GConst()) * -1;
                break;
                case 3:
                tmpf += node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->GConst() * 0.5;
  //               tmpf /= hmms->mInputVectorSize;
                break;
              }
  
              (*outProbOrMahDist)[i * nEmitingStates + state_counter] = tmpf;
            }
          }
        }
      }
  
      net->mTime--;
      TokenPropagationInNetwork(net,
                                beta + net->mNumberOfNetStates * net->mTime, NULL);
    }
  
    if (!IS_ACTIVE(*net->rNetwork().pFirst()->mpExitToken)) {
      TokenPropagationDone(net);
      net->mpOutPCache = outPCache;
      free(alfa);
      for (i=0; i < net->mNumberOfNetStates * nFrames; i++) {
        beta[i] = LOG_0;
      }
      return beta; // No token survivered
    }
  
    totalLike = net->rNetwork().pFirst()->mpExitToken->mLike;
  //  totalLike = 0;
    TokenPropagationDone(net);
  
    //Forward Pass
    net->mPropagDir = FORWARD;
    net->mTime = 0;
    TokenPropagationInit(net,
                        beta + net->mNumberOfNetStates * net->mTime,
                        alfa + net->mNumberOfNetStates * net->mTime);
  
    for (i = 0; i < nFrames; i++) {
      net->mpOutPCache = outPCache + hmms->mNStates * i;
      net->mTime++;
      TokenPropagationInModels(net,  obsMx + hmms->mInputVectorSize * i,
                              beta + net->mNumberOfNetStates * net->mTime,
                              alfa + net->mNumberOfNetStates * net->mTime);
  
      TokenPropagationInNetwork(net,
                                beta + net->mNumberOfNetStates * net->mTime,
                                alfa + net->mNumberOfNetStates * net->mTime);
  
      if (outProbOrMahDist != NULL && getMahalDist == 4) {
        for (j=0; j < nEmitingStates; j++) {
          (*outProbOrMahDist)[i * nEmitingStates + j] = totalLike;
        }
      }
    }
  
    net->mpOutPCache = outPCache;
  
    TokenPropagationDone(net);
  
    for (i = 1; i <= nFrames; i++) {
      int state_counter = 0;
  
      for (k=0; k < net->nnodes; k++) {
        NetworkType::NodeType *node = &net->mpNodes[k];
        if (node->mType & NT_MODEL) {
          for (j = 0; j < node->mpHmm->mNStates - 2; j++, state_counter++) {
            int idx = (net->mNumberOfNetStates * i) + node->mEmittingStateId + j + 1;
            FLOAT occp = beta[idx] + alfa[idx] - totalLike;
            occupProb[(i-1) * nEmitingStates + state_counter] = occp < LOG_MIN ?
                                                                LOG_0 : occp;
          }
        }
      }
    }
  
    free(alfa); free(beta);
    return occupProb;
  }*/
  
  

  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::TokenPropagationInit()
  {
    size_t                    i;
    NetworkType::NodeType*   p_node;
  
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
    mActiveTokens = 0;
    //  mTime = 0;
  
    if (mCollectAlphaBeta && InForwardPass()) 
    {
      if (mCompactRepresentation)
        Error("Fatal: CSTK format used for forward-backward");
        
      for (p_node = rNetwork().pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      {
        if (!(p_node->mType & NT_MODEL)) 
          continue;
          
        p_node->mpAlphaBetaList = p_node->mpAlphaBetaListReverse = NULL;
      }
    }
    
    // Needed to load last FWBWRs to mpAlphaBetaList
    if (mCollectAlphaBeta && !InForwardPass()) 
    {
      if (mCompactRepresentation)
        Error("Fatal: CSTK format used for forward-backward");

      for (p_node = rNetwork().pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      {
        if (!(p_node->mType & NT_MODEL)) 
          continue;
          
        BackwardPruning(mTime, p_node, p_node->mpHmm->mNStates-1);
      }
    }
    
    // Last p_node is not activated and deactivated in usual way; allocate ActiveNodeRecord here
    p_node = InForwardPass() ?  rNetwork().pLast() : rNetwork().pFirst();
    p_node->mpAnr = new ActiveNodeRecord(p_node);

    // First p_node is also not activated and deactivated in usual way...
    p_node = InForwardPass() ? rNetwork().pFirst() : rNetwork().pLast();
    mpActiveNodes = p_node;
    p_node->mpAnr = new ActiveNodeRecord(p_node);
    p_node->mpAnr->mpTokens[0].mLike = 0.0;
    p_node->mpAnr->mpTokens[0].mAcousticLike = 0.0;
    p_node->mpAnr->mpTokens[0].mAccuracy.logvalue = LOG_0;
    p_node->mpAnr->mpTokens[0].mAccuracy.negative = 0;
//    p_node->mpAnr->mpTokens[0].mpWlr = NULL; 
#ifdef bordel_staff
//    p_node->mpTokens[0].mpTWlr = NULL;
    p_node->mpAnr->mpTokens[0].mBestLike = 0.0;
#endif
  
    p_node->mpAnr->mpPrevActiveNode = p_node->mpAnr->mpNextActiveNode = NULL;
    
    //Is this needed?
    p_node->mpAnr->mActiveNodeFlag = 1;

    MarkWordNodesLeadingFrom(p_node);
    pActivateWordNodesLeadingFrom(p_node);
    TokenPropagationInNetwork();
    DeactivateWordNodesLeadingFrom(p_node);
    
//    p_node->mpAnr->mActiveNodeFlag = 0;
    
    if (p_node->mpAnr->mpPrevActiveNode) 
    {
      p_node->mpAnr->mpPrevActiveNode->mpAnr->mpNextActiveNode = p_node->mpAnr->mpNextActiveNode;
    }
    delete p_node->mpAnr;
    p_node->mpAnr = NULL;    
  } 
  // Decoder::TokenPropagationInit()
  //***************************************************************************
  


  //***************************************************************************
  //***************************************************************************
  /*
  void
  Decoder::
  GenerateToken(NetworkType::NodeType* pNode)
  {
    // activate the node
    ActivateNode(pNode);
    
    // generate the token
    pNode->mpAnr->mpTokens[0].mLike = 0.0;
    pNode->mpAnr->mpTokens[0].mAcousticLike = 0.0;
    pNode->mpAnr->mpTokens[0].mAccuracy.logvalue = LOG_0;
    pNode->mpAnr->mpTokens[0].mAccuracy.negative = 0;

#ifdef bordel_staff
    pNode->mpAnr->mpTokens[0].mBestLike = 0.0;
#endif
  
    // propagate the token further
    TokenPropagationInNetwork();

    DeactivateNode(pNode);
    
    return;
    DeactivateWordNodesLeadingFrom(pNode);
    
    
    if (pNode->mpAnr->mpPrevActiveNode) 
    {
      pNode->mpAnr->mpPrevActiveNode->mpAnr->mpNextActiveNode = pNode->mpAnr->mpNextActiveNode;
    }

    delete pNode->mpAnr;
    pNode->mpAnr = NULL;    
  }
  */
  // GenerateToken(NetworkType::NodeType* pNode)
  //***************************************************************************



  //***************************************************************************
  //***************************************************************************
  /*void 
  Decoder::
  AddLinkToLattice(NetworkType::NodeType *from, NetworkType::NodeType *to, FLOAT lmLike)
  {
    assert(from->mpAnr && from->mpAnr && from->mpAnr->mpExitToken->mpWlr);
    long long from_time = from->mpAnr->mpExitToken->mpWlr->mpNext 
                          ? from->mpAnr->mpExitToken->mpWlr->mpNext->mTime : 0;
    
    std::ostringstream ss;
    ss << ios::hex
       << from
       << from_time;
               
    NetworkType::NodeType * lattice_from = find_or_create_node(&mLatticeNodeHash, ss.str().c_str(), &mpLatticeLastNode);
    lattice_from->mType  = from->mType;
    lattice_from->mpName = from->mpName; // !!! Relys on union, stays for all mpName, mpHmm and mpPronun

    ss.str("");
    ss << ios::hex 
       << to
       << from->mpAnr->mpExitToken->mpWlr->mTime;

    NetworkType::NodeType * lattice_to = find_or_create_node(&mLatticeNodeHash, ss.str().c_str(), &mpLatticeLastNode);
    lattice_to->mType  = to->mType;
    lattice_to->mpName = to->mpName; // !!! Relys on union, stays for all mpName, mpHmm and mpPronun
    
    int nl = ++lattice_from->mNLinks;

    lattice_from->mpLinks = (NetworkType::LinkType *) realloc(lattice_from->mpLinks, nl * sizeof(NetworkType::LinkType));
    if (lattice_from->mpLinks == NULL) Error("Insufficient memory");
          
    lattice_from->mpLinks[nl-1].mpNode = lattice_to;
    lattice_from->mpLinks[nl-1].mLike = lmLike;
    lattice_from->mpLinks[nl-1].mpNode->mNBackLinks++;
    
  }*/
  // AddLinkToLattice(NetworkType::NodeType *from, NetworkType::NodeType *to, FLOAT lmLike)
  //***************************************************************************
  

  
  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  TokenPropagationInNetwork()
  {
    NetworkType::NodeType*   p_node;
    NetworkType::LinkType*   links;
    int                       i;
    int                       n_links;
  
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
                  ? mpBestToken->mLike - mPruningThresh 
                  : LOG_MIN;
  
    p_node = InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst();

    assert(p_node->mpAnr != NULL);
    KillToken(p_node->mpAnr->mpExitToken);

//  NetworkType::NodeType *Xnode = mpActiveNodes;
/*  for (p_node = InForwardPass() ? rNetwork().pFirst() : rNetwork().pLast();
      p_node != NULL;
      p_node = InForwardPass() ? p_node->mpNext : p_node->mpBackNext)  // */
    for (p_node = mpActiveNodes; p_node != NULL; p_node = p_node->mpAnr->mpNextActiveNode) 
    {
      assert(p_node->mpAnr != NULL);

      // If tee model is entered
      if ((p_node->mType & NT_TEE) && p_node->mpAnr->mpTokens[0].IsActive()) 
      {
//      assert(p_node->mActiveNodeFlag || (p_node->mType & NT_TEE && p_node->mIsActive));
//      for (Xnode = mpActiveNodes; Xnode && Xnode != p_node; Xnode = Xnode->mpNextActiveNode);
//      assert(Xnode);
    
        if (!(mCollectAlphaBeta && !InForwardPass()
        &&   BackwardPruning(mTime, p_node, 0))) // backward pruning after forward pass
        {
          Hmm*    p_hmm       = p_node->mpHmm;
          FLOAT   trans_prob  = p_hmm->mpTransition->mpMatrixO[p_hmm->mNStates - 1] *
            mTranScale;

#ifdef TRACE_TOKENS
          printf("Tee model State 0 -> Exit State ");
#endif
          if (mLatticeGeneration)
          //||  p_node->mpAnr->mpTokens[0].mLike > mBeamThresh)
          {
            p_node->mpAnr->mpTokens[0].AddWordLinkRecord(rNetwork().pFirst(), 
                -1, mTime);
          }
            
          PassTokenInModel(&p_node->mpAnr->mpTokens[0], 
              p_node->mpAnr->mpExitToken, trans_prob, trans_prob);
        }
      }
  
      if (p_node->mpAnr->mpExitToken->IsActive()) 
      {
//      assert(p_node->mActiveNodeFlag || (p_node->mType & NT_TEE && p_node->mIsActive));
//      for (Xnode = mpActiveNodes; Xnode && Xnode != p_node; Xnode = Xnode->mpNextActiveNode);
//      assert(Xnode);
        if (p_node->mType & NT_MODEL) 
        {
          if (mCollectAlphaBeta) 
          {
            if (InForwardPass()) 
              WriteAlpha(mTime, p_node, p_node->mpHmm->mNStates-1, p_node->mpAnr->mpExitToken);
            else 
              WriteBeta(mTime, p_node, 0, p_node->mpAnr->mpExitToken);
          }
          p_node->mpAnr->mpExitToken->mLike += mMPenalty;
        }
        else if (p_node->mType & NT_WORD && p_node->mpPronun != NULL) 
        {
          p_node->mpAnr->mpExitToken->mLike += mWPenalty +
                                             mPronScale * p_node->mpPronun->prob;
          
          /*if (p_node->mpExitToken->mLike < mWordThresh) {
            p_node->mpExitToken->mLike = LOG_0;
          }*/
        }
  
        if (p_node->mpAnr->mpExitToken->mLike > mBeamThresh) 
        {
          if (p_node->mType & NT_WORD 
          && (p_node->mpPronun != NULL || mLatticeGeneration)
          && mAlignment & WORD_ALIGNMENT) 
          {
            p_node->mpAnr->mpExitToken->AddWordLinkRecord(p_node, -1, mTime);
          } 
          else if (p_node->mType & NT_MODEL && mAlignment & MODEL_ALIGNMENT) 
          {
            p_node->mpAnr->mpExitToken->AddWordLinkRecord(p_node, -1, mTime);
          }
  
          n_links = InForwardPass() ? p_node->mNLinks : p_node->mNBackLinks;
          links  = InForwardPass() ? p_node->mpLinks  : p_node->mpBackLinks;
  
          for (i = 0; i < n_links; i++) 
          {
            FLOAT acoustic_like = links[i].AcousticLike() * mTranScale;
            FLOAT lm_like       = (links[i].LmLike() - acoustic_like) * mLmScale;
            
            if (p_node->mpAnr->mpExitToken->mLike + lm_like > mBeamThresh 
                && (mCompactRepresentation
                ||((/*links[i].pNode()->Start() == UNDEF_TIME ||*/
                      links[i].pNode()->Start() <= mTime) 
                
                && (  links[i].pNode()->Stop()  == UNDEF_TIME        
                ||    links[i].pNode()->Stop()  >= mTime)
                      
                && (mSearchPaths != SP_TRUE_ONLY
                || (p_node->mType & NT_TRUE)   
                || !(links[i].pNode()->mType & NT_MODEL)))))
            {
              if (links[i].pNode()->mType & NT_MODEL) 
              {
                ActivateModel(links[i].pNode());
              } 
              else 
              {
                assert(links[i].pNode()->mpAnr->mActiveNodeFlag ||
                       links[i].pNode() == (InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst()));
              }

              assert(links[i].pNode()->mpAnr != NULL);

#ifdef TRACE_TOKENS
              printf("Node %d -> Node %d ", p_node->mAux, links[i].pNode()->mAux);
#endif

              // Current lattice generation algorithm expect that word link record is created for all tokens leaving
              // any network node (including model and !NULL nodes).
//              if (mLatticeGeneration)
//                AddLinkToLattice(p_node, links[i].pNode(), lm_like);
              
              PassTokenInNetwork(p_node->mpAnr->mpExitToken,
                  &links[i].pNode()->mpAnr->mpTokens[0], lm_like, acoustic_like);

            }
          }
        } 
        else if(mLatticeGeneration )
        {
          p_node->mpAnr->mpExitToken->AddWordLinkRecord(p_node, -1, mTime);
        }
        
        if (!(p_node->mType & NT_STICKY)) 
          KillToken(p_node->mpAnr->mpExitToken);
      }  
    }
  
  
    if (mCollectAlphaBeta) 
    {
      if (InForwardPass()) 
      {
        for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mpAnr->mpNextActiveModel) 
        {
          assert(p_node->mpAnr);
          if (/*!(p_node->mType & NT_MODEL) ||*/ p_node->mpAnr->mpTokens[0].mLike < LOG_MIN) 
            continue;
            
          WriteAlpha(mTime, p_node, 0, &p_node->mpAnr->mpTokens[0]);
        }
      } 
      else 
      {
        for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mpAnr->mpNextActiveModel) 
        {
          assert(p_node->mpAnr);
          if (/*!(p_node->mType & NT_MODEL) ||*/ p_node->mpAnr->mpTokens[0].mLike < LOG_MIN
            || BackwardPruning(mTime, p_node, p_node->mpHmm->mNStates-1)) 
            continue;
            
          WriteBeta(mTime, p_node, p_node->mpHmm->mNStates-1, &p_node->mpAnr->mpTokens[0]);          
        }
      }
    }
    
    if (mLatticeGeneration)
    {
      for (p_node = mpActiveModels; p_node != NULL; p_node = p_node->mpAnr->mpNextActiveModel)
      {
        if (p_node->mpAnr->mpTokens[0].IsActive())
        {
          p_node->mpAnr->mpTokens[0].AddWordLinkRecord(rNetwork().pFirst(), -1, mTime);
        }
      }

      if (rNetwork().pLast()->mpAnr != NULL 
      &&rNetwork().pLast()->mpAnr->mpExitToken->IsActive())
      {
        rNetwork().pLast()->mpAnr->mpExitToken->AddWordLinkRecord(rNetwork().pLast(), -1, mTime);
      }
    }
  
    //  Go through newly activeted models and costruct list of active nodes
    //  for (p_node=mpActiveModels; p_node&&!p_node->mIsActive; p_node=p_node->mpNextActiveModel) {
    //    ActivateNode(net, p_node);
    //  }
    assert(!HasCycle());
  }
  // TokenPropagationInNetwork()
  //***************************************************************************

  
  
  //***************************************************************************
  //***************************************************************************
  void
  Decoder::
  TokenPropagationInModels(FLOAT* pObservation)
  {
    NetworkType::NodeType*    p_node;
    NetworkType::NodeType*    p_next_active;
    Hmm*                      p_hmm;
    size_t                    winingToken = 0;
    size_t                    i;
    size_t                    j;
    int                       from;
    int                       to;
    int                       state_idx;
    //  FLOAT                 threshOutProb = LOG_0;
  
    mpBestToken = NULL;
    /*  if (mpThreshState) {
      threshOutProb = OutputProbability(mpThreshState, pObservation, net);
    } */
  
    for (p_node = mpActiveModels; p_node != NULL; p_node = p_next_active) 
    {
      // Store pointer to mpNextActiveModel here, because p_node->mpAnr can be dealocated in this block
      p_next_active = p_node->mpAnr->mpNextActiveModel;
      
      //  for (p_node = rNetwork().pFirst(); p_node != NULL; p_node = p_node->mpNext) 
      //    if (!(p_node->mType & NT_MODEL)) continue;
  
      assert(p_node->mpAnr);
      p_hmm = p_node->mpHmm;
      
      if (!mCompactRepresentation && (
          (/*p_node->Start() != UNDEF_TIME &&*/p_node->Start() >= mTime)
          ||  (  p_node->Stop()  != UNDEF_TIME &&  p_node->Stop()  <  mTime)
          ||  mSearchPaths == SP_TRUE_ONLY && !(p_node->mType & NT_TRUE)))
      {
        for (i = 0; i < p_hmm->mNStates-1; i++) 
        {
          KillToken(&p_node->mpAnr->mpTokens[i]);
        }
        
        DeactivateModel(p_node);
        continue;
      }
  
      if (mAccumType == AT_MPE && InForwardPass() && p_node->mpAnr->mpTokens[0].IsActive()) 
      {
        FloatInLog fil_lmpa =
          {p_node->mpAnr->mpTokens[0].mLike + log(fabs(p_node->mPhoneAccuracy)),
           p_node->mPhoneAccuracy < 0};
  
        p_node->mpAnr->mpTokens[0].mAccuracy = FIL_Add(p_node->mpAnr->mpTokens[0].mAccuracy, fil_lmpa);
      }
  
      // make a copy of active tokens into mpAuxTokens
      for (i = 0; i < p_hmm->mNStates-1; i++) 
      {
        assert(!mpAuxTokens[i].IsActive());
        if (p_node->mpAnr->mpTokens[i].IsActive()) 
        {
          mpAuxTokens[i] = p_node->mpAnr->mpTokens[i];
          assert(mpAuxTokens[i].mpAltHyps == NULL);
          p_node->mpAnr->mpTokens[i].mLike = LOG_0;
          p_node->mpAnr->mpTokens[i].mAcousticLike = LOG_0;
          p_node->mpAnr->mpTokens[i].mpWlr  = NULL;
        }
      }
  
      bool keep_model_active = false;
  
      assert(!p_node->mpAnr->mpTokens[p_hmm->mNStates-1].IsActive());
  
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
  
          assert(!mpAuxTokens[i].IsActive() || p_node->mpAnr->mIsActiveModel);
  
          if (p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] > LOG_MIN 
          &&  mpAuxTokens[i].IsActive()) 
          {
            FLOAT trans_prob = 
              p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] * mTranScale;
  
#ifdef TRACE_TOKENS
            printf("Model %d State %d -> State %d ",  p_node->mAux, (int) i, (int) j);
#endif

            if (PassTokenInModel(&mpAuxTokens[i], &p_node->mpAnr->mpTokens[j], 
                  trans_prob, trans_prob)) 
            {
              winingToken = i;
            }
          }
        }
  
        // if (IS_ACTIVE(p_node->mpTokens[j])) 
        if (p_node->mpAnr->mpTokens[j].mLike > mBeamThresh) 
        {
          FLOAT out_prob = OutputProbability(p_hmm->mpState[state_idx-1],
                                            pObservation, this);
          
          out_prob *= mOutpScale;
  
          /*if (out_prob < threshOutProb) {
            out_prob = threshOutProb;
          }*/
  
          if (mCollectAlphaBeta && !InForwardPass())
            WriteBeta(mTime, p_node, state_idx, &p_node->mpAnr->mpTokens[j]);
          
          if (mAccumType == AT_MFE && p_node->mType & NT_TRUE) 
          {
            FloatInLog fil_like = {p_node->mpAnr->mpTokens[j].mLike, 0};
            p_node->mpAnr->mpTokens[j].mAccuracy = FIL_Add(p_node->mpAnr->mpTokens[j].mAccuracy, fil_like);
          }
  
          p_node->mpAnr->mpTokens[j].mAccuracy.logvalue += out_prob;
          p_node->mpAnr->mpTokens[j].mLike              += out_prob;
          p_node->mpAnr->mpTokens[j].mAcousticLike      += out_prob;
  
          if (mCollectAlphaBeta && InForwardPass()) 
            WriteAlpha(mTime, p_node, state_idx, &p_node->mpAnr->mpTokens[j]);
  
          if (mAlignment & STATE_ALIGNMENT && winingToken > 0 &&
            (winingToken != j || mAlignment & FRAME_ALIGNMENT)) 
          {
            p_node->mpAnr->mpTokens[j].AddWordLinkRecord(
                p_node,
                (InForwardPass() ? winingToken : p_hmm->mNStates-1 - winingToken)-1,
                mTime-1);            
          }
  
          mActiveTokens++;
          keep_model_active = true;
          assert(p_node->mpAnr->mIsActiveModel);
        } 
        else 
        {
          assert(p_node->mpAnr->mIsActiveModel || !p_node->mpAnr->mpTokens[j].IsActive());
          KillToken(&p_node->mpAnr->mpTokens[j]);
        }
      }
  
      for (i = 0; i < p_hmm->mNStates-1; i++) 
      {
        KillToken(&mpAuxTokens[i]);
      }
  
      if (!keep_model_active) 
        DeactivateModel(p_node);
  
      state_idx = (InForwardPass() ? p_hmm->mNStates - 1 : 0);
  
      if (!keep_model_active ||
        (mCollectAlphaBeta && !InForwardPass() &&
          BackwardPruning(mTime-1, p_node, state_idx))) 
      {
        // backward pruning after forward pass
        continue;
        //KillToken(&p_node->mpTokens[p_hmm->mNStates - 1]);
      }
  
      for (i = 1; i < p_hmm->mNStates-1; i++) 
      {
        from = InForwardPass() ? i : 0;
        to   = InForwardPass() ? p_hmm->mNStates-1 : p_hmm->mNStates-1 - i;
  
        if (p_node->mpAnr->mpTokens[i].IsActive()) 
        {
          if (!mpBestToken || mpBestToken->mLike < p_node->mpAnr->mpTokens[i].mLike) 
          {
            mpBestToken = &p_node->mpAnr->mpTokens[i];
            mpBestNode  = p_node;
          }
  
          if (p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] > LOG_MIN) 
          {
            FLOAT trans_prob = p_hmm->mpTransition->mpMatrixO[from * p_hmm->mNStates + to] *
              mTranScale;
#ifdef TRACE_TOKENS
            printf("Model %d State %d -> Exit State ",  p_node->mAux, (int) i);
#endif
            if (PassTokenInModel(&p_node->mpAnr->mpTokens[i],
                                 &p_node->mpAnr->mpTokens[p_hmm->mNStates - 1],
                                 trans_prob, trans_prob)) 
            {
              winingToken = i;
            }
          }
        }
      }
  
      if (p_node->mpAnr->mpTokens[p_hmm->mNStates - 1].IsActive()) 
      {
        if (mAccumType == AT_MPE && !InForwardPass()) 
        {
          FloatInLog fil_lmpa =
            {p_node->mpAnr->mpTokens[p_hmm->mNStates - 1].mLike + log(fabs(p_node->mPhoneAccuracy)),
             p_node->mPhoneAccuracy < 0};
  
          p_node->mpAnr->mpTokens[p_hmm->mNStates - 1].mAccuracy =
            FIL_Add(p_node->mpAnr->mpTokens[p_hmm->mNStates - 1].mAccuracy, fil_lmpa);
        }
  
    //    ActivateNode(net, p_node);
        if (mAlignment & STATE_ALIGNMENT) 
        {
          p_node->mpAnr->mpTokens[p_hmm->mNStates - 1].AddWordLinkRecord(
              p_node,
              (InForwardPass() ? winingToken : p_hmm->mNStates-1 - winingToken-1)-1,
              mTime);
        }
      }
    }
    assert(!HasCycle());
  }
  // TokenPropagationInModels(FLOAT* pObservation)
  //***************************************************************************
  
  

  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  TokenPropagationDone()
  // {{{
  {
    int     j;
    NetworkType::NodeType*   node = InForwardPass() ? rNetwork().pLast() : rNetwork().pFirst();
  
    // Last node is not activated and deactivated in usual way; deallocate ActiveNodeRecord here
    assert(node->mpAnr != NULL);    
    KillToken(node->mpAnr->mpExitToken);
    delete node->mpAnr;
    node->mpAnr = NULL;
  
    while(mpActiveModels != NULL) 
    {
      for (j=0; j < mpActiveModels->mpHmm->mNStates; j++) 
      {
        KillToken(&mpActiveModels->mpAnr->mpTokens[j]);
      }
      
      DeactivateModel(mpActiveModels);
    }
    
//    if(mLatticeGeneration)
//      my_hdestroy_r(&mLatticeNodeHash, 1);
//
  } // TokenPropagationDone() }}}                 
  //***************************************************************************

    

  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  FreeFWBWRecords()
  {
    NetworkType::NodeType *node;
    
    // go through all nodes in the network and free the alpha/beta list
    // and alpha/beta reverse list
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (!(node->mType & NT_MODEL)) 
        continue;
  
      while (node->mpAlphaBetaList) 
      {
        FWBWR *fwbwr = node->mpAlphaBetaList;
        node->mpAlphaBetaList = fwbwr->mpNext;
        free(fwbwr);
      }
  
      while (node->mpAlphaBetaListReverse) 
      {
        FWBWR *fwbwr = node->mpAlphaBetaListReverse;
        node->mpAlphaBetaListReverse = fwbwr->mpNext;
        free(fwbwr);
      }
    }
  }

    
/*
  //***************************************************************************
  //***************************************************************************
  int 
  cmplnk(const void *a, const void *b)
  {
    return ((Link<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *) a)->pNode()->mAux - ((Link<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR> *) b)->pNode()->mAux;
  }

  //***************************************************************************
  //***************************************************************************
  void
  Decoder:: 
  SortNodes()
  {
    int       i;
    int       j;
    NetworkType::NodeType *    chain;
    NetworkType::NodeType *    last;
    NetworkType::NodeType *    node;
  
    // Sort nodes for forward (Viterbi) propagation
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      node->mAux = node->mNBackLinks;
    }
  
    for (i = 0; i < rNetwork().pFirst()->mNLinks; i++) 
    {
      rNetwork().pFirst()->mpLinks[i].pNode()->mAux--;
    }
  
    last = rNetwork().pFirst();
    chain = rNetwork().pFirst()->mpNext;
  
    while (chain) 
    {
      bool    short_curcuit = true;
      NetworkType::NodeType ** curPtr = &chain;
      i = 0;
  
      while (*curPtr) 
      {
        if ((((*curPtr)->mType & NT_MODEL) && !((*curPtr)->mType & NT_TEE))
          || (*curPtr)->mAux == 0) 
        {
          for (j = 0; j < (*curPtr)->mNLinks; j++) 
          {
            (*curPtr)->mpLinks[j].pNode()->mAux--;
          }
  
          last = (last->mpNext = *curPtr);
          last->mAux = i++;
          *curPtr = (*curPtr)->mpNext;
          short_curcuit = false;
        } 
        else 
        {
          curPtr = &(*curPtr)->mpNext;
        }
      }
  
      if (short_curcuit) 
      {
  //      fprintf(stderr, "Nodes in loop: ");
  //      for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->next)
  //        fprintf(stderr, "%d %d", *curPtr - mpNodes, (*curPtr)->mType);
  //      fprintf(stderr, "\n");
        Error("Loop of non-emiting nodes found in network");
      }
    }
  
    last->mpNext = NULL;
  
    /// !!! What is this sorting links good for ???
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (node->mNLinks > 1)
        qsort(node->mpLinks, node->mNLinks, sizeof(NetworkType::LinkType), cmplnk);
    }
  
  // Sort nodes for backward propagation
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext)
      node->mAux = node->mNLinks;
  
    for (i = 0; i < rNetwork().pLast()->mNBackLinks; i++)
      rNetwork().pLast()->mpBackLinks[i].pNode()->mAux--;
  
    last = rNetwork().pLast();
    chain = rNetwork().pLast()->mpBackNext;
    i = 0;
  
    while (chain) 
    {
      bool short_curcuit = true;
      NetworkType::NodeType **curPtr = &chain;
  
      while (*curPtr) 
      {
        if ((((*curPtr)->mType & NT_MODEL) && !((*curPtr)->mType & NT_TEE))
          || (*curPtr)->mAux == 0) 
        {
          for (j = 0; j < (*curPtr)->mNBackLinks; j++) 
          {
            (*curPtr)->mpBackLinks[j].pNode()->mAux--;
          }
  
          last = (last->mpBackNext = *curPtr);
          last->mAux = i++;
          *curPtr = (*curPtr)->mpBackNext;
          short_curcuit = false;
        } 
        else 
        {
          curPtr = &(*curPtr)->mpBackNext;
        }
      }
  
      assert(!short_curcuit); // Shouldn't happen, since it didnot happen before
    }
  
    last->mpBackNext = NULL;
  
    /// !!! What is this sorting links good for ???
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (node->mNBackLinks > 1)
        qsort(node->mpBackLinks, node->mNBackLinks, sizeof(NetworkType::LinkType), cmplnk);
    }
  } // Decoder::SortNodes();
*/  
  

  //***************************************************************************
  //***************************************************************************
  Decoder::NetworkType::NodeType* 
  Decoder:: 
  PhoneNodesToModelNodes(ModelSet * pHmms, ModelSet *pHmmsToUpdate, int& maxStatesInModel)
  {
    NetworkType::iterator p_node;
    NetworkType::iterator p_last_node;
  
    if (pHmmsToUpdate == NULL) 
      pHmmsToUpdate = pHmms;

    mNumberOfNetStates = 0;
    maxStatesInModel   = 0;
  
    // for (p_node =  rNetwork().begin(); 
    // p_node != rNetwork().end(); 
    // p_last_node = p_node, p_node = mCompactRepresentation 
    // ? (p_node->mNLinks ? reinterpret_cast<NetworkType::NodeType*>(reinterpret_cast<NodeBasic<NODE_REGULAR, LINK_REGULAR>* >(p_node)+1) : NULL)
    // : p_node->mpNext)
    
    for (p_node =  rNetwork().begin(); 
         p_node != rNetwork().end(); 
         p_last_node = p_node, p_node++)
    {
      if (p_node->mType & NT_PHONE)  
      {
        Macro *macro;
        p_node->mType &= ~NT_PHONE;
        p_node->mType |= NT_MODEL;        
        macro = FindMacro(&pHmms->mHmmHash, p_node->mpName);
    
        if (macro == NULL) 
        {
          Error("Model %s not defined in %sHMM set", p_node->mpName,
                pHmmsToUpdate != pHmms ? "alignment " : "");
        }
        p_node->mpHmm =  (Hmm *) macro->mpData;
        
        if (!mCompactRepresentation) 
        {
          if (pHmmsToUpdate != pHmms) 
          {
            macro = FindMacro(&pHmmsToUpdate->mHmmHash, p_node->mpName);
            if (macro == NULL) {
              Error("Model %s not defined in HMM set", p_node->mpName,
                    pHmmsToUpdate != pHmms ? "" : "target ");
            }
            p_node->mpHmmToUpdate = (Hmm *) macro->mpData;
          } 
          else
          {
            p_node->mpHmmToUpdate = p_node->mpHmm;
          }
          p_node->mpHmmToUpdate->mpMacro->mOccurances++;
        }
        
        int nstates = p_node->mpHmm->mNStates;
        
        if (p_node->mpHmm->mpTransition->mpMatrixO[nstates - 1] > LOG_MIN) 
        {
          p_node->mType |= NT_TEE;
        }
#ifndef NDEBUG
        p_node->mAux2 = 0;
        p_node->mEmittingStateId = mNumberOfNetStates;
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
  void
  Decoder:: 
  Init(NetworkType::NodeType* pFirstNode, ModelSet* pHmms, ModelSet* pHmmsToUpdate, 
      bool compactRepresentation) 
  {
    int                       maxStatesInModel;
    NetworkType::NodeType*   p_last_node;

    mCompactRepresentation = compactRepresentation;
    rNetwork().SetFirst(pFirstNode);
  
    p_last_node = PhoneNodesToModelNodes(pHmms, pHmmsToUpdate, maxStatesInModel);
    rNetwork().SetLast(p_last_node);
  
    // Need to contain at least 2 states for propper token passing
    if (maxStatesInModel < 2)
    {
      Error("The network does not contain any model with 2 or more states (possibly no dictionary suplied for word-based network)");
    }

    mpAuxTokens = new Token[maxStatesInModel-1];
    mpOutPCache = (Cache*) malloc(pHmms->mNStates      * sizeof(Cache));
    mpMixPCache = (Cache*) malloc(pHmms->mNMixtures    * sizeof(Cache));
    
    if (mpOutPCache == NULL || mpMixPCache == NULL) 
    {
      Error("Insufficient memory");
    }
  
    mWPenalty          = 0.0;
    mMPenalty          = 0.0;
    mPronScale         = 1.0;
    mTranScale         = 1.0;
    mOutpScale         = 1.0;
    mOcpScale          = 1.0;
    mLmScale           = 1.0;
    mLatticeGeneration = false;
    
    OutputProbability =
      pHmms->mOutPdfKind == KID_DiagC     ? &::DiagCGaussianMixtureDensity :
      pHmms->mOutPdfKind == KID_PDFObsVec ? &::FromObservationAtStateId    : NULL;
    
    PassTokenInNetwork  = &PassTokenMax;
    PassTokenInModel    = &PassTokenMax;
    
    mPropagDir          = FORWARD;
    mAlignment          = WORD_ALIGNMENT;
    mpThreshState       = NULL;
    mPruningThresh      = -LOG_0;
    mpModelSet          = pHmms;
    mpModelSetToUpdate  = pHmmsToUpdate;
    mCollectAlphaBeta   = 0;
    mAccumType          = AT_ML;            
    mSearchPaths        = SP_ALL;

    mKeepExitToken      = false;
  }
  //***************************************************************************



  //***************************************************************************
  //***************************************************************************
  void
  Decoder:: 
  Init(ModelSet* pHmms, ModelSet* pHmmsToUpdate, bool compactRepresentation) 
  {
    int maxStatesInModel;
    
    mCompactRepresentation = compactRepresentation;
    PhoneNodesToModelNodes(pHmms, pHmmsToUpdate, maxStatesInModel);
  
    mpAuxTokens = new Token[maxStatesInModel-1];
    mpOutPCache = (Cache*) malloc(pHmms->mNStates      * sizeof(Cache));
    mpMixPCache = (Cache*) malloc(pHmms->mNMixtures    * sizeof(Cache));
    
    if (mpOutPCache == NULL || mpMixPCache == NULL) 
    {
      Error("Insufficient memory");
    }
  
    mWPenalty          = 0.0;
    mMPenalty          = 0.0;
    mPronScale         = 1.0;
    mTranScale         = 1.0;
    mOutpScale         = 1.0;
    mOcpScale          = 1.0;
    mLmScale           = 1.0;
    mLatticeGeneration = false;
    
    OutputProbability =
      pHmms->mOutPdfKind == KID_DiagC     ? &::DiagCGaussianMixtureDensity :
      pHmms->mOutPdfKind == KID_PDFObsVec ? &::FromObservationAtStateId    : NULL;
    
    PassTokenInNetwork  = &PassTokenMax;
    PassTokenInModel    = &PassTokenMax;
    
    mPropagDir          = FORWARD;
    mAlignment          = WORD_ALIGNMENT;
    mpThreshState       = NULL;
    mPruningThresh      = -LOG_0;
    mpModelSet          = pHmms;
    mpModelSetToUpdate  = pHmmsToUpdate;
    mCollectAlphaBeta   = 0;
    mAccumType          = AT_ML;            
    mSearchPaths        = SP_ALL;

    mKeepExitToken      = false;
  }

  
  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  Clear()
  {
    rNetwork().Clear();
    
    delete [] mpAuxTokens;
    free(mpOutPCache);
    free(mpMixPCache);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  Decoder::
  ViterbiInit()
  {
    PassTokenInModel    = &PassTokenMax;
    PassTokenInNetwork  = mLatticeGeneration ? &PassTokenMaxForLattices : &PassTokenMax;
    mPropagDir          = FORWARD;
    mpModelSet->ResetXformInstances();
  
    mTime = 0; // Must not be set to -mpModelSet->mTotalDelay yet
               // otherwise token cannot enter first model node
               // with start set to 0
               
    TokenPropagationInit();
    mTime = -mpModelSet->mTotalDelay;
  }  
  
  
  //***************************************************************************
  //***************************************************************************
  void
  Decoder:: 
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
  FLOAT 
  Decoder::
  ViterbiDone(Label** pLabels, Lattice* pLattice)
  {
    FLOAT tot_like = LOG_0;
    
    if (pLabels  != NULL) *pLabels  = NULL;
    if (pLattice != NULL) pLattice->SetFirst(NULL); 
    
    // Although this peace of code is also in TokenPropagationDone(),
    // we must make sure that no token other then the token in the last node
    // is active before calling pGetLattice(). This is because pGetLattice() rely
    // on mNReferences and expect that (except the last node) WLR is referenced
    // only by other WRLs and not by any token.
    
    // !!! It would be better to rewrite pGetLattice() so that it do not depend on mNReferences
    
    while(mpActiveModels != NULL)
    {
      for (size_t j=0; j < mpActiveModels->mpHmm->mNStates; j++)
        KillToken(&mpActiveModels->mpAnr->mpTokens[j]);
      DeactivateModel(mpActiveModels);
    }

    if (rNetwork().pLast()->mpAnr &&
        rNetwork().pLast()->mpAnr->mpExitToken && 
        rNetwork().pLast()->mpAnr->mpExitToken->IsActive()) 
    {
      tot_like = rNetwork().pLast()->mpAnr->mpExitToken->mLike;
      
      if (pLabels != NULL)
        *pLabels = rNetwork().pLast()->mpAnr->mpExitToken->pGetLabels();
        
      if (mLatticeGeneration && pLattice != NULL)
      {
        //pLattice->SetFirst(rNetwork().pLast()->mpAnr->mpExitToken->pGetLattice());
        pLattice->BuildFromWlr(rNetwork().pLast()->mpAnr->mpExitToken->mpWlr);
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
  /*
  FLOAT 
  Decoder::
  ViterbiDone(Label** pLabels, NetworkType::NodeType** pLattice)
  {
    FLOAT totLike = LOG_0;
    
    if (pLabels  != NULL) *pLabels  = NULL;
    if (pLattice != NULL) *pLattice = NULL; 
    
    // Although this peace of code is also in TokenPropagationDone(),
    // we must make sure that no token other then the token in the last node
    // is active before calling pGetLattice(). This is because pGetLattice() rely
    // on mNReferences and expect that (except the last node) WLR is referenced
    // only by other WRLs and not by any token.
    
    // !!! It would be better to rewrite pGetLattice() so that it do not depend on mNReferences
    
    while(mpActiveModels != NULL)
    {
      for (size_t j=0; j < mpActiveModels->mpHmm->mNStates; j++)
        KillToken(&mpActiveModels->mpAnr->mpTokens[j]);
      DeactivateModel(mpActiveModels);
    }

    if (rNetwork().pLast()->mpAnr &&
        rNetwork().pLast()->mpAnr->mpExitToken && 
        rNetwork().pLast()->mpAnr->mpExitToken->IsActive()) 
    {
      totLike = rNetwork().pLast()->mpAnr->mpExitToken->mLike;
      
      if (pLabels != NULL)
        *pLabels = rNetwork().pLast()->mpAnr->mpExitToken->pGetLabels();
        
      if (mLatticeGeneration && pLattice != NULL)
        *pLattice = rNetwork().pLast()->mpAnr->mpExitToken->pGetLattice();
    }
    
    TokenPropagationDone();
  
  #ifdef DEBUG_MSGS
    printf("Number of WordLinkRecord Unreleased: "); PrintNumOfWLR();
    printf("Number of output prob. computations: %d\n", gaus_computaions);
  #endif
    return totLike;
  }
  */


  
  //***************************************************************************
  //***************************************************************************
  Decoder::FWBWRet
  Decoder::
  ForwardBackward(const Matrix<FLOAT>& rFeatureMatrix, size_t nFrames)
  {
    int         i;
    Cache*      p_out_p_cache;
    FWBWRet     ret;
    ModelSet*   hmms = mpModelSet;
    
    p_out_p_cache = (Cache*) malloc(nFrames * hmms->mNStates * sizeof(Cache));
  
    if (p_out_p_cache == NULL) 
    {
      Error("Insufficient memory");
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
  
    // tady nekde je chyba
    // !!!!
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
  
    if (rNetwork().pLast()->mpAnr == NULL 
    || !rNetwork().pLast()->mpAnr->mpExitToken->IsActive())  
    { // No token survivered
      TokenPropagationDone();
      FreeFWBWRecords();
      mpOutPCache = p_out_p_cache;
      ret.totLike = LOG_0;
      return ret;
    }
  
    ret.totLike = rNetwork().pLast()->mpAnr->mpExitToken->mLike; //  totalLikelihood;
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
  
    if (rNetwork().pFirst()->mpAnr == NULL || !rNetwork().pFirst()->mpAnr->mpExitToken->IsActive()) 
    { // No token survivered
      TokenPropagationDone();
      FreeFWBWRecords();
      ret.totLike = LOG_0;
      return ret;
    }
  
    ret.totLike = HIGHER_OF(ret.totLike, rNetwork().pFirst()->mpAnr->mpExitToken->mLike); //  totalLikelihood;
    // Backward pass P can differ from forward pass P because of the precision
    // problems. Take the higher one to decrease the possibility of getting
    // an occupation probability (when normalizing by P) higher that one.
  
    FloatInLog fil_ret_totLike = {ret.totLike, 0};
    ret.avgAccuracy  = FIL_Div(rNetwork().pFirst()->mpAnr->mpExitToken->mAccuracy, fil_ret_totLike);
    TokenPropagationDone();
  
    // There may be remaining records in mpAlphaBetaListReverse unused in
    // backward pass. Free them and set mpAlphaBetaListReverse's to NULLs;
    NetworkType::NodeType* node;
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (!(node->mType & NT_MODEL)) 
        continue;
  
      while (node->mpAlphaBetaListReverse) 
      {
        FWBWR* fwbwr = node->mpAlphaBetaListReverse;
        node->mpAlphaBetaListReverse = fwbwr->mpNext;
        free(fwbwr);
      }
    }
  
    return ret;
  }
  // ForwardBackward(const Matrix<FLOAT>& rFeatureMatrix, size_t nFrames)
  //***************************************************************************

  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  Decoder::
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
    NetworkType::NodeType*             node;
  
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
  
    TP = rNetwork().pLast()->mpAnr == NULL ? LOG_0 : rNetwork().pLast()->mpAnr->mpExitToken->mLike;
    
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
      F = exp(-sigSlope * F);
      F = (sigSlope*F) / SQR(1+F);
      printf("weight: %g\n", F);
      weight *= F;
    }

    if (P < LOG_MIN) return LOG_0;
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) {
      if (node->mType & NT_MODEL && node->mpAlphaBetaList != NULL &&
        node->mpAlphaBetaList->mTime == 0) {
        node->mpAlphaBetaListReverse = node->mpAlphaBetaList;
        node->mpAlphaBetaList = node->mpAlphaBetaList->mpNext;
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
  
      for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
      { //for every model
  //    for (k=0; k < nnodes; k++) {
  //      NetworkType::NodeType *node = &mpNodes[k];
        if (node->mType & NT_MODEL &&
          node->mpAlphaBetaList != NULL &&
          node->mpAlphaBetaList->mTime == mTime+1) {
  
          struct AlphaBetaMPE *st;
          int Nq       = node->mpHmm->mNStates;
          st = node->mpAlphaBetaList->mpState;
  
          for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
            if (st[j].mAlpha + st[j].mBeta - P > MIN_LOG_WEGIHT) {
              assert(node->mpAlphaBetaListReverse->mTime == mTime);
  
              ReestState(
                  node,
                  j-1,
                  (st[j].mAlpha + st[j].mBeta - P)  * mOcpScale,
                  -weight,
                  obs,
                  obs2);
            }
          }
          
          if (node->mpAlphaBetaListReverse) 
            free(node->mpAlphaBetaListReverse);
          
          node->mpAlphaBetaListReverse = node->mpAlphaBetaList;
          node->mpAlphaBetaList = node->mpAlphaBetaList->mpNext;
        }
      }
    }
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) {
      if (node->mpAlphaBetaListReverse != NULL)
        free(node->mpAlphaBetaListReverse);
    }
  
  
    ////////////////// Numerator accumulation //////////////////
    mSearchPaths = SP_TRUE_ONLY;
  
  
    ForwardBackward(rObsMx, nFrames);
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) {
      if (node->mType & NT_MODEL && node->mpAlphaBetaList != NULL &&
        node->mpAlphaBetaList->mTime == 0) {
        node->mpAlphaBetaListReverse = node->mpAlphaBetaList;
        node->mpAlphaBetaList = node->mpAlphaBetaList->mpNext;
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
  
      for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
      { //for every model
        if (node->mType & NT_MODEL 
        && node->mpAlphaBetaList != NULL 
        && node->mpAlphaBetaList->mTime == mTime+1) 
        {
  
          struct AlphaBetaMPE *st;
          int Nq       = node->mpHmm->mNStates;
          FLOAT *aq    = node->mpHmm->        mpTransition->mpMatrixO;
          FLOAT *aqacc = node->mpHmmToUpdate->mpTransition->mpMatrixO + SQR(Nq);
  //        int qt_1 = (mNumberOfNetStates * mTime) + node->mEmittingStateId;
  //        int qt = qt_1 + mNumberOfNetStates;
  
          st = node->mpAlphaBetaList->mpState;
  
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
              FLOAT bjtO =mpOutPCache[p_hmms_alig->mNStates * mTime +
                                        node->mpHmm->mpState[j-1]->mID].mValue;
              // ForwardBackward() set mpOutPCache to contain out prob. for all frames
  
              assert(node->mpAlphaBetaListReverse->mTime == mTime);
  
  //            if (!mmi_den_pass) {
              for (i = 0; i < Nq - 1; i++) 
              {
                LOG_INC(aqacc[i * Nq + j],
                        aq[i * Nq + j]    * mTranScale +
                        (node->mpAlphaBetaListReverse->mpState[i].mAlpha +
                        bjtO              * mOutpScale +
                        st[j].mBeta - TP)   * mOcpScale);
              }
  //            }
  
              if (origAccumType == AT_MMI)
              {
                ReestState(
                    node,
                    j-1,
                    (st[j].mAlpha + st[j].mBeta - TP)  * mOcpScale,
                    weight, 
                    obs, 
                    obs2);
              }
              else // origAccumType == AT_MCE
              {
                ReestState(node, j-1,
                    (st[j].mAlpha + st[j].mBeta - TP + LogAdd(TP,P) - P)  * mOcpScale,
                    weight, 
                    obs, 
                    obs2);
              }
            }
          }
          
          if (node->mpAlphaBetaListReverse) 
            free(node->mpAlphaBetaListReverse);
            
          node->mpAlphaBetaListReverse = node->mpAlphaBetaList;
          node->mpAlphaBetaList = node->mpAlphaBetaList->mpNext;
        }
      }
    }
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (node->mpAlphaBetaListReverse != NULL)
        free(node->mpAlphaBetaListReverse);
    }
  
    mAccumType = origAccumType;
    return TP;
  } 
  // MCEReest
  //***************************************************************************
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  Decoder::
  BaumWelchReest(const Matrix<FLOAT>& rObsMx, const Matrix<FLOAT>& rObsMx2, int nFrames, FLOAT weight)
  {
    struct FWBWRet          fwbw;
    FLOAT                   P;
    FLOAT                   update_dir;
    int                     i;
    int                     j;
    int                     k;
    ModelSet*               p_hmms_alig = mpModelSet;
    ModelSet*               p_hmms_upd = mpModelSetToUpdate;
    NetworkType::NodeType*                   node;
  
    fwbw = ForwardBackward(rObsMx, nFrames);
    P    = fwbw.totLike;
    
    if (P < LOG_MIN) 
      return LOG_0;
  
#ifdef MOTIF
    FLOAT* ocprob = (FLOAT*) malloc(mNumberOfNetStates * (nFrames+1) * sizeof(FLOAT));
    for (i=0; i<mNumberOfNetStates * (nFrames+1); i++) ocprob[i] = 0;
#endif
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (node->mType & NT_MODEL && node->mpAlphaBetaList != NULL &&
        node->mpAlphaBetaList->mTime == 0) 
      {
        node->mpAlphaBetaListReverse = node->mpAlphaBetaList;
        node->mpAlphaBetaList = node->mpAlphaBetaList->mpNext;
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
    
    if (mpMixPCache == NULL)
    { 
      Error("Insufficient memory");
    }
    
    for (i = 0; i < k; i++)
    { 
      mpMixPCache[i].mTime = UNDEF_TIME;
    }
    
  // Update accumulators
    for (mTime = 0; mTime < nFrames; mTime++) 
    { //for every frame
      //FLOAT* obs  = pObsMx +p_hmms_alig->mInputVectorSize*(mTime+p_hmms_alig->mTotalDelay);
      //FLOAT* obs2 = pObsMx2+p_hmms_upd->mInputVectorSize*(mTime+p_hmms_upd->mTotalDelay);
      
      FLOAT* obs  = rObsMx [mTime + p_hmms_alig->mTotalDelay];
      FLOAT* obs2 = rObsMx2[mTime + p_hmms_upd->mTotalDelay];
      
      p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
      
      if (p_hmms_alig != p_hmms_upd)
        p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
  
      for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
      { //for every model
  //    for (k=0; k < nnodes; k++) {
  //      NetworkType::NodeType *node = &mpNodes[k];
        if (node->mType & NT_MODEL &&
          node->mpAlphaBetaList != NULL &&
          node->mpAlphaBetaList->mTime == mTime+1) 
        {
          struct AlphaBetaMPE *st;
          int Nq       = node->mpHmm->mNStates;
          FLOAT* aq    = node->mpHmm->        mpTransition->mpMatrixO;
          FLOAT* aqacc = node->mpHmmToUpdate->mpTransition->mpMatrixO + SQR(Nq);
//        int qt_1 = (mNumberOfNetStates * mTime) + node->mEmittingStateId;
//        int qt = qt_1 + mNumberOfNetStates;
  
          st = node->mpAlphaBetaList->mpState;
  
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
            ocprob[mNumberOfNetStates * (mTime+1) + node->mEmittingStateId + j]
//            = node->mPhoneAccuracy;
              = exp(st[j].mAlpha+st[j].mBeta-P) *
                ((1-2*static_cast<int>(st[j].mAlphaAccuracy.negative)) * exp(st[j].mAlphaAccuracy.logvalue - st[j].mAlpha) +
                (1-2*static_cast<int>(st[j].mBetaAccuracy.negative))  * exp(st[j].mBetaAccuracy.logvalue  - st[j].mBeta)
                - (1-2*static_cast<int>(fwbw.avgAccuracy.negative)) * exp(fwbw.avgAccuracy.logvalue)
                );
  
#endif
//            int qt_1   = qt - mNumberOfNetStates;
              FLOAT bjtO =mpOutPCache[p_hmms_alig->mNStates * mTime +
                                        node->mpHmm->mpState[j-1]->mID].mValue;
              // ForwardBackward() set mpOutPCache to contain out prob. for all frames
  
              assert(node->mpAlphaBetaListReverse->mTime == mTime);
  
//            if (!mmi_den_pass) {
              for (i = 0; i < Nq - 1; i++) 
              {
                LOG_INC(aqacc[i * Nq + j],
                        aq[i * Nq + j] * mTranScale +
                        (node->mpAlphaBetaListReverse->mpState[i].mAlpha +
                        bjtO * mOutpScale + st[j].mBeta - P) * mOcpScale);
              }
//            }
  
              if (mAccumType == AT_MFE || mAccumType == AT_MPE) 
              {
                update_dir = (1-2*static_cast<int>(st[j].mAlphaAccuracy.negative)) * exp(st[j].mAlphaAccuracy.logvalue - st[j].mAlpha) +
                             (1-2*static_cast<int>(st[j].mBetaAccuracy.negative))  * exp(st[j].mBetaAccuracy.logvalue  - st[j].mBeta)  -
                             (1-2*static_cast<int>(fwbw.avgAccuracy.negative))     * exp(fwbw.avgAccuracy.logvalue);
              } 
              else 
              {
                update_dir = 1.0;
              }
  
              ReestState(
                node, 
                j - 1,
                (st[j].mAlpha + st[j].mBeta - P)  * mOcpScale,
                update_dir*weight, 
                obs, 
                obs2);
            }
          }
          
          if (node->mpAlphaBetaListReverse) 
            free(node->mpAlphaBetaListReverse);
            
          node->mpAlphaBetaListReverse = node->mpAlphaBetaList;
          node->mpAlphaBetaList = node->mpAlphaBetaList->mpNext;
        }
      }
    }
  
    for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
    {
      if (node->mpAlphaBetaListReverse != NULL)
        free(node->mpAlphaBetaListReverse);
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
  FLOAT 
  Decoder::
  ViterbiReest(
    const Matrix<FLOAT>&  rObsMx, 
    const Matrix<FLOAT>&  rObsMx2, 
    int                   nFrames, 
    FLOAT                 weight)
  {
    int                     t;
    WordLinkRecord *        wlr;
    NetworkType::NodeType* prevnode =        NULL;
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
  
    if (rNetwork().pLast()->mpAnr == NULL || !rNetwork().pLast()->mpAnr->mpExitToken->IsActive()) 
    {
      ViterbiDone(NULL);
      return LOG_0;
    }
  
    p_hmms_alig->ResetXformInstances();
  
    if (p_hmms_alig != p_hmms_upd)
      p_hmms_upd->ResetXformInstances();
  
    // invert order of WRLs
    wlr = rNetwork().pLast()->mpAnr->mpExitToken->mpWlr;
    
    while (wlr->mpNext != NULL) 
    {
      WordLinkRecord *tmp = wlr->mpNext->mpNext;
      wlr->mpNext->mpNext = rNetwork().pLast()->mpAnr->mpExitToken->mpWlr;
      rNetwork().pLast()->mpAnr->mpExitToken->mpWlr = wlr->mpNext;
      wlr->mpNext = tmp;
    }
  
    for (mTime = -p_hmms_alig->mTotalDelay; mTime < 0; mTime++) 
    {
      //FLOAT *obs = rObsMx+p_hmms_alig->mInputVectorSize*(mTime+p_hmms_alig->mTotalDelay);
      FLOAT *obs = rObsMx[mTime+p_hmms_alig->mTotalDelay];
      p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
    }
  
    if (p_hmms_alig != p_hmms_upd) 
    {
      //FLOAT *obs2 = pObsMx2+p_hmms_upd->mInputVectorSize*(mTime+p_hmms_upd->mTotalDelay);
      FLOAT *obs2 = rObsMx2[mTime+p_hmms_upd->mTotalDelay];
      
      for (mTime = -p_hmms_upd->mTotalDelay; mTime < 0; mTime++) 
      {
        p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
      }
    }
  
  // Update accumulators
    for (wlr = rNetwork().pLast()->mpAnr->mpExitToken->mpWlr; wlr != NULL; wlr = wlr->mpNext) 
    {
      NetworkType::NodeType *node   = wlr->mpNode;
      int Nq       = node->mpHmmToUpdate->mNStates;
      FLOAT *aqacc = node->mpHmmToUpdate->mpTransition->mpMatrixO + SQR(Nq);
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
      LOG_INC(aqacc[currstate * Nq + currstate], log(duration-1));
      LOG_INC(aqacc[currstate * Nq + nextstate], 0 /*ln(1)*/);
  //    }
  
      for (; mTime < wlr->mTime; mTime++) 
      {
        //for every frame of segment
        //FLOAT *obs  = pObsMx  + p_hmms_alig->mInputVectorSize*(mTime+p_hmms_alig->mTotalDelay);
        //FLOAT *obs2 = pObsMx2 + p_hmms_upd->mInputVectorSize*(mTime+p_hmms_upd->mTotalDelay);
        FLOAT *obs  = rObsMx [mTime+p_hmms_alig->mTotalDelay];
        FLOAT *obs2 = rObsMx2[mTime+p_hmms_upd->mTotalDelay];
        
  
        p_hmms_alig->UpdateStacks(obs, mTime, FORWARD);
        
        if (p_hmms_alig != p_hmms_upd)
          p_hmms_upd->UpdateStacks(obs2, mTime, FORWARD);
        
        ReestState(node, currstate-1, 0.0, 1.0*weight, obs, obs2);
      }
    }
  
    P = rNetwork().pLast()->mpAnr->mpExitToken->mpWlr->mLike;
    //ViterbiDone(net, NULL);
    ViterbiDone(NULL);
    return P;
  }
    
    
    
  //***************************************************************************
  //***************************************************************************
  // Token section
  //***************************************************************************
  //***************************************************************************
  void 
  Token::
  AddWordLinkRecord(Node<NodeBasicContent, LinkContent, NODE_REGULAR, LINK_REGULAR>* pNode, int stateIndex, int time)
  {
    WordLinkRecord* wlr;
  
    if ((wlr = (WordLinkRecord*) malloc(sizeof(WordLinkRecord))) == NULL)
      Error("Insufficient memory");
  
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
    
# ifdef bordel_staff
    if (!mpTWlr) 
      mpTWlr = wlr;
# endif // bordel_staff
  
# ifdef DEBUG_MSGS
    wlr->mpTmpNext = firstWLR;
    firstWLR = wlr;
    wlr->mIsFreed = false;
# endif
  }
  //***************************************************************************

  
  //***************************************************************************
  //***************************************************************************
  static void 
  MakeLatticeNodesForWordLinkRecords(WordLinkRecord* pWlr, 
      Lattice::NodeType*& rpFirst)
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
                         
    // If rpFirst == NULL, pWlr is the wery last WLR referenced by token in last state
    pNode->mNLinks     = rpFirst == NULL ? 0 : pWlr->mNReferences;


    // Allocate space for links and backlinks      
    pNode->mpLinks     = (Lattice::LinkType *) malloc(pNode->mNLinks     * sizeof(Lattice::LinkType));
    pNode->mpBackLinks = (Lattice::LinkType *) malloc(pNode->mNBackLinks * sizeof(Lattice::LinkType));
    if (pNode->mpLinks == NULL || pNode->mpBackLinks == NULL) 
      Error("Insufficient memory");
        
    // Push new node to lattice node list
    pNode->mpNext = rpFirst;
    pNode->mpBackNext = NULL;
    rpFirst = pNode;
    
    if (pNode->mpNext != NULL) 
      pNode->mpNext->mpBackNext = pNode;
      
      
    // Set pWlr->pNode() to the new lattice node, so we can latter easily
    // establish  links between the nodes in the lattice nodes
    pWlr->mpNode = pNode;
      
    // All successors has been already processed, so continue recursively with Wlr's predecessors
    if(pWlr->mpNext != NULL)
    {
      MakeLatticeNodesForWordLinkRecords(pWlr->mpNext, rpFirst);
        
      if(pWlr->mpAltHyps != NULL)
      {
        for(AltHypList::iterator i = pWlr->mpAltHyps->begin(); i != pWlr->mpAltHyps->end(); i++)
          MakeLatticeNodesForWordLinkRecords(i->mpWlr, rpFirst);
      }
    }
  }
  // MakeLatticeNodesForWordLinkRecords(WordLinkRecord* pWlr, Node*& rpFirst)
  //***************************************************************************

  

  //***************************************************************************
  //***************************************************************************
  static void 
  EstablishLinksBetweenLatticeNodes(WordLinkRecord* pWlr)
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
      EstablishLinksBetweenLatticeNodes(pWlr->mpNext);
        
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
          EstablishLinksBetweenLatticeNodes(i->mpWlr);
        }
      }
    }
  }
  // EstablishLinksBetweenLatticeNodes(WordLinkRecord* pWlr)
  //***************************************************************************
  

  //***************************************************************************
  //***************************************************************************
  Token::NodeType*
  Token::
  pGetLattice()
  {
    NodeType* pFirst = NULL;
    MakeLatticeNodesForWordLinkRecords(mpWlr, pFirst);
    EstablishLinksBetweenLatticeNodes(mpWlr);
    
    return pFirst;
  }

  
  //***************************************************************************
  //***************************************************************************
  Label *
  Token::
  pGetLabels()
  {
    WordLinkRecord *  wlr;
    Label *           tmp;
    Label *           level[3] = {NULL, NULL, NULL};
    int               li = 0;
      
    if (!this->IsActive()) 
      return NULL;
    
    for (wlr = this->mpWlr; wlr != NULL; wlr = wlr->mpNext) 
    {
      if ((tmp = (Label *) malloc(sizeof(Label))) == NULL)
        Error("Insufficient memory");
  
      tmp->mScore = wlr->mLike;
      tmp->mStop  = wlr->mTime;
      tmp->mId    = wlr->mStateIdx;
  
      if (wlr->pNode()->mType & NT_MODEL) 
      {
        li = wlr->mStateIdx >= 0 ? 0 : 1;
        tmp->mpData = wlr->pNode()->mpHmm;
        tmp->mpName = wlr->pNode()->mpHmm->mpMacro->mpName;
        tmp->mpNextLevel = level[li+1] ? level[li+1] : level[2];
      } 
      else //if (wlr->pNode()->mpPronun->outSymbol) 
      {
        li = 2;
        tmp->mpData = wlr->pNode()->mpPronun ? wlr->pNode()->mpPronun->mpWord    : NULL;
        tmp->mpName = wlr->pNode()->mpPronun ? wlr->pNode()->mpPronun->outSymbol : NULL;
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
  
        level[li]->mScore -= tmp->mScore;
        level[li]->mStart  = tmp->mStop;
      }
  
      tmp->mpNext = level[li];
      level[li] = tmp;
    }
    
    for (li = 0; li < 3; li++) 
    {
      if (level[li]) 
        level[li]->mStart = 0;
    }
    
    return level[0] ? level[0] : level[1] ? level[1] : level[2];
  }

}; // namespace STK
