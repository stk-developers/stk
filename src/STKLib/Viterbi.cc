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

#include "Viterbi.h"
#include "labels.h"
#include "common.h"
#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
#include <cstring>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include <iostream>
using namespace std;

#ifdef MOTIF
#include "imagesc.h"
#endif

//#define TRACE_TOKENS
#define SQR(x) ((x) * (x))
//#define IS_ACTIVE(token) ((token).mLike > LOG_MIN)
#define FORWARD_PASS (net->mPropagDir == FORWARD)

namespace STK
{
  void AddWordLinkRecord(Token *token, Node *node, int state_idx, int time);
  void TokenPropagationInNetwork(Network *net);
  void ReestState(Network *net, Node *node, int state_idx,
                  FLOAT logPriorProb, FLOAT updateDir, FLOAT *obs, FLOAT *obs2);
  int  BackwardPruning(int time, Node *node, int state);
  
  #ifdef DEBUG_MSGS
  WLR *firstWLR;
  #endif
  
  void SortNodes(Network *net);
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  PhoneNodesToModelNodes(Node *first, ModelSet *hmms, ModelSet *hmmsToUpdate)
  {
    Node *node;
  
    if (hmmsToUpdate == NULL) hmmsToUpdate = hmms;
  
    for (node = first; node != NULL; node = node->mpNext) {
      if (node->mType & NT_Phone) {
        Macro *macro;
  
        node->mType &= ~NT_Phone;
        node->mType |= NT_Model;
        macro = FindMacro(&hmms->mHmmHash, node->mpName);
        if (macro == NULL) {
          Error("Model %s not defined in %sHMM set", node->mpName,
                hmmsToUpdate != hmms ? "alignment " : "");
        }
        node->mpHmm = node->mpHmmToUpdate = (Hmm *) macro->mpData;
  
        if (hmmsToUpdate != hmms) {
          macro = FindMacro(&hmmsToUpdate->mHmmHash, node->mpName);
          if (macro == NULL) {
            Error("Model %s not defined in HMM set", node->mpName,
                  hmmsToUpdate != hmms ? "" : "target ");
          }
          node->mpHmmToUpdate = (Hmm *) macro->mpData;
        }
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  InitNetwork(Network *net, Node *first, ModelSet *hmms, ModelSet *hmmsToUpdate) {
    Node *node;
    int maxStatesInModel = 0;
    int i;
  
    net->mpFirst = net->mpLast  = first;
  
    PhoneNodesToModelNodes(first, hmms, hmmsToUpdate);
  
    //Allocate tokens and count emiting states
    net->mNumberOfNetStates = 0;
    
    for (node = first; node != NULL; net->mpLast = node, node = node->mpNext) 
    {
      
  #ifndef NDEBUG
      node->aux2 = 0;
  #endif
      int numOfTokens = 1;
      if (node->mType & NT_Model) 
      {
        numOfTokens = node->mpHmm->mNStates;
        node->mpHmmToUpdate->mpMacro->mOccurances++;
        if (node->mpHmm->mpTransition->mpMatrixO[numOfTokens - 1] > LOG_MIN) {
          node->mType |= NT_Tee;
        }
      } 
      else if (node->mType & NT) 
      {
        numOfTokens = 1;
      } 
      else 
      {
        Error("Fatal: Incorect node type");
      }
  
      node->tokens = (Token *) malloc(numOfTokens * sizeof(Token));
      if (node->tokens == NULL) 
        Error("Insufficient memory");
  
      node->exitToken = &node->tokens[numOfTokens-1];
      node->estate_id = net->mNumberOfNetStates;
      
      if (node->mType & NT_Model) {
        int nstates = node->mpHmm->mNStates;
  
        if (maxStatesInModel < nstates) maxStatesInModel = nstates;
        assert(nstates >= 2); // two non-emiting states
        net->mNumberOfNetStates += nstates;
      }
    }
  
    SortNodes(net);
  
    
    net->mpAuxTokens = (Token *) malloc((maxStatesInModel-1) * sizeof(Token));
    net->mpOutPCache = (Cache *) malloc(hmms->mNStates   * sizeof(Cache));
    net->mpMixPCache = (Cache *) malloc(hmms->mNMixtures * sizeof(Cache));
  
    if (net->mpAuxTokens == NULL ||
      net->mpOutPCache == NULL || net->mpMixPCache == NULL) {
      Error("Insufficient memory");
    }
  
    for (i = 0; i < maxStatesInModel-1; i++) {
      net->mpAuxTokens[i].mLike = LOG_0;
      net->mpAuxTokens[i].wlr = NULL;
    }
  
    net->mWPenalty          = 0.0;
    net->mMPenalty          = 0.0;
    net->mPronScale         = 1.0;
    net->mTranScale         = 1.0;
    net->mOutpScale         = 1.0;
    net->mOcpScale          = 1.0;
    net->mLmScale           = 1.0;
    net->OutputProbability =
      hmms->mOutPdfKind == KID_DiagC     ? &DiagCGaussianMixtureDensity :
      hmms->mOutPdfKind == KID_PDFObsVec ? &FromObservationAtStateId    : NULL;
  
    net->PassTokenInNetwork= &PassTokenMax;
    net->PassTokenInModel  = &PassTokenMax;
    net->mPropagDir         = FORWARD;
    net->mAlignment         = WORD_ALIGNMENT;
    net->mpThreshState       = NULL;
    net->mPruningThresh     = -LOG_0;
    net->mpModelSet            = hmms;
    net->mpModelSetToUpdate    = hmmsToUpdate;
    net->mCollectAlphaBeta  = 0;
  //  net->mmi_den_pass      = 0;
    net->mAccumType          = AT_ML;            
    net->mSearchPaths        = SP_ALL;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ReleaseNetwork(Network *net)
  {
    Node *node;
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) free(node->tokens);
    FreeNetwork(net->mpFirst);
    free(net->mpAuxTokens);
    free(net->mpOutPCache);
    free(net->mpMixPCache);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  cmplnk(const void *a, const void *b)
  {
    return ((Link *) a)->mpNode->mAux - ((Link *) b)->mpNode->mAux;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  SortNodes(Network *net)
  {
    int i, j;
    Node *chain, *last, *node;
  
  // Sort nodes for forward (Viterbi) propagation
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      node->mAux = node->mNBackLinks;
    }
  
    for (i = 0; i < net->mpFirst->mNLinks; i++) {
      net->mpFirst->mpLinks[i].mpNode->mAux--;
    }
  
    last = net->mpFirst;
    chain = net->mpFirst->mpNext;
  
    while (chain) {
      BOOL short_curcuit = TRUE;
      Node **curPtr = &chain;
      i = 0;
  
      while (*curPtr) {
        if ((((*curPtr)->mType & NT_Model) && !((*curPtr)->mType & NT_Tee))
          || (*curPtr)->mAux == 0) {
          for (j = 0; j < (*curPtr)->mNLinks; j++) {
            (*curPtr)->mpLinks[j].mpNode->mAux--;
          }
  
          last = (last->mpNext = *curPtr);
          last->mAux = i++;
          *curPtr = (*curPtr)->mpNext;
          short_curcuit = FALSE;
        } else {
          curPtr = &(*curPtr)->mpNext;
        }
      }
  
      if (short_curcuit) {
  //      fprintf(stderr, "Nodes in loop: ");
  //      for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->next)
  //        fprintf(stderr, "%d %d", *curPtr - net->mpNodes, (*curPtr)->mType);
  //      fprintf(stderr, "\n");
        Error("Loop of non-emiting nodes found in network");
      }
    }
  
    last->mpNext = NULL;
  
    /// !!! What is this sorting links good for ???
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
        if (node->mNLinks > 1)
        qsort(node->mpLinks, node->mNLinks, sizeof(Link), cmplnk);
    }
  
  // Sort nodes for backward propagation
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      node->mAux = node->mNLinks;
    }
  
    for (i = 0; i < net->mpLast->mNBackLinks; i++) {
      net->mpLast->mpBackLinks[i].mpNode->mAux--;
    }
  
    last = net->mpLast;
    chain = net->mpLast->mpBackNext;
    i = 0;
  
    while (chain) {
      BOOL short_curcuit = TRUE;
      Node **curPtr = &chain;
  
      while (*curPtr) {
        if ((((*curPtr)->mType & NT_Model) && !((*curPtr)->mType & NT_Tee))
          || (*curPtr)->mAux == 0) {
          for (j = 0; j < (*curPtr)->mNBackLinks; j++) {
            (*curPtr)->mpBackLinks[j].mpNode->mAux--;
          }
  
          last = (last->mpBackNext = *curPtr);
          last->mAux = i++;
          *curPtr = (*curPtr)->mpBackNext;
          short_curcuit = FALSE;
        } else {
          curPtr = &(*curPtr)->mpBackNext;
        }
      }
  
      /*if (short_curcuit) {
        fprintf(stderr, "Nodes in loop: ");
        for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->mpBackNext)
          fprintf(stderr, "%d ", *curPtr - net->mpNodes);
        fprintf(stderr, "\n");
        Error("Loop of non-emiting nodes found in net");
      }*/
  
      assert(!short_curcuit); // Shouldn't happen, since it didnot happen before
    }
  
    last->mpBackNext = NULL;
  
    /// !!! What is this sorting links good for ???
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
        if (node->mNBackLinks > 1)
        qsort(node->mpBackLinks, node->mNBackLinks, sizeof(Link), cmplnk);
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  WriteAlpha(int time, Node *node, int state, Token *token)
  {
    if (node->alphaBetaListReverse == NULL ||
      node->alphaBetaListReverse->mTime != time) 
    {
      size_t    i;
      FWBWR *   newrec;
      
      newrec  = (FWBWR*) malloc(sizeof(FWBWR) +
                                sizeof(newrec->mpState[0]) * (node->mpHmm->mNStates-1));
      if (newrec == NULL) 
        Error("Insufficient memory");
        
      newrec->mpNext = node->alphaBetaListReverse;
      newrec->mTime = time;
      
      for (i=0; i<node->mpHmm->mNStates; i++) 
      {
        newrec->mpState[i].alpha = newrec->mpState[i].beta = LOG_0;
      }
      
      node->alphaBetaListReverse = newrec;
    }
    node->alphaBetaListReverse->mpState[state].alpha = token->mLike;
    node->alphaBetaListReverse->mpState[state].alphaAccuracy = token->mAccuracy;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  WriteBeta(int time, Node *node, int state, Token *token)
  {
  
    // Record for current time must be already moved to
    // alphaBetaList by function BackwardPruning
    assert(node->alphaBetaListReverse == NULL ||
          node->alphaBetaListReverse->mTime < time);
  
    if (node->alphaBetaList != NULL && node->alphaBetaList->mTime == time) {
      node->alphaBetaList->mpState[state].beta = token->mLike;
      node->alphaBetaList->mpState[state].betaAccuracy = token->mAccuracy;
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  BackwardPruning(int time, Node *node, int state)
  {
    while (node->alphaBetaListReverse != NULL &&
          node->alphaBetaListReverse->mTime > time) {
      FWBWR *fwbwr = node->alphaBetaListReverse;
      node->alphaBetaListReverse = fwbwr->mpNext;
      free(fwbwr);
    }
  
    if (node->alphaBetaListReverse != NULL &&
      node->alphaBetaListReverse->mTime == time) {
      FWBWR *fwbwr = node->alphaBetaListReverse;
      node->alphaBetaListReverse = fwbwr->mpNext;
      fwbwr->mpNext = node->alphaBetaList;
      node->alphaBetaList = fwbwr;
    }
  
    return !(node->alphaBetaList != NULL &&
            node->alphaBetaList->mTime == time &&
            node->alphaBetaList->mpState[state].alpha > LOG_MIN);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  FreeFWBWRecords(Network *net)
  {
    Node *node;
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (!(node->mType & NT_Model)) continue;
  
      while (node->alphaBetaList) {
        FWBWR *fwbwr = node->alphaBetaList;
        node->alphaBetaList = fwbwr->mpNext;
        free(fwbwr);
      }
  
      while (node->alphaBetaListReverse) {
        FWBWR *fwbwr = node->alphaBetaListReverse;
        node->alphaBetaListReverse = fwbwr->mpNext;
        free(fwbwr);
      }
    }
  }
  
  
  #ifndef NDEBUG
  int test_for_cycle = 0;
  int HasCycleCounter = 100000;
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  HasCycle(Network *net) 
  {
    Node *node;
    HasCycleCounter++;
    if (!test_for_cycle) return 0;
    for (node = net->mpActiveNodes; node; node= node->nextActiveNode) {
      int i, nlinks = FORWARD_PASS ? node->mNLinks : node->mNBackLinks;
      Link *links   = FORWARD_PASS ? node->mpLinks  : node->mpBackLinks;
      if (node->aux2 == HasCycleCounter) {
        printf("Cycle in list of active nodes\n");
        return 1;
      }
      node->aux2 = HasCycleCounter;
  
      for (i=0; i <nlinks; i++)
        if (links[i].mpNode->aux2 == HasCycleCounter &&
          (!(links[i].mpNode->mType & NT_Model)
          || links[i].mpNode->mType & NT_Tee)) {
          printf("Active node %d listed after his non-model succesor %d\n",
                node - net->mpFirst, links[i].mpNode - net->mpFirst);
  
          return 2;
        }
    }
    return 0;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  AllWordSuccessorsAreActive(Network *net) 
  {
    Node *node;
    if (!test_for_cycle) return 1;
  
    for (node = net->mpActiveNodes; node; node= node->nextActiveNode) {
      int i, nlinks = FORWARD_PASS ? node->mNLinks : node->mNBackLinks;
      Link *links   = FORWARD_PASS ? node->mpLinks  : node->mpBackLinks;
  
      for (i=0; i <nlinks; i++)
  
        if (links[i].mpNode->aux2 != HasCycleCounter &&
          links[i].mpNode != (FORWARD_PASS ? net->mpLast : net->mpFirst) &&
          (!(links[i].mpNode->mType & NT_Model)
          || links[i].mpNode->mType & NT_Tee)) {
          printf("Active node %d has nonactive non-model succesor %d\n",
                node - net->mpFirst, links[i].mpNode - net->mpFirst);
          return 0;
        }
    }
    return 1;
  }
  #endif
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  MarkWordNodesLeadingFrom(Network *net, Node *node)
  {
    int i, nlinks = FORWARD_PASS ? node->mNLinks : node->mNBackLinks;
    Link *links   = FORWARD_PASS ? node->mpLinks  : node->mpBackLinks;
  
    for (i = 0; i < nlinks; i++) {
      Node *lnode = links[i].mpNode;
      if ((lnode->mType & NT_Model && !(lnode->mType & NT_Tee))
        || (lnode == (FORWARD_PASS ? net->mpLast : net->mpFirst))) continue;
      if (lnode->isActiveNode > 0) continue;
  
      if (lnode->isActive) {
        assert(lnode->mType & NT_Tee);
        continue;
      }
      
      if(lnode->isActiveNode-- == 0) {
        lnode->mAux = 0;
        MarkWordNodesLeadingFrom(net, lnode);
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Node *
  ActivateWordNodesLeadingFrom(Network *net, Node *node)
  {
    int i, nlinks = FORWARD_PASS ? node->mNLinks : node->mNBackLinks;
    Link *links   = FORWARD_PASS ? node->mpLinks  : node->mpBackLinks;
  
    for (i = 0; i < nlinks; i++) {
      Node *lnode = links[i].mpNode;
      if ((lnode->mType & NT_Model && !(lnode->mType & NT_Tee))
        || (lnode == (FORWARD_PASS ? net->mpLast : net->mpFirst))) continue;
      if (lnode->isActiveNode++ > 0) continue;
  
      if (lnode->isActive) {
        assert(lnode->mType & NT_Tee);
        continue;
      }
  
      lnode->mAux++;
      if (lnode->isActiveNode < 0) continue;
  
      assert(lnode->isActiveNode == 0);
      lnode->isActiveNode = lnode->mAux;
  
      lnode->nextActiveNode = node->nextActiveNode;
      lnode->prevActiveNode = node;
      if (node->nextActiveNode) node->nextActiveNode->prevActiveNode = lnode;
      node->nextActiveNode  = lnode;
      node = ActivateWordNodesLeadingFrom(net, lnode);
    }
    assert(!HasCycle(net));
    return node;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ActivateModel(Network *net, Node *node)
  {
    if (node->isActive) return;
    node->isActive = 1;
    node->prevActiveModel = NULL;
    node->nextActiveModel = net->mpActiveModels;
    if (net->mpActiveModels != NULL) {
      net->mpActiveModels->prevActiveModel = node;
    }
    net->mpActiveModels = node;
  
    if (node->isActiveNode) {
      assert(node->mType & NT_Tee);
      return;
    }
  
    node->isActiveNode = 1;
    node->prevActiveNode = NULL;
    node->nextActiveNode = net->mpActiveNodes;
    if (net->mpActiveNodes != NULL) {
      net->mpActiveNodes->prevActiveNode = node;
    }
    net->mpActiveNodes = node;
  
    assert(!HasCycle(net));
  
    MarkWordNodesLeadingFrom(net, node);
    ActivateWordNodesLeadingFrom(net, node);
  
    assert(AllWordSuccessorsAreActive(net));
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  DeactivateWordNodesLeadingFrom(Network *net, Node *node)
  {
    int i, nlinks = FORWARD_PASS ? node->mNLinks : node->mNBackLinks;
    Link *links   = FORWARD_PASS ? node->mpLinks  : node->mpBackLinks;
  
    for (i = 0; i < nlinks; i++) {
      Node *lnode = links[i].mpNode;
      if (lnode->mType & NT_Model && !(lnode->mType & NT_Tee)) continue;
      assert(!(lnode->mType & NT_Tee) || lnode->isActiveNode);
      if (--lnode->isActiveNode) continue;
  
      if (lnode->mType & NT_Tee && lnode->isActive) return;
  
      DeactivateWordNodesLeadingFrom(net, lnode);
  
      assert(lnode->prevActiveNode);
      lnode->prevActiveNode->nextActiveNode = lnode->nextActiveNode;
      if (lnode->nextActiveNode) {
        lnode->nextActiveNode->prevActiveNode = lnode->prevActiveNode;
      }
    }
    assert(!HasCycle(net));
    assert(AllWordSuccessorsAreActive(net));
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  DeactivateModel(Network *net, Node *node)
  {
    if (!node->isActive) return;
    node->isActive = 0;
  
    if (node->nextActiveModel != NULL) {
      node->nextActiveModel->prevActiveModel = node->prevActiveModel;
    }
  
    if (node->prevActiveModel != NULL) {
      node->prevActiveModel->nextActiveModel = node->nextActiveModel;
    } else {
      assert(net->mpActiveModels == node);
      net->mpActiveModels = node->nextActiveModel;
    }
  
    assert(!HasCycle(net));
    if (node->mType & NT_Tee && node->isActiveNode) return;
  
    node->isActiveNode = 0;
    DeactivateWordNodesLeadingFrom(net, node);
  
    if (node->nextActiveNode != NULL) {
      node->nextActiveNode->prevActiveNode = node->prevActiveNode;
    }
  
    if (node->prevActiveNode != NULL) {
      node->prevActiveNode->nextActiveNode = node->nextActiveNode;
    } else {
      assert(net->mpActiveNodes == node);
      net->mpActiveNodes = node->nextActiveNode;
    }
  
    assert(!HasCycle(net));
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  TokenPropagationInit(Network *net)
  {
    int i, j;
    Node *node;
  
    InitLogMath();
  
  //  for (i=0; i < net->nnodes; i++) {
  //  node = &net->mpNodes[i];
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      int numOfTokens = (node->mType & NT_Model) ? node->mpHmm->mNStates : 1;
  
      for (j=0; j < numOfTokens; j++) {
        node->tokens[j].mLike = LOG_0;
        node->tokens[j].wlr = NULL;
  #ifdef bordel_staff
        node->tokens[j].twlr = NULL;
  #endif // bordel_staff
      }
      node->isActive = 0;
      node->isActiveNode = 0;
    }
  
    if (net->mpOutPCache != NULL) {
      for (i = 0; i < net->mpModelSet->mNStates; i++) {
        net->mpOutPCache[i].mTime = UNDEF_TIME;
      }
    }
  
    if (net->mpMixPCache != NULL) {
      for (i = 0; i < net->mpModelSet->mNMixtures; i++) {
        net->mpMixPCache[i].mTime = UNDEF_TIME;
      }
    }
  
    net->mpBestToken  = NULL;
    net->mWordThresh = LOG_MIN;
    net->mpActiveModels = NULL;
    net->mpActiveNodes = NULL;
    net->mActiveTokens = 0;
  //  net->mTime = 0;
    node = FORWARD_PASS ? net->mpFirst : net->mpLast;
    node->tokens[0].mLike = 0;
    node->tokens[0].mAccuracy.logvalue = LOG_0;
    node->tokens[0].mAccuracy.negative = 0;
    node->tokens[0].wlr = NULL;
  #ifdef bordel_staff
    node->tokens[0].twlr = NULL;
    node->tokens[0].bestlike = 0;
  #endif
  
    if (net->mCollectAlphaBeta && FORWARD_PASS) {
      for (node = net->mpFirst; node != NULL; node = node->mpNext) {
        if (!(node->mType & NT_Model)) continue;
        node->alphaBetaList = node->alphaBetaListReverse = NULL;
      }
    }
    // Needed to load last FWBWRs to alphaBetaList
    if (net->mCollectAlphaBeta && !FORWARD_PASS) {
      for (node = net->mpFirst; node != NULL; node = node->mpNext) {
        if (!(node->mType & NT_Model)) continue;
        BackwardPruning(net->mTime, node, node->mpHmm->mNStates-1);
      }
    }
  
  
    node = FORWARD_PASS ? net->mpFirst : net->mpLast;
    net->mpActiveNodes = node;
    node->prevActiveNode = node->nextActiveNode = NULL;
    node->isActiveNode = 1;
  
    MarkWordNodesLeadingFrom(net, node);
    ActivateWordNodesLeadingFrom(net, node);
  
    TokenPropagationInNetwork(net);
  
    DeactivateWordNodesLeadingFrom(net, node);
    node->isActiveNode = 0;
    if (node->prevActiveNode) {
      node->prevActiveNode->nextActiveNode = node->nextActiveNode;
    }
  
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  TokenPropagationInNetwork(Network *net)
  {
    Node *  node;
    Link *  links;
    int     i;
    int     nlinks;
  
    // Beam pruning is not active in backward pass. First, it is not necessary
    // since only token that fit into the forward pass beam are allowed (backward
    // pruning after forward pass). Second, it could be dangerous. For example,
    // when training from lattices where LM weights are pushed to the begining,
    // the lowest weight (the biggest penalty) will be on the wery first link.
    // If this weight was lower than minus pruning treshold, it would result in
    // always killing token in the first node, when pruning during the backward pass.
    //                                      |
    net->mBeamThresh = FORWARD_PASS && // <--'
                      net->mpBestToken &&
                      net->mpBestToken->mLike - net->mPruningThresh > LOG_MIN
                      ? net->mpBestToken->mLike - net->mPruningThresh : LOG_MIN;
  
    node = FORWARD_PASS ? net->mpLast : net->mpFirst;
    
    KillToken(node->exitToken);
  
  //  Node *Xnode = net->mpActiveNodes;
  /*  for (node = FORWARD_PASS ? net->mpFirst : net->mpLast;
        node != NULL;
        node = FORWARD_PASS ? node->mpNext : node->mpBackNext) { //*/
    for (node = net->mpActiveNodes; node != NULL; node = node->nextActiveNode) 
    {
      if ((node->mType & NT_Tee) && node->tokens[0].IsActive()) 
      {
  //        assert(node->isActiveNode || (node->mType & NT_Tee && node->isActive));
  //        for (Xnode = net->mpActiveNodes; Xnode && Xnode != node; Xnode = Xnode->nextActiveNode);
  //        assert(Xnode);
  
        if (!(net->mCollectAlphaBeta && !FORWARD_PASS && // backward pruning
            BackwardPruning(net->mTime, node, 0))) 
        {   // after forward pass
          Hmm *hmm = node->mpHmm;
          FLOAT transP = hmm->mpTransition->mpMatrixO[hmm->mNStates - 1];
    #ifdef TRACE_TOKENS
          printf("Tee model State 0 -> Exit State ");
    #endif
          net->PassTokenInModel(&node->tokens[0], node->exitToken,
                                          transP * net->mTranScale);
        }
      }
  
      if (node->exitToken->IsActive()) 
      {
  //      assert(node->isActiveNode || (node->mType & NT_Tee && node->isActive));
  //      for (Xnode = net->mpActiveNodes; Xnode && Xnode != node; Xnode = Xnode->nextActiveNode);
  //      assert(Xnode);
  
        if (node->mType & NT_Model) 
        {
          if (net->mCollectAlphaBeta) 
          {
            if (FORWARD_PASS) 
              WriteAlpha(net->mTime, node, node->mpHmm->mNStates-1, node->exitToken);
            else 
              WriteBeta(net->mTime, node, 0, node->exitToken);
          }
          node->exitToken->mLike += net->mMPenalty;
        } 
        else if (node->mType & NT && node->mpPronun != NULL) 
        {
          node->exitToken->mLike += net->mWPenalty +
                                  net->mPronScale * node->mpPronun->prob;
          /*if (node->exitToken->mLike < net->mWordThresh) {
            node->exitToken->mLike = LOG_0;
          }*/
        }
  
        if (node->exitToken->mLike > net->mBeamThresh) 
        {
          if (node->mType & NT && node->mpPronun != NULL &&
            net->mAlignment & WORD_ALIGNMENT) 
          {
            AddWordLinkRecord(node->exitToken, node, -1, net->mTime);
          } 
          else if (node->mType & NT_Model && net->mAlignment & MODEL_ALIGNMENT) 
          {
            AddWordLinkRecord(node->exitToken, node, -1, net->mTime);
          }
  
          nlinks = FORWARD_PASS ? node->mNLinks : node->mNBackLinks;
          links  = FORWARD_PASS ? node->mpLinks  : node->mpBackLinks;
  
          for (i = 0; i < nlinks; i++) 
          {
            FLOAT lmLike = links[i].mLike * net->mLmScale;
            
            if (node->exitToken->mLike + lmLike > net->mBeamThresh 
                && (/*links[i].mpNode->mStart == UNDEF_TIME ||*/
                      links[i].mpNode->mStart <= net->mTime) 
                      
                && (  links[i].mpNode->mStop  == UNDEF_TIME        ||
                      links[i].mpNode->mStop  >= net->mTime)
                      
                && (  net->mSearchPaths != SP_TRUE_ONLY || 
                      (node->mType & NT_True)                   || 
                      !(links[i].mpNode->mType & NT_Model))) 
            {
  #           ifdef TRACE_TOKENS
              printf("Node %d -> Node %d ", node->mAux, links[i].mpNode->mAux);
  #           endif
              net->PassTokenInNetwork(node->exitToken,
                                      &links[i].mpNode->tokens[0], lmLike);
              if (links[i].mpNode->mType & NT_Model) 
              {
                ActivateModel(net, links[i].mpNode);
              } 
              else 
              {
                assert(links[i].mpNode->isActiveNode ||
                      links[i].mpNode == (FORWARD_PASS ? net->mpLast : net->mpFirst));
              }
            }
          }
        }
        
        if (!(node->mType & NT_Sticky)) 
          KillToken(node->exitToken);
      }
  
    }
  
  
    if (net->mCollectAlphaBeta) {
      if (FORWARD_PASS) {
        for (node = net->mpActiveModels; node != NULL; node = node->nextActiveModel) {
          if (/*!(node->mType & NT_Model) ||*/ node->tokens[0].mLike < LOG_MIN) continue;
          WriteAlpha(net->mTime, node, 0, &node->tokens[0]);
        }
      } else {
        for (node = net->mpActiveModels; node != NULL; node = node->nextActiveModel) {
          if (/*!(node->mType & NT_Model) ||*/ node->tokens[0].mLike < LOG_MIN
            || BackwardPruning(net->mTime, node, node->mpHmm->mNStates-1)) continue;
          WriteBeta(net->mTime, node, node->mpHmm->mNStates-1, &node->tokens[0]);
        }
      }
    }
  
  //  Go through newly activeted models and costruct list of active nodes
  //  for (node=net->mpActiveModels; node&&!node->isActive; node=node->nextActiveModel) {
  //    ActivateNode(net, node);
  //  }
    assert(!HasCycle(net));
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  TokenPropagationInModels(Network *net, FLOAT *observation)
  {
    Node *  node;
    Hmm *   hmm;
    size_t  winingToken = 0;
    size_t  i;
    size_t  j;
    int     from;
    int     to;
  //  int estate_id;
    int     state_idx;
  //  FLOAT threshOutProb = LOG_0;
  
    net->mpBestToken = NULL;
  /*  if (net->mpThreshState) {
      threshOutProb = net->OutputProbability(net->mpThreshState, observation, net);
    } */
  
    for (node = net->mpActiveModels; node != NULL; node = node->nextActiveModel) {
  //  for (node = net->mpFirst; node != NULL; node = node->mpNext) {
  //    if (!(node->mType & NT_Model)) continue;
  
      hmm = node->mpHmm;
  
      if (    (/*node->mStart != UNDEF_TIME &&*/node->mStart >= net->mTime)
          ||  (  node->mStop  != UNDEF_TIME &&  node->mStop  <  net->mTime)
          ||  net->mSearchPaths == SP_TRUE_ONLY && !(node->mType & NT_True)) 
      {
        for (i = 0; i < hmm->mNStates-1; i++) 
        {
          KillToken(&node->tokens[i]);
        }
        
        DeactivateModel(net, node);
        continue;
      }
  
      if (net->mAccumType == AT_MPE && FORWARD_PASS && node->tokens[0].IsActive()) 
      {
        FloatInLog fil_lmpa =
          {node->tokens[0].mLike + log(fabs(node->phoneAccuracy)),
          node->phoneAccuracy < 0};
  
        node->tokens[0].mAccuracy = FIL_Add(node->tokens[0].mAccuracy, fil_lmpa);
      }
  
      for (i = 0; i < hmm->mNStates-1; i++) 
      {
        assert(!net->mpAuxTokens[i].IsActive());
        if (node->tokens[i].IsActive()) 
        {
          net->mpAuxTokens[i] = node->tokens[i];
          node->tokens[i].mLike = LOG_0;
          node->tokens[i].wlr  = NULL;
        }
      }
  
      int keepModelActive = FALSE;
  
      assert(!node->tokens[hmm->mNStates-1].IsActive());
  
      for (j = 1; j < hmm->mNStates-1; j++) 
      {
        state_idx = (FORWARD_PASS ? j : hmm->mNStates-1 - j);
  
        if (net->mCollectAlphaBeta &&
          !FORWARD_PASS &&
          BackwardPruning(net->mTime, node, state_idx)) 
        {
          continue; // backward pruning after forward pass
        }
  
        for (i = 0; i < hmm->mNStates-1; i++) 
        {
          from = FORWARD_PASS ? i : hmm->mNStates-1 - j;
          to   = FORWARD_PASS ? j : hmm->mNStates-1 - i;
  
          assert(!net->mpAuxTokens[i].IsActive() || node->isActive);
  
          if (hmm->mpTransition->mpMatrixO[from * hmm->mNStates + to] > LOG_MIN &&
              net->mpAuxTokens[i].IsActive()) 
          {
            FLOAT transP = hmm->mpTransition->mpMatrixO[from * hmm->mNStates + to];
  
  #ifdef TRACE_TOKENS
            printf("Model %d State %d -> State %d ",  node->mAux, i, j);
  #endif
            if (net->PassTokenInModel(&net->mpAuxTokens[i], &node->tokens[j],
                                    transP * net->mTranScale)) {
              winingToken = i;
            }
          }
        }
  
        // if (IS_ACTIVE(node->tokens[j])) {
        if (node->tokens[j].mLike > net->mBeamThresh) 
        {
          FLOAT outProb = net->OutputProbability(hmm->mpState[state_idx-1],
                                                    observation, net);
          outProb *= net->mOutpScale;
  
          /*if (outProb < threshOutProb) {
            outProb = threshOutProb;
          }*/
  
          if (net->mCollectAlphaBeta && !FORWARD_PASS)
            WriteBeta(net->mTime, node, state_idx, &node->tokens[j]);
          
          if (net->mAccumType == AT_MFE && node->mType & NT_True) 
          {
            FloatInLog fil_like = {node->tokens[j].mLike, 0};
            node->tokens[j].mAccuracy = FIL_Add(node->tokens[j].mAccuracy, fil_like);
          }
  
          node->tokens[j].mAccuracy.logvalue += outProb;
          node->tokens[j].mLike              += outProb;
  
          if (net->mCollectAlphaBeta && FORWARD_PASS) 
            WriteAlpha(net->mTime, node, state_idx, &node->tokens[j]);
  
          if (net->mAlignment & STATE_ALIGNMENT && winingToken > 0 &&
            (winingToken != j || net->mAlignment & FRAME_ALIGNMENT)) 
          {
            AddWordLinkRecord(&node->tokens[j], node,
                              (FORWARD_PASS ? winingToken
                                            : hmm->mNStates-1 - winingToken)-1,
                              net->mTime-1);
          }
  
          net->mActiveTokens++;
          keepModelActive = TRUE;
          assert(node->isActive);
        } 
        else 
        {
          assert(node->isActive || !node->tokens[j].IsActive());
          KillToken(&node->tokens[j]);
        }
      }
  
      for (i = 0; i < hmm->mNStates-1; i++) 
      {
        KillToken(&net->mpAuxTokens[i]);
      }
  
      if (!keepModelActive) 
        DeactivateModel(net, node);
  
      state_idx = (FORWARD_PASS ? hmm->mNStates - 1 : 0);
      assert(!node->tokens[hmm->mNStates - 1].IsActive());
  
      if (!keepModelActive ||
        (net->mCollectAlphaBeta && !FORWARD_PASS &&
          BackwardPruning(net->mTime-1, node, state_idx))) 
      {
        // backward pruning after forward pass
        continue;
        //KillToken(&node->tokens[hmm->mNStates - 1]);
      }
  
      for (i = 1; i < hmm->mNStates-1; i++) 
      {
        from = FORWARD_PASS ? i : 0;
        to   = FORWARD_PASS ? hmm->mNStates-1 : hmm->mNStates-1 - i;
  
        if (node->tokens[i].IsActive()) 
        {
          if (!net->mpBestToken || net->mpBestToken->mLike < node->tokens[i].mLike) 
          {
            net->mpBestToken = &node->tokens[i];
            net->mpBestNode  = node;
          }
  
          if (hmm->mpTransition->mpMatrixO[from * hmm->mNStates + to] > LOG_MIN) 
          {
            FLOAT transP = hmm->mpTransition->mpMatrixO[from * hmm->mNStates + to];
  #ifdef TRACE_TOKENS
            printf("Model %d State %d -> Exit State ",  node->mAux, i);
  #endif
            if (net->PassTokenInModel(&node->tokens[i],
                                    &node->tokens[hmm->mNStates - 1],
                                    transP * net->mTranScale)) 
            {
              winingToken = i;
            }
          }
        }
      }
  
      if (node->tokens[hmm->mNStates - 1].IsActive()) 
      {
        if (net->mAccumType == AT_MPE && !FORWARD_PASS) 
        {
          FloatInLog fil_lmpa =
            {node->tokens[hmm->mNStates - 1].mLike + log(fabs(node->phoneAccuracy)),
            node->phoneAccuracy < 0};
  
          node->tokens[hmm->mNStates - 1].mAccuracy =
            FIL_Add(node->tokens[hmm->mNStates - 1].mAccuracy, fil_lmpa);
        }
  
    //    ActivateNode(net, node);
        if (net->mAlignment & STATE_ALIGNMENT) 
        {
          AddWordLinkRecord(&node->tokens[hmm->mNStates - 1], node,
                            (FORWARD_PASS ? winingToken :
                                            hmm->mNStates-1 - winingToken-1)-1,
                            net->mTime);
        }
      }
    }
    assert(!HasCycle(net));
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  TokenPropagationDone(Network *net)
  {
    int j;
    Node *node;
  
    KillToken(FORWARD_PASS ? net->mpLast->exitToken
                          : net->mpFirst->exitToken);
  
  //  for (i=0; i < net->nnodes; i++) {
  //    node = &net->mpNodes[i];
      for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      int numOfTokens;
  
      if (!(node->mType & NT_Model)) {
        assert(!node->exitToken->IsActive());
        continue;
      }
  
      numOfTokens = node->mpHmm->mNStates;
  
      for (j=0; j < numOfTokens; j++) {
        KillToken(&node->tokens[j]);
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  PassTokenMax(Token *from, Token *to, FLOAT mLike)
  {
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
      if (to->wlr) 
      {
        to->wlr->refs++;
      }
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
  
    if (!IS_ACTIVE(*to) || from->mLike + mLike > to->bestlike) {
      KillToken(to);
  
      ret = 1;
      *to = *from;
      to->bestlike = from->mLike + mLike;
      if (to->wlr) {
        to->wlr->refs++;
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
  PassTokenSum(Token *from, Token *to, FLOAT mLike)
  {
    double tl;
    FloatInLog fe, fil_from_like = {mLike, 0};
    int ret = 0;
  #ifdef TRACE_TOKENS
    printf("(%.2f + %.2f -> %.2f = ", from->mLike, mLike, to->mLike);
  #endif
    if (to->IsActive()) 
    {
      tl = LogAdd(to->mLike, from->mLike + mLike);
      fe = FIL_Add(to->mAccuracy, FIL_Mul(from->mAccuracy, fil_from_like));
    } 
    else 
    {
      tl = from->mLike + mLike;
      fe = FIL_Mul(from->mAccuracy, fil_from_like);
    }
  
    if (!to->IsActive() || from->mLike + mLike > to->bestlike) 
    {
      KillToken(to);
  
      ret = 1;
      *to = *from;
      to->bestlike = from->mLike + mLike;
      
      if (to->wlr) 
      {
        to->wlr->refs++;
      }
    }
  
    to->mLike   = tl;
    to->mAccuracy = fe;
  #ifdef TRACE_TOKENS
    printf("%.2f)\n", to->mLike);
  #endif
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  FreeWordLinkRecords(WLR *wlr)
  {
    if (wlr != NULL) {
      --wlr->refs;
      assert(wlr->refs >= 0);
  
      if (wlr->refs == 0) {
        FreeWordLinkRecords(wlr->mpNext);
  #ifdef DEBUG_MSGS
        assert(wlr->freed == 0);
        wlr->freed = 1;
  #else
        free(wlr);
  #endif
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  KillToken(Token *token)
  {
    token->mLike = LOG_0;
    FreeWordLinkRecords(token->wlr);
    token->wlr  = NULL;
    token->twlr = NULL;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  AddWordLinkRecord(Token *token, Node *node, int state_idx, int time)
  {
    WLR *wlr;
  
    if ((wlr = (WLR *) malloc(sizeof(WLR))) == NULL) {
      Error("Insufficient memory");
    }
  
    wlr->mStateIdx = state_idx;
    wlr->mLike  = token->mLike;
    wlr->mpNode  = node;
    wlr->mTime  = time;
    wlr->mpNext  = token->wlr;
    wlr->refs  = 1;
    token->wlr  = wlr;
  #ifdef bordel_staff
    if (!token->twlr) token->twlr = wlr;
  #endif // bordel_staff
  #ifdef DEBUG_MSGS
    wlr->mpTmpNext = firstWLR;
    firstWLR = wlr;
    wlr->freed = 0;
  #endif
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  DiagCGaussianDensity(Mixture *mix, FLOAT  *obs, Network *net)
  {
    FLOAT mLike = 0.0;
    int j;
  
    if (net && net->mpMixPCache[mix->mID].mTime == net->mTime) {
      return net->mpMixPCache[mix->mID].mValue;
    }
  
    for (j = 0; j < mix->mpMean->mVectorSize; j++) {
      mLike += SQR(obs[j] - mix->mpMean->mpVectorO[j]) * mix->mpVariance->mpVectorO[j];
    }
  
    mLike = -0.5 * (mix->mGConst + mLike);
  
    if (net) {
      net->mpMixPCache[mix->mID].mTime  = net->mTime;
      net->mpMixPCache[mix->mID].mValue = mLike;
    }
  
    return mLike;
  }
  
  #ifdef DEBUG_MSGS
  int gaus_computaions = 0;
  #endif
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  DiagCGaussianMixtureDensity(State *state, FLOAT *obs, Network *net)
  {
    size_t  i;
    FLOAT   mLike = LOG_0;
  
    assert(state->mOutPdfKind == KID_DiagC);
  
    if (net && net->mpOutPCache[state->mID].mTime == net->mTime) {
      return net->mpOutPCache[state->mID].mValue;
    }
  
    for (i = 0; i < state->mNumberOfMixtures; i++) {
      FLOAT glike;
      Mixture *mix = state->mpMixture[i].mpEstimates;
  
      obs = XFormPass(mix->mpInputXForm, obs,
                      net ? net->mTime : UNDEF_TIME,
                      net ? net->mPropagDir : FORWARD);
  
      assert(obs != NULL);
      glike = DiagCGaussianDensity(mix, obs, net);
      mLike  = LogAdd(mLike, glike + state->mpMixture[i].mWeight);
    }
  
  #ifdef DEBUG_MSGS
    gaus_computaions++;
  #endif
  
    if (net) {
      net->mpOutPCache[state->mID].mTime = net->mTime;
      net->mpOutPCache[state->mID].mValue = mLike;
    }
    return mLike;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  FromObservationAtStateId(State *state, FLOAT *obs, Network *net)
  {
    obs = XFormPass(net->mpModelSet->mpInputXForm, obs,
                    net ? net->mTime : UNDEF_TIME,
                    net ? net->mPropagDir : FORWARD);
    assert(obs != NULL);
    return obs[state->PDF_obs_coef];
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Label *
  GetLabels(Token *token)
  {
    WLR *wlr;
    Label *tmp, *level[3] = {NULL, NULL, NULL};
    int li = 0;
  
    if (!token || !token->IsActive())
      return NULL;
    
    for (wlr = token->wlr; wlr != NULL; wlr = wlr->mpNext) 
    {
      if ((tmp = (Label *) malloc(sizeof(Label))) == NULL)
        Error("Insufficient memory");
  
      tmp->mScore = wlr->mLike;
      tmp->mStop  = wlr->mTime;
      tmp->mId    = wlr->mStateIdx;
  
      if (wlr->mpNode->mType & NT_Model) 
      {
        li = wlr->mStateIdx >= 0 ? 0 : 1;
        tmp->mpData = wlr->mpNode->mpHmm;
        tmp->mpName = wlr->mpNode->mpHmm->mpMacro->mpName;
        tmp->mpNextLevel = level[li+1] ? level[li+1] : level[2];
      } 
      else //if (wlr->mpNode->mpPronun->outSymbol) 
      {
        li = 2;
        tmp->mpData = wlr->mpNode->mpPronun->word;
        tmp->mpName = wlr->mpNode->mpPronun->outSymbol;
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
            level[li-1]->mStart  = tmp2->mStop;
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
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ViterbiInit(Network *net)
  {
    net->PassTokenInModel   = &PassTokenMax;
    net->PassTokenInNetwork = &PassTokenMax;
    net->mPropagDir = FORWARD;
    net->mpModelSet->ResetXFormInstances();
  
    net->mTime = 0; // Must not be set to -net->mpModelSet->mTotalDelay yet
                  // otherwise token cannot enter first model node
                  // with start set to 0
    TokenPropagationInit(net);
    net->mTime = -net->mpModelSet->mTotalDelay;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  #ifdef DEBUG_MSGS
  void 
  PrintNumOfWLR() 
  {
    WLR *wlr = firstWLR;
    int nwlrs = 0, nwlrsnf = 0;
    while (wlr) {
      nwlrs++;
      if (!wlr->freed) {
        nwlrsnf++;
      }
      wlr = wlr->mpTmpNext;
    }
    printf("%d Released: %d\n", nwlrsnf, nwlrs - nwlrsnf);
  }
  #endif
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ViterbiStep(Network *net, FLOAT *observation)
  {
    net->mTime++;
    net->mpModelSet->UpdateStacks(observation, net->mTime, net->mPropagDir);
  
    if (net->mTime <= 0) {
      return;
    }
  
    TokenPropagationInModels(net, observation);
    TokenPropagationInNetwork(net);
  #ifdef DEBUG_MSGS
    printf("Frame: %ld Nuberm of WLR Active: ", net->mTime); PrintNumOfWLR();
  #endif
  
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  ViterbiDone(Network *net, Label **labels)
  {
    FLOAT totLike = LOG_0;
    if (labels) 
    {
      if (net->mpLast->exitToken && 
          net->mpLast->exitToken->IsActive()) 
      {
        totLike = net->mpLast->exitToken->mLike;
        *labels = net->mpLast->exitToken->pGetLabels();
      } 
      else 
      {
        *labels = NULL;
      }
    }
    
    TokenPropagationDone(net);
  
  #ifdef DEBUG_MSGS
    printf("Number of WLR Unreleased: "); PrintNumOfWLR();
    printf("Number of output prob. computations: %d\n", gaus_computaions);
  #endif
    return totLike;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  /*void BaumWelchInit(ModelSet *hmms) 
  {
    //ResetAccumsForHMMSet(hmms);
    hmms.ResetAccums();
  }*/
  
  struct FWBWRet {
    double totLike;
    FloatInLog avgAccuracy;
  };
  
  
  //***************************************************************************
  //***************************************************************************
  struct FWBWRet 
  ForwardBackward(Network *net, FLOAT *obsMx, int nFrames)
  {
    int         i;
    Cache *     outPCache;
    struct      FWBWRet ret;
    ModelSet *  hmms = net->mpModelSet;
  
    outPCache = (Cache *) malloc(nFrames * hmms->mNStates * sizeof(Cache));
  
    if (outPCache == NULL) {
      Error("Insufficient memory");
    }
  
    for (i = 0; i < nFrames * hmms->mNStates; i++) {
      outPCache[i].mTime  = UNDEF_TIME;
      outPCache[i].mValue = LOG_0;
    }
  
    free(net->mpOutPCache);
    net->mpOutPCache = NULL;
  
    net->PassTokenInModel   = &PassTokenSum;
    net->PassTokenInNetwork = &PassTokenSum;
    net->mAlignment          = NO_ALIGNMENT;
  
  
    //Forward Pass
    net->mPropagDir = FORWARD;
    net->mCollectAlphaBeta = 1;
    net->mTime = -hmms->mTotalDelay;
  
    net->mpModelSet->ResetXFormInstances();
    for (i = 0; i < hmms->mTotalDelay; i++) {
      net->mTime++;
      hmms->UpdateStacks(obsMx + hmms->mInputVectorSize * i, net->mTime, FORWARD);
    }
  
    
    // tady nekde je chyba
    // !!!!
    TokenPropagationInit(net);
  
    for (i = hmms->mTotalDelay; i < nFrames+hmms->mTotalDelay; i++) {
      net->mTime++;
      net->mpOutPCache = outPCache + hmms->mNStates * (net->mTime-1);
  
      net->mpModelSet->UpdateStacks(obsMx + hmms->mInputVectorSize * i,
                                net->mTime, net->mPropagDir);
  
      TokenPropagationInModels(net,  obsMx + hmms->mInputVectorSize * i);
      TokenPropagationInNetwork(net);
    }
  
    if (!net->mpLast->exitToken->IsActive()) { // No token survivered
      TokenPropagationDone(net);
      FreeFWBWRecords(net);
      net->mpOutPCache = outPCache;
      ret.totLike = LOG_0;
      return ret;
    }
  
    ret.totLike = net->mpLast->exitToken->mLike; //  totalLikelihood;
    TokenPropagationDone(net);
  
    //Backward Pass
    net->mPropagDir = BACKWARD;
    net->mTime = nFrames+hmms->mTotalDelay;
  
    for (i = nFrames + hmms->mTotalDelay - 1; i >= nFrames; i--) {
  //  We do not need any features, in backward prop. All output probab. are cached.
  //  UpdateStacks(hmms, obsMx + hmms->mInputVectorSize * i,
  //                 net->mTime, BACKWARD);
      net->mTime--;
    }
  
    net->mpOutPCache = NULL;  // TokenPropagationInit would reset valid 1st frm cache
    TokenPropagationInit(net);
  
    for (i = nFrames-1; i >= 0; i--) {
      net->mpOutPCache = outPCache + hmms->mNStates * (net->mTime-1);
  
  //  We do not need any features, in backward prop. All output probab. are cached.
  //  UpdateStacks(net->mpModelSet, obsMx + hmms->mInputVectorSize * i,     |
  //                 net->mTime, net->mPropagDir);                   |
  //                                                               V
      TokenPropagationInModels(net,  NULL); //obsMx + hmms->mInputVectorSize * i);
      net->mTime--;
      TokenPropagationInNetwork(net);
    }
  
    net->mpOutPCache = outPCache;
  
    if (!net->mpFirst->exitToken->IsActive()) { // No token survivered
      TokenPropagationDone(net);
      FreeFWBWRecords(net);
      ret.totLike = LOG_0;
      return ret;
    }
  
    ret.totLike = HIGHER_OF(ret.totLike, net->mpFirst->exitToken->mLike); //  totalLikelihood;
    // Backward pass P can differ from forward pass P because of the precision
    // problems. Take the higher one to decrease the possibility of getting
    // an occupation probability (when normalizing by P) higher that one.
  
    FloatInLog fil_ret_totLike = {ret.totLike, 0};
    ret.avgAccuracy  = FIL_Div(net->mpFirst->exitToken->mAccuracy, fil_ret_totLike);
    TokenPropagationDone(net);
  
    // There may be remaining records in alphaBetaListReverse unused in
    // backward pass. Free them and set alphaBetaListReverse's to NULLs;
  
    Node *node;
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (!(node->mType & NT_Model)) continue;
  
      while (node->alphaBetaListReverse) {
        FWBWR *fwbwr = node->alphaBetaListReverse;
        node->alphaBetaListReverse = fwbwr->mpNext;
        free(fwbwr);
      }
    }
  
    return ret;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT 
  MCEReest(Network *net, FLOAT *obsMx, FLOAT *obsMx2, int nFrames, FLOAT weight, FLOAT sigSlope)
  {
    struct FWBWRet fwbw;
    FLOAT TP, P, F;
  
  //  FLOAT updateDir;
    int i, j, k;
    int t;
    ModelSet *hmmsAlig = net->mpModelSet;
    ModelSet *hmmsUpdt = net->mpModelSetToUpdate;
    Node *node;
  
    net->mAccumType = AT_ML;
    net->PassTokenInModel   = &PassTokenSum;
    net->PassTokenInNetwork = &PassTokenSum;
    net->mPropagDir          = FORWARD;
    net->mAlignment          = NO_ALIGNMENT;
  
    net->mSearchPaths        = SP_TRUE_ONLY;
    net->mpModelSet->ResetXFormInstances();
  
    net->mTime = 0; // Must not be set to -net->mpModelSet->totalDelay yet
                  // otherwise token cannot enter first model node
                  // with start set to 0
    TokenPropagationInit(net);
    net->mTime = -net->mpModelSet->mTotalDelay;
  
    for (t = 0; t < nFrames+hmmsAlig->mTotalDelay; t++) {
      ViterbiStep(net, obsMx + hmmsAlig->mInputVectorSize * t);
    }
  
    TP = net->mpLast->exitToken->mLike;
    ViterbiDone(net, NULL);
  
    if (TP <= LOG_MIN) return LOG_0;
  
  
    ////////////////// Denominator accumulation //////////////////
    net->mSearchPaths = SP_ALL;
  
    fwbw = ForwardBackward(net, obsMx, nFrames);
    P = fwbw.totLike;
  
    assert(P >= TP);
  
    F = TP - LogSub(P, TP);
  printf("MCE distance: %g; ", F);
    F = exp(-sigSlope * F);
    F = (sigSlope*F) / SQR(1+F);
  printf("weight: %g\n", F);
    weight *= F;
  
    if (P < LOG_MIN) return LOG_0;
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (node->mType & NT_Model && node->alphaBetaList != NULL &&
        node->alphaBetaList->mTime == 0) {
        node->alphaBetaListReverse = node->alphaBetaList;
        node->alphaBetaList = node->alphaBetaList->mpNext;
      }
    }
  
    hmmsAlig->ResetXFormInstances();
  
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->ResetXFormInstances();
    }
  
    for (i = 0; i < hmmsAlig->mTotalDelay; i++) {
      hmmsAlig->UpdateStacks(obsMx+hmmsAlig->mInputVectorSize*i,
                  i-hmmsAlig->mTotalDelay, FORWARD);
    }
  
    if (hmmsAlig != hmmsUpdt) {
      for (i = 0; i < hmmsUpdt->mTotalDelay; i++) {
        hmmsUpdt->UpdateStacks(obsMx2+hmmsUpdt->mInputVectorSize*i,
                    i-hmmsUpdt->mTotalDelay, FORWARD);
      }
    }
  
  
    // net->mpMixPCache might be used to cache likelihoods of mixtures of target
    // models. Reallocate the cache to fit mixtures of both models and reset it.
    k = HIGHER_OF(hmmsUpdt->mNMixtures, hmmsAlig->mNMixtures);
    net->mpMixPCache = (Cache *) realloc(net->mpMixPCache, k * sizeof(Cache));
    if (net->mpMixPCache == NULL) Error("Insufficient memory");
  
    for (i = 0; i < k; i++) net->mpMixPCache[i].mTime = UNDEF_TIME;
  
  // Update accumulators
    for (net->mTime = 0; net->mTime < nFrames; net->mTime++) {//for every frame
      FLOAT *obs  =obsMx +hmmsAlig->mInputVectorSize*(net->mTime+hmmsAlig->mTotalDelay);
      FLOAT *obs2 =obsMx2+hmmsUpdt->mInputVectorSize*(net->mTime+hmmsUpdt->mTotalDelay);
      hmmsAlig->UpdateStacks(obs, net->mTime, FORWARD);
      if (hmmsAlig != hmmsUpdt) {
        hmmsUpdt->UpdateStacks(obs2, net->mTime, FORWARD);
      }
  
      for (node = net->mpFirst; node != NULL; node = node->mpNext) { //for every model
  //    for (k=0; k < net->nnodes; k++) {
  //      Node *node = &net->mpNodes[k];
        if (node->mType & NT_Model &&
          node->alphaBetaList != NULL &&
          node->alphaBetaList->mTime == net->mTime+1) {
  
          struct AlphaBeta *st;
          int Nq       = node->mpHmm->mNStates;
          st = node->alphaBetaList->mpState;
  
          for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
            if (st[j].alpha + st[j].beta - P > MIN_LOG_WEGIHT) {
              assert(node->alphaBetaListReverse->mTime == net->mTime);
  
              ReestState(net, node, j-1,
                        (st[j].alpha + st[j].beta - P)  * net->mOcpScale,
                          -weight, obs, obs2);
            }
          }
          if (node->alphaBetaListReverse) free(node->alphaBetaListReverse);
          node->alphaBetaListReverse = node->alphaBetaList;
          node->alphaBetaList = node->alphaBetaList->mpNext;
        }
      }
    }
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (node->alphaBetaListReverse != NULL)
        free(node->alphaBetaListReverse);
    }
  
  
    ////////////////// Numerator accumulation //////////////////
    net->mSearchPaths = SP_TRUE_ONLY;
  
  
    ForwardBackward(net, obsMx, nFrames);
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (node->mType & NT_Model && node->alphaBetaList != NULL &&
        node->alphaBetaList->mTime == 0) {
        node->alphaBetaListReverse = node->alphaBetaList;
        node->alphaBetaList = node->alphaBetaList->mpNext;
      }
    }
  
    hmmsAlig->ResetXFormInstances();
  
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->ResetXFormInstances();
    }
  
    for (i = 0; i < hmmsAlig->mTotalDelay; i++) {
      hmmsAlig->UpdateStacks(obsMx+hmmsAlig->mInputVectorSize*i,
                  i-hmmsAlig->mTotalDelay, FORWARD);
    }
  
    if (hmmsAlig != hmmsUpdt) {
      for (i = 0; i < hmmsUpdt->mTotalDelay; i++) {
        hmmsUpdt->UpdateStacks(obsMx2+hmmsUpdt->mInputVectorSize*i,
                    i-hmmsUpdt->mTotalDelay, FORWARD);
      }
    }
  
  
    // net->mpMixPCache might be used to cache likelihoods of mixtures of target
    // models. Reallocate the cache to fit mixtures of both models and reset it.
    k = HIGHER_OF(hmmsUpdt->mNMixtures, hmmsAlig->mNMixtures);
    net->mpMixPCache = (Cache *) realloc(net->mpMixPCache, k * sizeof(Cache));
    if (net->mpMixPCache == NULL) Error("Insufficient memory");
  
    for (i = 0; i < k; i++) net->mpMixPCache[i].mTime = UNDEF_TIME;
  
  // Update accumulators
    for (net->mTime = 0; net->mTime < nFrames; net->mTime++) {//for every frame
      FLOAT *obs  =obsMx +hmmsAlig->mInputVectorSize*(net->mTime+hmmsAlig->mTotalDelay);
      FLOAT *obs2 =obsMx2+hmmsUpdt->mInputVectorSize*(net->mTime+hmmsUpdt->mTotalDelay);
      hmmsAlig->UpdateStacks(obs, net->mTime, FORWARD);
      if (hmmsAlig != hmmsUpdt) {
        hmmsUpdt->UpdateStacks(obs2, net->mTime, FORWARD);
      }
  
      for (node = net->mpFirst; node != NULL; node = node->mpNext) 
      { //for every model
        if (node->mType & NT_Model &&
          node->alphaBetaList != NULL &&
          node->alphaBetaList->mTime == net->mTime+1) 
        {
  
          struct AlphaBeta *st;
          int Nq       = node->mpHmm->mNStates;
          FLOAT *aq    = node->mpHmm->        mpTransition->mpMatrixO;
          FLOAT *aqacc = node->mpHmmToUpdate->mpTransition->mpMatrixO + SQR(Nq);
  //        int qt_1 = (net->mNumberOfNetStates * net->mTime) + node->estate_id;
  //        int qt = qt_1 + net->mNumberOfNetStates;
  
          st = node->alphaBetaList->mpState;
  
          if (//!net->mmi_den_pass &&
            st[Nq-1].alpha + st[Nq-1].beta - TP > MIN_LOG_WEGIHT) 
          {
            for (i = 0; i < Nq - 1; i++) 
            {
              LOG_INC(aqacc[i * Nq + Nq-1], aq[i * Nq + Nq-1]  * net->mTranScale +
                                          (st[i].alpha                         +
                                            st[Nq-1].beta - TP) * net->mOcpScale);
            }
          }
  
          for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
            if (st[j].alpha + st[j].beta - TP > MIN_LOG_WEGIHT) {
              FLOAT bjtO =net->mpOutPCache[hmmsAlig->mNStates * net->mTime +
                                        node->mpHmm->mpState[j-1]->mID].mValue;
              // ForwardBackward() set net->mpOutPCache to contain out prob. for all frames
  
              assert(node->alphaBetaListReverse->mTime == net->mTime);
  
  //            if (!net->mmi_den_pass) {
              for (i = 0; i < Nq - 1; i++) {
                LOG_INC(aqacc[i * Nq + j],
                        aq[i * Nq + j]    * net->mTranScale +
                        (node->alphaBetaListReverse->mpState[i].alpha +
                        bjtO              * net->mOutpScale +
                        st[j].beta - TP)   * net->mOcpScale);
              }
  //            }
  
              ReestState(net, node, j-1,
                        (st[j].alpha + st[j].beta - TP)  * net->mOcpScale,
                          weight, obs, obs2);
  
  // For True MCE
  //            ReestState(net, node, j-1,
  //                       (st[j].alpha + st[j].beta - TP + LogAdd(TP,P) - P)  * net->mOcpScale,
  //                        weight, obs, obs2);
  
            }
          }
          
          if (node->alphaBetaListReverse) 
            free(node->alphaBetaListReverse);
            
          node->alphaBetaListReverse = node->alphaBetaList;
          node->alphaBetaList = node->alphaBetaList->mpNext;
        }
      }
    }
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (node->alphaBetaListReverse != NULL)
        free(node->alphaBetaListReverse);
    }
  
    net->mAccumType = AT_MCE;
    return TP;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT BaumWelchReest(Network *net, FLOAT *obsMx, FLOAT *obsMx2, int nFrames, FLOAT weight)
  {
    struct FWBWRet fwbw;
    FLOAT P, updateDir;
    int i, j, k;
    ModelSet *hmmsAlig = net->mpModelSet;
    ModelSet *hmmsUpdt = net->mpModelSetToUpdate;
    Node *node;
  
    fwbw = ForwardBackward(net, obsMx, nFrames);
    P = fwbw.totLike;
    if (P < LOG_MIN) return LOG_0;
  
  #ifdef MOTIF
    FLOAT *ocprob = (FLOAT *) malloc(net->mNumberOfNetStates * (nFrames+1) * sizeof(FLOAT));
    for (i=0; i<net->mNumberOfNetStates * (nFrames+1); i++) ocprob[i] = 0;
  #endif
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (node->mType & NT_Model && node->alphaBetaList != NULL &&
        node->alphaBetaList->mTime == 0) {
        node->alphaBetaListReverse = node->alphaBetaList;
        node->alphaBetaList = node->alphaBetaList->mpNext;
      }
    }
  
    hmmsAlig->ResetXFormInstances();
  
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->ResetXFormInstances();
    }
  
    for (i = 0; i < hmmsAlig->mTotalDelay; i++) {
      hmmsAlig->UpdateStacks(obsMx+hmmsAlig->mInputVectorSize*i,
                            i-hmmsAlig->mTotalDelay, FORWARD);
    }
  
    if (hmmsAlig != hmmsUpdt) {
      for (i = 0; i < hmmsUpdt->mTotalDelay; i++) {
        hmmsUpdt->UpdateStacks(obsMx2+hmmsUpdt->mInputVectorSize*i,
                    i-hmmsUpdt->mTotalDelay, FORWARD);
      }
    }
  
  
    // net->mpMixPCache might be used to cache likelihoods of mixtures of target
    // models. Reallocate the cache to fit mixtures of both models and reset it.
    k = HIGHER_OF(hmmsUpdt->mNMixtures, hmmsAlig->mNMixtures);
    net->mpMixPCache = (Cache *) realloc(net->mpMixPCache, k * sizeof(Cache));
    
    if (net->mpMixPCache == NULL) 
      Error("Insufficient memory");
  
    for (i = 0; i < k; i++) net->mpMixPCache[i].mTime = UNDEF_TIME;
  
  // Update accumulators
    for (net->mTime = 0; net->mTime < nFrames; net->mTime++) {//for every frame
      FLOAT *obs  =obsMx +hmmsAlig->mInputVectorSize*(net->mTime+hmmsAlig->mTotalDelay);
      FLOAT *obs2 =obsMx2+hmmsUpdt->mInputVectorSize*(net->mTime+hmmsUpdt->mTotalDelay);
      hmmsAlig->UpdateStacks(obs, net->mTime, FORWARD);
      if (hmmsAlig != hmmsUpdt) {
        hmmsUpdt->UpdateStacks(obs2, net->mTime, FORWARD);
      }
  
      for (node = net->mpFirst; node != NULL; node = node->mpNext) { //for every model
  //    for (k=0; k < net->nnodes; k++) {
  //      Node *node = &net->mpNodes[k];
        if (node->mType & NT_Model &&
          node->alphaBetaList != NULL &&
          node->alphaBetaList->mTime == net->mTime+1) {
  
          struct AlphaBeta *st;
          int Nq       = node->mpHmm->mNStates;
          FLOAT *aq    = node->mpHmm->        mpTransition->mpMatrixO;
          FLOAT *aqacc = node->mpHmmToUpdate->mpTransition->mpMatrixO + SQR(Nq);
  //        int qt_1 = (net->mNumberOfNetStates * net->mTime) + node->estate_id;
  //        int qt = qt_1 + net->mNumberOfNetStates;
  
          st = node->alphaBetaList->mpState;
  
          if (//!net->mmi_den_pass &&
            st[Nq-1].alpha + st[Nq-1].beta - P > MIN_LOG_WEGIHT) {
            for (i = 0; i < Nq - 1; i++) {
              LOG_INC(aqacc[i * Nq + Nq-1], aq[i * Nq + Nq-1]  * net->mTranScale +
                                          (st[i].alpha                         +
                                            st[Nq-1].beta - P) * net->mOcpScale);
            }
          }
  
          for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
            if (st[j].alpha + st[j].beta - P > MIN_LOG_WEGIHT) {
  #ifdef MOTIF
            ocprob[net->mNumberOfNetStates * (net->mTime+1) + node->estate_id + j]
  //           = node->phoneAccuracy;
              = exp(st[j].alpha+st[j].beta-P) *
                ((1-2*st[j].alphaAccuracy.negative) * exp(st[j].alphaAccuracy.logvalue - st[j].alpha) +
                (1-2*st[j].betaAccuracy.negative)  * exp(st[j].betaAccuracy.logvalue  - st[j].beta)
                - (1-2*fwbw.avgAccuracy.negative) * exp(fwbw.avgAccuracy.logvalue)
                );
  
  #endif
  //            int qt_1   = qt - net->mNumberOfNetStates;
              FLOAT bjtO =net->mpOutPCache[hmmsAlig->mNStates * net->mTime +
                                        node->mpHmm->mpState[j-1]->mID].mValue;
              // ForwardBackward() set net->mpOutPCache to contain out prob. for all frames
  
              assert(node->alphaBetaListReverse->mTime == net->mTime);
  
  //            if (!net->mmi_den_pass) {
              for (i = 0; i < Nq - 1; i++) {
                LOG_INC(aqacc[i * Nq + j],
                        aq[i * Nq + j]    * net->mTranScale +
                        (node->alphaBetaListReverse->mpState[i].alpha +
                        bjtO              * net->mOutpScale +
                        st[j].beta - P)   * net->mOcpScale);
              }
  //            }
  
              if (net->mAccumType == AT_MFE || net->mAccumType == AT_MPE) {
                updateDir = (1-2*st[j].alphaAccuracy.negative) * exp(st[j].alphaAccuracy.logvalue - st[j].alpha) +
                            (1-2*st[j].betaAccuracy.negative)  * exp(st[j].betaAccuracy.logvalue  - st[j].beta)  -
                            (1-2*fwbw.avgAccuracy.negative)    * exp(fwbw.avgAccuracy.logvalue);
              } else {
                updateDir = 1.0;
              }
  
              ReestState(net, node, j-1,
                        (st[j].alpha + st[j].beta - P)  * net->mOcpScale,
                          updateDir*weight, obs, obs2);
            }
          }
          if (node->alphaBetaListReverse) free(node->alphaBetaListReverse);
          node->alphaBetaListReverse = node->alphaBetaList;
          node->alphaBetaList = node->alphaBetaList->mpNext;
        }
      }
    }
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) {
      if (node->alphaBetaListReverse != NULL)
        free(node->alphaBetaListReverse);
    }
  
  #ifdef MOTIF
    FLOAT max = LOG_0;
    printf("mTranScale: %f\noutpScale: %f\n",net->mTranScale, net->mOutpScale);
    for (i = 0; i < net->mNumberOfNetStates * (nFrames+1); i++) max = HIGHER_OF(max, ocprob[i]);
  
    imagesc(ocprob, net->mNumberOfNetStates, (nFrames+1),
            sizeof(FLOAT) == 4 ? "float" : "double",
            NULL,  cm_color, "OccProb");
  
    for (j=0; j<nFrames+1; j++) {
      for (i=0; i<net->mNumberOfNetStates; i++) {
        printf("%6.3g ", (ocprob[net->mNumberOfNetStates * j + i]));
      }
      printf("\n");
    }
  //  for (i = 0; i < net->mNumberOfNetStates * (nFrames+1); i++) {
  //    if (ocprob[i] < LOG_MIN) ocprob[i] = max;
  //  }
  //  imagesc(ocprob, net->mNumberOfNetStates, (nFrames+1), STR(FLOAT), NULL, cm_color, "Log OccProb");
    free(ocprob);
  #endif
  
  
    return P;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  FLOAT ViterbiReest(Network *net, FLOAT *obsMx, FLOAT *obsMx2, int nFrames, FLOAT weight)
  {
    int t;
    WLR *wlr;
    Node *prevnode = NULL;
    FLOAT P;
    Cache *outPCache;
    ModelSet *hmmsAlig = net->mpModelSet;
    ModelSet *hmmsUpdt = net->mpModelSetToUpdate;
  
    outPCache = (Cache *) malloc(nFrames * hmmsAlig->mNStates * sizeof(Cache));
    if (outPCache == NULL) Error("Insufficient memory");
  
    for (t = 0; t < nFrames * hmmsAlig->mNStates; t++) {
      outPCache[t].mTime  = UNDEF_TIME;
      outPCache[t].mValue = LOG_0;
    }
  
    free(net->mpOutPCache);
    net->mpOutPCache = NULL;
    net->mAlignment = STATE_ALIGNMENT;
    ViterbiInit(net);
    nFrames += hmmsAlig->mTotalDelay;
  
    for (t = 0; t < nFrames; t++) {
      if (t >= hmmsAlig->mTotalDelay) {
        net->mpOutPCache = outPCache + hmmsAlig->mNStates*(t-hmmsAlig->mTotalDelay);
      }
      ViterbiStep(net, obsMx + hmmsAlig->mInputVectorSize * t);
    }
  
    net->mpOutPCache = outPCache;
  
    if (!net->mpLast->exitToken->IsActive()) {
      ViterbiDone(net, NULL);
      return LOG_0;
    }
  
    hmmsAlig->ResetXFormInstances();
  
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->ResetXFormInstances();
    }
  
    // invert order of WRLs
    wlr = net->mpLast->exitToken->wlr;
    
    while (wlr->mpNext != NULL) 
    {
      WLR *tmp = wlr->mpNext->mpNext;
      wlr->mpNext->mpNext = net->mpLast->exitToken->wlr;
      net->mpLast->exitToken->wlr = wlr->mpNext;
      wlr->mpNext = tmp;
    }
  
    for (net->mTime = -hmmsAlig->mTotalDelay; net->mTime < 0; net->mTime++) 
    {
      FLOAT *obs = obsMx+hmmsAlig->mInputVectorSize*(net->mTime+hmmsAlig->mTotalDelay);
      hmmsAlig->UpdateStacks(obs, net->mTime, FORWARD);
    }
  
    if (hmmsAlig != hmmsUpdt) 
    {
      FLOAT *obs2= obsMx2+hmmsUpdt->mInputVectorSize*(net->mTime+hmmsUpdt->mTotalDelay);
      for (net->mTime = -hmmsUpdt->mTotalDelay; net->mTime < 0; net->mTime++) 
      {
        hmmsUpdt->UpdateStacks(obs2, net->mTime, FORWARD);
      }
    }
  
  // Update accumulators
    for (wlr = net->mpLast->exitToken->wlr; wlr != NULL; wlr = wlr->mpNext) 
    {
      Node *node   = wlr->mpNode;
      int Nq       = node->mpHmmToUpdate->mNStates;
      FLOAT *aqacc = node->mpHmmToUpdate->mpTransition->mpMatrixO + SQR(Nq);
      int currstate = wlr->mStateIdx+1;
      int nextstate = (wlr->mpNext && node == wlr->mpNext->mpNode)
                      ? wlr->mpNext->mStateIdx+1 : Nq-1;
      int duration  = wlr->mTime - net->mTime;
  
  
  
      if (prevnode != node) {
  //      if (!net->mmi_den_pass)
        LOG_INC(aqacc[currstate], 0 /*ln(1)*/);
        prevnode = node;
      }
  
  //    if (!net->mmi_den_pass) { // So far we dont do any MMI estimation of trasitions
      LOG_INC(aqacc[currstate * Nq + currstate], log(duration-1));
      LOG_INC(aqacc[currstate * Nq + nextstate], 0 /*ln(1)*/);
  //    }
  
      for (; net->mTime < wlr->mTime; net->mTime++) {
        //for every frame of segment
        FLOAT *obs =obsMx +hmmsAlig->mInputVectorSize*(net->mTime+hmmsAlig->mTotalDelay);
        FLOAT *obs2=obsMx2+hmmsUpdt->mInputVectorSize*(net->mTime+hmmsUpdt->mTotalDelay);
  
        hmmsAlig->UpdateStacks(obs, net->mTime, FORWARD);
        if (hmmsAlig != hmmsUpdt) {
          hmmsUpdt->UpdateStacks(obs2, net->mTime, FORWARD);
        }
        ReestState(net, node, currstate-1, 0.0, 1.0*weight, obs, obs2);
      }
    }
  
    P = net->mpLast->exitToken->wlr->mLike;
    ViterbiDone(net, NULL);
    return P;
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  UpdateXFormStatCache(XFormStatCache *xfsc,
                            XForm *topXForm,       //to locate positions in input vector
                            FLOAT *input) 
  {
    size_t    i;
    size_t    j;
    size_t    size = xfsc->mpXForm->mInSize;
  
    if (topXForm == NULL)
      return;
  
    if (topXForm == xfsc->mpXForm) 
    {
      if (xfsc->mNorm > 0) 
      {
        for (i=0; i < size; i++) 
        {
          xfsc->mpStats[i] += input[i];
          for (j=0; j <= i; j++) {
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
    
    else if (topXForm->mXFormType == XT_COMPOSITE) 
    {
      CompositeXForm *cxf = (CompositeXForm *) topXForm;
      for (i=0; i<cxf->mpLayer[0].mNBlocks; i++) 
      {
        UpdateXFormStatCache(xfsc, cxf->mpLayer[0].mpBlock[i], input);
        input += cxf->mpLayer[0].mpBlock[i]->mInSize;
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  UpdateXFormInstanceStatCaches(XFormInstance *xformInstance,
                                FLOAT *observation, int time)
  {
    int i, j;
    FLOAT *obs;
  
    if (xformInstance == NULL || xformInstance->mStatCacheTime == time) return;
  
    xformInstance->mStatCacheTime = time;
  
    if (xformInstance->mNumberOfXFormStatCaches == 0) return;
  
    if (xformInstance->mpInput) {
      UpdateXFormInstanceStatCaches(xformInstance->mpInput, observation, time);
    }
  
    for (i = 0; i < xformInstance->mNumberOfXFormStatCaches; i++) {
      XFormStatCache *xfsc = &xformInstance->mpXFormStatCache[i];
  
      if (xfsc->mpUpperLevelStats != NULL &&
        xfsc->mpUpperLevelStats->mpStats == xfsc->mpStats) { //just link to upper level?
        xfsc->mNorm = xfsc->mpUpperLevelStats->mNorm;
        continue;
      }
  
      obs = XFormPass(xformInstance->mpInput, observation, time, FORWARD);
      xfsc->mNorm = 0;
      UpdateXFormStatCache(xfsc, xformInstance->mpXForm, obs);
      if (xfsc->mpUpperLevelStats != NULL) {
        int size = xfsc->mpXForm->mInSize;
        for (j = 0; j < size + size*(size+1)/2; j++) {
          xfsc->mpStats[j] += xfsc->mpUpperLevelStats->mpStats[j];
        }
        xfsc->mNorm += xfsc->mpUpperLevelStats->mNorm;
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ReestState(Network *net, Node *node,
             int state_idx, FLOAT logPriorProb, FLOAT updateDir,
             FLOAT *obs, FLOAT *obs2) 
  {
    int i, j, k, m;
    State *state  = node->mpHmm->        mpState[state_idx];
    State *state2 = node->mpHmmToUpdate->mpState[state_idx];
    FLOAT bjtO    = LOG_0;
    int nmixtures;
  
    if (!net->mpModelSetToUpdate->mGaussLvl2ModelReest
    && net->mpModelSet != net->mpModelSetToUpdate) {
      // Occupation probabilities of mixtures are computed using target model
      state = state2; obs = obs2;
    } else if (state->mNumberOfMixtures <= state2->mNumberOfMixtures) {
      bjtO = net->mpOutPCache[net->mpModelSet->mNStates * net->mTime + state->mID].mValue;
    }
  
    nmixtures = LOWER_OF(state->mNumberOfMixtures, state2->mNumberOfMixtures);
  
    if (bjtO < LOG_MIN) {
      // State likelihood was not available in cache because
      // - occupation probabilities of mixtures are computed using target model
      // - not all mixtures of mAlignment model are used for computation of
      //   occupation probabilities (state2->num_mix < state->num_mix)
      for (m = 0; m < nmixtures; m++) {
        Mixture *mix = state->mpMixture[m].mpEstimates;
        FLOAT cjm    = state->mpMixture[m].mWeight;
        FLOAT *xobs  = XFormPass(mix->mpInputXForm, obs, net->mTime, FORWARD);
        FLOAT bjmtO  = DiagCGaussianDensity(mix, xobs, net);
        bjtO  = LogAdd(bjtO, cjm + bjmtO);
      }
    }
  
    for (m = 0; m < nmixtures; m++) {                  //for every emitting state
      Mixture *mix = state->mpMixture[m].mpEstimates;
      FLOAT cjm    = state->mpMixture[m].mWeight;
      FLOAT *xobs  = XFormPass(mix->mpInputXForm, obs, net->mTime, FORWARD);
      FLOAT bjmtO  = DiagCGaussianDensity(mix, xobs, net);
      FLOAT Lqjmt  = logPriorProb - bjtO + cjm + bjmtO;
  
      if (Lqjmt > MIN_LOG_WEGIHT) {
        Mixture *mix       = state2->mpMixture[m].mpEstimates;
        int vec_size       = mix->mpMean->mVectorSize;
        FLOAT *mnvec       = mix->mpMean->mpVectorO;
        FLOAT *mnacc       = mix->mpMean->mpVectorO     +     vec_size;
        FLOAT *vvacc       = mix->mpVariance->mpVectorO +     vec_size;
        FLOAT *vmacc       = mix->mpVariance->mpVectorO + 2 * vec_size;
        XFormInstance *ixf = mix->mpInputXForm;
        FLOAT *xobs        = XFormPass(ixf, obs2, net->mTime, FORWARD);
  
  /*      if (net->mmi_den_pass) {
          mnacc += vec_size + 1;
          vvacc += 2 * vec_size + 1;
          vmacc += 2 * vec_size + 1;
        }*/
  
        Lqjmt = exp(Lqjmt) * updateDir;
  
        for (i = 0; i < vec_size; i++) {                 // Update
          mnacc[i] += Lqjmt * xobs[i];                  // mean
          if (net->mpModelSetToUpdate->mUpdateMask & UM_OLDMEANVAR) {
            vvacc[i] += Lqjmt * SQR(xobs[i]-mnvec[i]);  // var
          } else {
            vvacc[i] += Lqjmt * SQR(xobs[i]);           // scatter
            vmacc[i] += Lqjmt * xobs[i];                //mean for var.
          }
        }
  
        mnacc[vec_size] += Lqjmt; //norms for mean
        vmacc[vec_size] += Lqjmt; //norms for variance
  
  //      if (net->mmi_den_pass) {
  //        state2->mpMixture[m].mWeightAccumDen += Lqjmt;
  //      } else {
        state2->mpMixture[m].mWeightAccum     += Lqjmt; //  Update weight accum
  //      }
  
        if (Lqjmt < 0) {
          state2->mpMixture[m].mWeightAccumDen -= Lqjmt;
        }
  
        if (ixf == NULL || ixf->mNumberOfXFormStatCaches == 0) continue;
  
        UpdateXFormInstanceStatCaches(ixf, obs2, net->mTime);
  
        for (i = 0; i < ixf->mNumberOfXFormStatCaches; i++) {
          XFormStatCache *xfsc = &ixf->mpXFormStatCache[i];
          Variance *var  = state2->mpMixture[m].mpEstimates->mpVariance;
          Mean     *mean = state2->mpMixture[m].mpEstimates->mpMean;
  
          for (j = 0; j < var->mNumberOfXFormStatAccums; j++) {
            XFormStatAccum *xfsa = &var->mpXFormStatAccum[j];
            if (xfsa->mpXForm == xfsc->mpXForm) {
              int size = xfsc->mpXForm->mInSize;
              for (k = 0; k < size+size*(size+1)/2; k++) {
                xfsa->mpStats[k] += xfsc->mpStats[k] * Lqjmt;
              }
              xfsa->mNorm += xfsc->mNorm * Lqjmt;
              break;
            }
          }
  
          for (j = 0; j < mean->mNumberOfXFormStatAccums; j++) {
            XFormStatAccum *xfsa = &mean->mpXFormStatAccum[j];
            if (xfsa->mpXForm == xfsc->mpXForm) {
              int size = xfsc->mpXForm->mInSize;
              for (k = 0; k < size; k++) {
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
  StateOccupationProbability(Network *net, FLOAT *obsMx, ModelSet *hmms,
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
      if (net->mpNodes[i].type & NT_Model) nNetModels++;
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
          Node *node = &net->mpNodes[k];
          if (node->mType & NT_Model) {
            for (j = 0; j < node->mpHmm->mNStates - 2; j++, state_counter++) {
              FLOAT tmpf = net->OutputProbability(node->mpHmm->mpState[j],
                                        obsMx + hmms->mInputVectorSize * i, net);
              switch (getMahalDist) {
                case 0:
                  break;
                case 1:
                tmpf = log((tmpf / -0.5) - node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->mGConst) * -0.5;
                break;
                case 2:
                tmpf = log((tmpf / -0.5) - node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->mGConst) * -1;
                break;
                case 3:
                tmpf += node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->mGConst * 0.5;
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
  
    if (!IS_ACTIVE(*net->mpFirst->exitToken)) {
      TokenPropagationDone(net);
      net->mpOutPCache = outPCache;
      free(alfa);
      for (i=0; i < net->mNumberOfNetStates * nFrames; i++) {
        beta[i] = LOG_0;
      }
      return beta; // No token survivered
    }
  
    totalLike = net->mpFirst->exitToken->mLike;
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
        Node *node = &net->mpNodes[k];
        if (node->mType & NT_Model) {
          for (j = 0; j < node->mpHmm->mNStates - 2; j++, state_counter++) {
            int idx = (net->mNumberOfNetStates * i) + node->estate_id + j + 1;
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
  WLR *
  TimePruning(Network *net, int frame_delay)
  {
    size_t    i;
    Node *    node;
    Token *   token = net->mpBestToken;
    WLR *     twlr;
    WLR *     rwlr=NULL;
  
    if (frame_delay > net->mTime-1 || !token) 
      return NULL;
  
    if (token->twlr != NULL &&
      token->twlr->mTime == net->mTime-1 - frame_delay) 
    {
      rwlr = token->twlr;
    }
  
    for (node = net->mpFirst; node != NULL; node = node->mpNext) 
    {
      if (!(node->mType & NT_Model)) 
        continue;
  
      for (i = 0; i < node->mpHmm->mNStates-1; i++) 
      {
        if (node->tokens[i].IsActive()) 
        {
          Token *token = &node->tokens[i];
  
          if (token->twlr != NULL &&
            token->twlr->mTime == net->mTime-1 - frame_delay) 
          {
            if (rwlr != token->twlr) 
            {
              KillToken(token);
            } 
            else 
            {
              if (token->wlr == token->twlr) 
              {
                token->twlr = NULL;
              } 
              else 
              {
                for (twlr = token->wlr; twlr->mpNext != token->twlr; twlr = twlr->mpNext)
                  ;
                token->twlr = twlr;
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
  Network:: 
  Init(Node * pFirstNode, ModelSet * pHmms, ModelSet *pHmmsToUpdate) 
  {
    Node *node;
    int maxStatesInModel = 0;
    int i;
  
    mpFirst = mpLast = pFirstNode;
  
    PhoneNodesToModelNodes(pFirstNode, pHmms, pHmmsToUpdate);
  
    //Allocate tokens and count emiting states
    mNumberOfNetStates = 0;
    
    for (node = pFirstNode; node != NULL; mpLast = node, node = node->mpNext) 
    {
      
  #ifndef NDEBUG
      node->aux2 = 0;
  #endif
      int numOfTokens = 1;
      if (node->mType & NT_Model) 
      {
        numOfTokens = node->mpHmm->mNStates;
        node->mpHmmToUpdate->mpMacro->mOccurances++;
        if (node->mpHmm->mpTransition->mpMatrixO[numOfTokens - 1] > LOG_MIN) {
          node->mType |= NT_Tee;
        }
      } 
      else if (node->mType & NT) 
      {
        numOfTokens = 1;
      } 
      else 
      {
        Error("Fatal: Incorect node type");
      }
  
      node->tokens = (Token *) malloc(numOfTokens * sizeof(Token));
      
      if (node->tokens == NULL) 
        Error("Insufficient memory");
  
      node->exitToken = &node->tokens[numOfTokens-1];
      node->estate_id = mNumberOfNetStates;
      
      if (node->mType & NT_Model) {
        int nstates = node->mpHmm->mNStates;
  
        if (maxStatesInModel < nstates) maxStatesInModel = nstates;
        assert(nstates >= 2); // two non-emiting states
        mNumberOfNetStates += nstates;
      }
    }
  
    SortNodes(this);  
    
    mpAuxTokens = (Token *) malloc((maxStatesInModel-1) * sizeof(Token));
    mpOutPCache = (Cache *) malloc(pHmms->mNStates   * sizeof(Cache));
    mpMixPCache = (Cache *) malloc(pHmms->mNMixtures * sizeof(Cache));
  
    if (mpAuxTokens == NULL ||
      mpOutPCache == NULL || mpMixPCache == NULL) {
      Error("Insufficient memory");
    }
  
    for (i = 0; i < maxStatesInModel-1; i++) {
      mpAuxTokens[i].mLike = LOG_0;
      mpAuxTokens[i].wlr = NULL;
    }
  
    mWPenalty          = 0.0;
    mMPenalty          = 0.0;
    mPronScale         = 1.0;
    mTranScale         = 1.0;
    mOutpScale         = 1.0;
    mOcpScale          = 1.0;
    mLmScale           = 1.0;
    OutputProbability =
      pHmms->mOutPdfKind == KID_DiagC     ? &DiagCGaussianMixtureDensity :
      pHmms->mOutPdfKind == KID_PDFObsVec ? &FromObservationAtStateId    : NULL;
  
    PassTokenInNetwork= &PassTokenMax;
    PassTokenInModel  = &PassTokenMax;
    mPropagDir         = FORWARD;
    mAlignment         = WORD_ALIGNMENT;
    mpThreshState       = NULL;
    mPruningThresh     = -LOG_0;
    mpModelSet            = pHmms;
    mpModelSetToUpdate    = pHmmsToUpdate;
    mCollectAlphaBeta  = 0;
  //  mmi_den_pass      = 0;
    mAccumType          = AT_ML;            
    mSearchPaths        = SP_ALL;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  Network::
  Release()
  {
    Node *node;
  
    for (node = mpFirst; node != NULL; node = node->mpNext) 
      free(node->tokens);
      
    FreeNetwork(mpFirst);
    
    free(mpAuxTokens);
    free(mpOutPCache);
    free(mpMixPCache);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  // Token section
  //***************************************************************************
  //***************************************************************************
  
  //***************************************************************************
  //***************************************************************************
  void 
  Token::
  Kill()
  {
    this->mLike = LOG_0;
    FreeWordLinkRecords(this->wlr);
    this->wlr  = NULL;
    this->twlr = NULL;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Label *
  Token::
  pGetLabels()
  {
    WLR *     wlr;
    Label *   tmp;
    Label *   level[3] = {NULL, NULL, NULL};
    int       li = 0;
      
    if (!this->IsActive())
      return NULL;
    
    for (wlr = this->wlr; wlr != NULL; wlr = wlr->mpNext) 
    {
      if ((tmp = (Label *) malloc(sizeof(Label))) == NULL)
        Error("Insufficient memory");
  
      tmp->mScore = wlr->mLike;
      tmp->mStop  = wlr->mTime;
      tmp->mId    = wlr->mStateIdx;
  
      if (wlr->mpNode->mType & NT_Model) 
      {
        li = wlr->mStateIdx >= 0 ? 0 : 1;
        tmp->mpData = wlr->mpNode->mpHmm;
        tmp->mpName = wlr->mpNode->mpHmm->mpMacro->mpName;
        tmp->mpNextLevel = level[li+1] ? level[li+1] : level[2];
      } 
      else //if (wlr->mpNode->mpPronun->outSymbol) 
      {
        li = 2;
        tmp->mpData = wlr->mpNode->mpPronun->word;
        tmp->mpName = wlr->mpNode->mpPronun->outSymbol;
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
            level[li-1]->mStart  = tmp2->mStop;
          }
  
          tmp2->mpNext  = level[li-1];
          level[li-1] = tmp2;
        }*/
  
        level[li]->mScore -= tmp->mScore;
        level[li]->mStart = tmp->mStop;
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
  
