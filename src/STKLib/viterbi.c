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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include "viterbi.h"
#include "common.h"
#include "labels.h"

#ifdef MOTIF
#include "imagesc.h"
#endif

//#define TRACE_TOKENS
#define SQR(x) ((x) * (x))
#define IS_ACTIVE(token) ((token).like > LOG_MIN)
#define FORWARD_PASS (net->propagDir == FORWARD)

using namespace STK;
void AddWordLinkRecord(Token *token, Node *node, int state_idx, int time);
void TokenPropagationInNetwork(Network *net);
void ReestState(Network *net, Node *node, int state_idx,
                FLOAT logPriorProb, FLOAT updateDir, FLOAT *obs, FLOAT *obs2);
int BackwardPruning(int time, Node *node, int state);

#ifdef DEBUG_MSGS
WLR *firstWLR;
#endif

void SortNodes(Network *net);

void PhoneNodesToModelNodes(Node *first, ModelSet *hmms, ModelSet *hmmsToUpdate)
{
  Node *node;

  if (hmmsToUpdate == NULL) hmmsToUpdate = hmms;

  for (node = first; node != NULL; node = node->next) {
    if (node->mType & NT_Phone) {
      Macro *macro;

      node->mType &= ~NT_Phone;
      node->mType |= NT_Model;
      macro = FindMacro(&hmms->mHmmHash, node->mpName);
      if (macro == NULL) {
        Error("Model %s not defined in %sHMM set", node->mpName,
              hmmsToUpdate != hmms ? "alignment " : "");
      }
      node->hmm = node->hmmToUpdate = (Hmm *) macro->mpData;

      if (hmmsToUpdate != hmms) {
        macro = FindMacro(&hmmsToUpdate->mHmmHash, node->mpName);
        if (macro == NULL) {
          Error("Model %s not defined in HMM set", node->mpName,
                hmmsToUpdate != hmms ? "" : "target ");
        }
        node->hmmToUpdate = (Hmm *) macro->mpData;
      }
    }
  }
}

void InitNetwork(Network *net, Node *first, ModelSet *hmms, ModelSet *hmmsToUpdate) {
  Node *node;
  int maxStatesInModel = 0;
  int i;

  net->first = net->last  = first;

  PhoneNodesToModelNodes(first, hmms, hmmsToUpdate);

  //Allocate tokens and count emiting states
  net->nNetStates = 0;
  for (node = first; node != NULL; net->last = node, node = node->next) {
#ifndef NDEBUG
    node->aux2 = 0;
#endif
    int numOfTokens = 1;
    if (node->mType & NT_Model) {
      numOfTokens = node->hmm->mNStates;
      node->hmmToUpdate->mpMacro->mOccurances++;
      if (node->hmm->mpTransition->matrix[numOfTokens - 1] > LOG_MIN) {
        node->mType |= NT_Tee;
      }
    } else if (node->mType & NT) {
      numOfTokens = 1;
    } else Error("Fatal: Incorect node type");

    node->tokens = (Token *) malloc(numOfTokens * sizeof(Token));
    if (node->tokens == NULL) Error("Insufficient memory");

    node->exitToken = &node->tokens[numOfTokens-1];
    node->estate_id = net->nNetStates;

    if (node->mType & NT_Model) {
      int nstates = node->hmm->mNStates;

      if (maxStatesInModel < nstates) maxStatesInModel = nstates;
      assert(nstates >= 2); // two non-emiting states
      net->nNetStates += nstates;
    }
  }

  SortNodes(net);

  net->auxTokens = (Token *) malloc((maxStatesInModel-1) * sizeof(Token));
  net->outPCache = (Cache *) malloc(hmms->mNStates   * sizeof(Cache));
  net->mixPCache = (Cache *) malloc(hmms->mNStates * sizeof(Cache));

  if (net->auxTokens == NULL ||
     net->outPCache == NULL || net->mixPCache == NULL) {
    Error("Insufficient memory");
  }

  for (i = 0; i < maxStatesInModel-1; i++) {
    net->auxTokens[i].like = LOG_0;
    net->auxTokens[i].wlr = NULL;
  }

  net->wPenalty          = 0.0;
  net->mPenalty          = 0.0;
  net->pronScale         = 1.0;
  net->tranScale         = 1.0;
  net->outpScale         = 1.0;
  net->ocpScale          = 1.0;
  net->lmScale           = 1.0;
  net->OutputProbability =
    hmms->mOutPdfKind == KID_DiagC     ? &DiagCGaussianMixtureDensity :
    hmms->mOutPdfKind == KID_PDFObsVec ? &FromObservationAtStateId    : NULL;

  net->PassTokenInNetwork= &PassTokenMax;
  net->PassTokenInModel  = &PassTokenMax;
  net->propagDir         = FORWARD;
  net->alignment         = WORD_ALIGNMENT;
  net->threshState       = NULL;
  net->pruningThresh     = -LOG_0;
  net->hmmSet            = hmms;
  net->hmmSetToUpdate    = hmmsToUpdate;
  net->collectAlphaBeta  = 0;
//  net->mmi_den_pass      = 0;
  net->accumType          = AT_ML;
}

void ReleaseNetwork(Network *net)
{
  Node *node;

  for (node = net->first; node != NULL; node = node->next) free(node->tokens);
  FreeNetwork(net->first);
  free(net->auxTokens);
  free(net->outPCache);
  free(net->mixPCache);
}

int cmplnk(const void *a, const void *b)
{
  return ((Link *) a)->node->aux - ((Link *) b)->node->aux;
}

void SortNodes(Network *net)
{
  int i, j;
  Node *chain, *last, *node;

// Sort nodes for forward (Viterbi) propagation

  for (node = net->first; node != NULL; node = node->next) {
    node->aux = node->nbacklinks;
  }

  for (i = 0; i < net->first->nlinks; i++) {
    net->first->links[i].node->aux--;
  }

  last = net->first;
  chain = net->first->next;

  while (chain) {
    BOOL short_curcuit = TRUE;
    Node **curPtr = &chain;
    i = 0;

    while (*curPtr) {
      if ((((*curPtr)->mType & NT_Model) && !((*curPtr)->mType & NT_Tee))
        || (*curPtr)->aux == 0) {
        for (j = 0; j < (*curPtr)->nlinks; j++) {
          (*curPtr)->links[j].node->aux--;
        }

        last = (last->next = *curPtr);
        last->aux = i++;
        *curPtr = (*curPtr)->next;
        short_curcuit = FALSE;
      } else {
        curPtr = &(*curPtr)->next;
      }
    }

    if (short_curcuit) {
//      fprintf(stderr, "Nodes in loop: ");
//      for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->next)
//        fprintf(stderr, "%d %d", *curPtr - net->nodes, (*curPtr)->mType);
//      fprintf(stderr, "\n");
      Error("Loop of non-emiting nodes found in network");
    }
  }

  last->next = NULL;

  /// !!! What is this sorting links good for ???
  for (node = net->first; node != NULL; node = node->next) {
      if (node->nlinks > 1)
      qsort(node->links, node->nlinks, sizeof(Link), cmplnk);
  }

// Sort nodes for backward propagation

  for (node = net->first; node != NULL; node = node->next) {
    node->aux = node->nlinks;
  }

  for (i = 0; i < net->last->nbacklinks; i++) {
    net->last->backlinks[i].node->aux--;
  }

  last = net->last;
  chain = net->last->backnext;
  i = 0;

  while (chain) {
    BOOL short_curcuit = TRUE;
    Node **curPtr = &chain;

    while (*curPtr) {
      if ((((*curPtr)->mType & NT_Model) && !((*curPtr)->mType & NT_Tee))
        || (*curPtr)->aux == 0) {
        for (j = 0; j < (*curPtr)->nbacklinks; j++) {
          (*curPtr)->backlinks[j].node->aux--;
        }

        last = (last->backnext = *curPtr);
        last->aux = i++;
        *curPtr = (*curPtr)->backnext;
        short_curcuit = FALSE;
      } else {
        curPtr = &(*curPtr)->backnext;
      }
    }

    /*if (short_curcuit) {
      fprintf(stderr, "Nodes in loop: ");
      for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->backnext)
        fprintf(stderr, "%d ", *curPtr - net->nodes);
      fprintf(stderr, "\n");
      Error("Loop of non-emiting nodes found in net");
    }*/

    assert(!short_curcuit); // Shouldn't happen, since it didnot happen before
  }

  last->backnext = NULL;

  /// !!! What is this sorting links good for ???
  for (node = net->first; node != NULL; node = node->next) {
      if (node->nbacklinks > 1)
      qsort(node->backlinks, node->nbacklinks, sizeof(Link), cmplnk);
  }
}

void WriteAlpha(int time, Node *node, int state, Token *token)
{
  if (node->alphaBetaListReverse == NULL ||
     node->alphaBetaListReverse->time != time) 
  {
    size_t    i;
    FWBWR *   newrec;
    
    newrec  = (FWBWR*) malloc(sizeof(FWBWR) +
                              sizeof(newrec->state[0]) * (node->hmm->mNStates-1));
    if (newrec == NULL) 
      Error("Insufficient memory");
      
    newrec->next = node->alphaBetaListReverse;
    newrec->time = time;
    
    for (i=0; i<node->hmm->mNStates; i++) 
    {
      newrec->state[i].alpha = newrec->state[i].beta = LOG_0;
    }
    
    node->alphaBetaListReverse = newrec;
  }
  node->alphaBetaListReverse->state[state].alpha = token->like;
  node->alphaBetaListReverse->state[state].alphaAccuracy = token->accuracy;
}

void WriteBeta(int time, Node *node, int state, Token *token)
{

  // Record for current time must be already moved to
  // alphaBetaList by function BackwardPruning
  assert(node->alphaBetaListReverse == NULL ||
         node->alphaBetaListReverse->time < time);

  if (node->alphaBetaList != NULL && node->alphaBetaList->time == time) {
    node->alphaBetaList->state[state].beta = token->like;
    node->alphaBetaList->state[state].betaAccuracy = token->accuracy;
  }
}

int BackwardPruning(int time, Node *node, int state)
{
  while (node->alphaBetaListReverse != NULL &&
        node->alphaBetaListReverse->time > time) {
    FWBWR *fwbwr = node->alphaBetaListReverse;
    node->alphaBetaListReverse = fwbwr->next;
    free(fwbwr);
  }

  if (node->alphaBetaListReverse != NULL &&
     node->alphaBetaListReverse->time == time) {
    FWBWR *fwbwr = node->alphaBetaListReverse;
    node->alphaBetaListReverse = fwbwr->next;
    fwbwr->next = node->alphaBetaList;
    node->alphaBetaList = fwbwr;
  }

  return !(node->alphaBetaList != NULL &&
           node->alphaBetaList->time == time &&
           node->alphaBetaList->state[state].alpha > LOG_MIN);
}

void FreeFWBWRecords(Network *net)
{
  Node *node;
  for (node = net->first; node != NULL; node = node->next) {
    if (!(node->mType & NT_Model)) continue;

    while (node->alphaBetaList) {
      FWBWR *fwbwr = node->alphaBetaList;
      node->alphaBetaList = fwbwr->next;
      free(fwbwr);
    }

    while (node->alphaBetaListReverse) {
      FWBWR *fwbwr = node->alphaBetaListReverse;
      node->alphaBetaListReverse = fwbwr->next;
      free(fwbwr);
    }
  }
}


#ifndef NDEBUG
int test_for_cycle = 0;
int HasCycleCounter = 100000;

int HasCycle(Network *net) {
  Node *node;
  HasCycleCounter++;
  if (!test_for_cycle) return 0;
  for (node = net->activeNodes; node; node= node->nextActiveNode) {
    int i, nlinks = FORWARD_PASS ? node->nlinks : node->nbacklinks;
    Link *links   = FORWARD_PASS ? node->links  : node->backlinks;
    if (node->aux2 == HasCycleCounter) {
      printf("Cycle in list of active nodes\n");
      return 1;
    }
    node->aux2 = HasCycleCounter;

    for (i=0; i <nlinks; i++)
      if (links[i].node->aux2 == HasCycleCounter &&
         (!(links[i].node->mType & NT_Model)
         || links[i].node->mType & NT_Tee)) {
        printf("Active node %d listed after his non-model succesor %d\n",
               node - net->first, links[i].node - net->first);

        return 2;
      }
  }
  return 0;
}

int AllWordSuccessorsAreActive(Network *net) {
  Node *node;
  if (!test_for_cycle) return 1;

  for (node = net->activeNodes; node; node= node->nextActiveNode) {
    int i, nlinks = FORWARD_PASS ? node->nlinks : node->nbacklinks;
    Link *links   = FORWARD_PASS ? node->links  : node->backlinks;

    for (i=0; i <nlinks; i++)

      if (links[i].node->aux2 != HasCycleCounter &&
         links[i].node != (FORWARD_PASS ? net->last : net->first) &&
         (!(links[i].node->mType & NT_Model)
         || links[i].node->mType & NT_Tee)) {
        printf("Active node %d has nonactive non-model succesor %d\n",
               node - net->first, links[i].node - net->first);
        return 0;
      }
  }
  return 1;
}
#endif

void MarkWordNodesLeadingFrom(Network *net, Node *node)
{
  int i, nlinks = FORWARD_PASS ? node->nlinks : node->nbacklinks;
  Link *links   = FORWARD_PASS ? node->links  : node->backlinks;

  for (i = 0; i < nlinks; i++) {
    Node *lnode = links[i].node;
    if ((lnode->mType & NT_Model && !(lnode->mType & NT_Tee))
       || (lnode == (FORWARD_PASS ? net->last : net->first))) continue;
    if (lnode->isActiveNode > 0) continue;

    if (lnode->isActive) {
      assert(lnode->mType & NT_Tee);
      continue;
    }
    
    if(lnode->isActiveNode-- == 0) {
      lnode->aux = 0;
      MarkWordNodesLeadingFrom(net, lnode);
    }
  }
}


Node *ActivateWordNodesLeadingFrom(Network *net, Node *node)
{
  int i, nlinks = FORWARD_PASS ? node->nlinks : node->nbacklinks;
  Link *links   = FORWARD_PASS ? node->links  : node->backlinks;

  for (i = 0; i < nlinks; i++) {
    Node *lnode = links[i].node;
    if ((lnode->mType & NT_Model && !(lnode->mType & NT_Tee))
       || (lnode == (FORWARD_PASS ? net->last : net->first))) continue;
    if (lnode->isActiveNode++ > 0) continue;

    if (lnode->isActive) {
      assert(lnode->mType & NT_Tee);
      continue;
    }

    lnode->aux++;
    if (lnode->isActiveNode < 0) continue;

    assert(lnode->isActiveNode == 0);
    lnode->isActiveNode = lnode->aux;

    lnode->nextActiveNode = node->nextActiveNode;
    lnode->prevActiveNode = node;
    if (node->nextActiveNode) node->nextActiveNode->prevActiveNode = lnode;
    node->nextActiveNode  = lnode;
    node = ActivateWordNodesLeadingFrom(net, lnode);
  }
  assert(!HasCycle(net));
  return node;
}

void ActivateModel(Network *net, Node *node)
{
  if (node->isActive) return;
  node->isActive = 1;
  node->prevActiveModel = NULL;
  node->nextActiveModel = net->activeModels;
  if (net->activeModels != NULL) {
    net->activeModels->prevActiveModel = node;
  }
  net->activeModels = node;

  if (node->isActiveNode) {
    assert(node->mType & NT_Tee);
    return;
  }

  node->isActiveNode = 1;
  node->prevActiveNode = NULL;
  node->nextActiveNode = net->activeNodes;
  if (net->activeNodes != NULL) {
    net->activeNodes->prevActiveNode = node;
  }
  net->activeNodes = node;

  assert(!HasCycle(net));

  MarkWordNodesLeadingFrom(net, node);
  ActivateWordNodesLeadingFrom(net, node);

  assert(AllWordSuccessorsAreActive(net));
}

void DeactivateWordNodesLeadingFrom(Network *net, Node *node)
{
  int i, nlinks = FORWARD_PASS ? node->nlinks : node->nbacklinks;
  Link *links   = FORWARD_PASS ? node->links  : node->backlinks;

  for (i = 0; i < nlinks; i++) {
    Node *lnode = links[i].node;
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

void DeactivateModel(Network *net, Node *node)
{
  if (!node->isActive) return;
  node->isActive = 0;

  if (node->nextActiveModel != NULL) {
    node->nextActiveModel->prevActiveModel = node->prevActiveModel;
  }

  if (node->prevActiveModel != NULL) {
    node->prevActiveModel->nextActiveModel = node->nextActiveModel;
  } else {
    assert(net->activeModels == node);
    net->activeModels = node->nextActiveModel;
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
    assert(net->activeNodes == node);
    net->activeNodes = node->nextActiveNode;
  }

  assert(!HasCycle(net));
}

void TokenPropagationInit(Network *net)
{
  int i, j;
  Node *node;

  InitLogMath();

//  for (i=0; i < net->nnodes; i++) {
//  node = &net->nodes[i];
  for (node = net->first; node != NULL; node = node->next) {
    int numOfTokens = (node->mType & NT_Model) ? node->hmm->mNStates : 1;

    for (j=0; j < numOfTokens; j++) {
      node->tokens[j].like = LOG_0;
      node->tokens[j].wlr = NULL;
#ifdef bordel_staff
      node->tokens[j].twlr = NULL;
#endif // bordel_staff
    }
    node->isActive = 0;
    node->isActiveNode = 0;
  }

  if (net->outPCache != NULL) {
    for (i = 0; i < net->hmmSet->mNStates; i++) {
      net->outPCache[i].time = UNDEF_TIME;
    }
  }

  if (net->mixPCache != NULL) {
    for (i = 0; i < net->hmmSet->mNStates; i++) {
      net->mixPCache[i].time = UNDEF_TIME;
    }
  }

  net->bestToken  = NULL;
  net->wordThresh = LOG_MIN;
  net->activeModels = NULL;
  net->activeNodes = NULL;
  net->activeTokens = 0;
//  net->time = 0;
  node = FORWARD_PASS ? net->first : net->last;
  node->tokens[0].like = 0;
  node->tokens[0].accuracy.logvalue = LOG_0;
  node->tokens[0].accuracy.negative = 0;
  node->tokens[0].wlr = NULL;
#ifdef bordel_staff
  node->tokens[0].twlr = NULL;
  node->tokens[0].bestlike = 0;
#endif

  if (net->collectAlphaBeta && FORWARD_PASS) {
    for (node = net->first; node != NULL; node = node->next) {
      if (!(node->mType & NT_Model)) continue;
      node->alphaBetaList = node->alphaBetaListReverse = NULL;
    }
  }
  // Needed to load last FWBWRs to alphaBetaList
  if (net->collectAlphaBeta && !FORWARD_PASS) {
    for (node = net->first; node != NULL; node = node->next) {
      if (!(node->mType & NT_Model)) continue;
      BackwardPruning(net->time, node, node->hmm->mNStates-1);
    }
  }


  node = FORWARD_PASS ? net->first : net->last;
  net->activeNodes = node;
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

void TokenPropagationInNetwork(Network *net)
{
  Node *node;
  Link *links;
  int i, nlinks;

  // Beam pruning is not active in backward pass. First, it is not necessary
  // since only token that fit into the forward pass beam are allowed (backward
  // pruning after forward pass). Second, it could be dangerous. For example,
  // when training from lattices where LM weights are pushed to the begining,
  // the lowest weight (the biggest penalty) will be on the wery first link.
  // If this weight was lower than minus pruning treshold, it would result in
  // always killing token in the first node, when pruning during the backward pass.
  //                                      |
  net->beamThresh = FORWARD_PASS && // <--'
                    net->bestToken &&
                    net->bestToken->like - net->pruningThresh > LOG_MIN
                    ? net->bestToken->like - net->pruningThresh : LOG_MIN;

  node = FORWARD_PASS ? net->last : net->first;
  KillToken(node->exitToken);

//  Node *Xnode = net->activeNodes;
/*  for (node = FORWARD_PASS ? net->first : net->last;
      node != NULL;
      node = FORWARD_PASS ? node->next : node->backnext) { //*/
  for (node = net->activeNodes; node != NULL; node = node->nextActiveNode) {
    if ((node->mType & NT_Tee) && IS_ACTIVE(node->tokens[0])) {
//        assert(node->isActiveNode || (node->mType & NT_Tee && node->isActive));
//        for (Xnode = net->activeNodes; Xnode && Xnode != node; Xnode = Xnode->nextActiveNode);
//        assert(Xnode);

        if (!(net->collectAlphaBeta && !FORWARD_PASS && // backward pruning
          BackwardPruning(net->time, node, 0))) {   // after forward pass
        Hmm *hmm = node->hmm;
        FLOAT transP = hmm->mpTransition->matrix[hmm->mNStates - 1];
  #ifdef TRACE_TOKENS
        printf("Tee model State 0 -> Exit State ");
  #endif
        net->PassTokenInModel(&node->tokens[0], node->exitToken,
                                        transP * net->tranScale);
      }
    }

    if (IS_ACTIVE(*node->exitToken)) {
//      assert(node->isActiveNode || (node->mType & NT_Tee && node->isActive));
//      for (Xnode = net->activeNodes; Xnode && Xnode != node; Xnode = Xnode->nextActiveNode);
//      assert(Xnode);

      if (node->mType & NT_Model) {
        if (net->collectAlphaBeta) {
          if (FORWARD_PASS) {
            WriteAlpha(net->time, node, node->hmm->mNStates-1, node->exitToken);
          } else {
            WriteBeta(net->time, node, 0, node->exitToken);
          }
        }
        node->exitToken->like += net->mPenalty;
      } else if (node->mType & NT && node->pronun != NULL) {
        node->exitToken->like += net->wPenalty +
                                 net->pronScale * node->pronun->prob;
        /*if (node->exitToken->like < net->wordThresh) {
          node->exitToken->like = LOG_0;
        }*/
      }

      if (node->exitToken->like > net->beamThresh) {
        if (node->mType & NT && node->pronun != NULL &&
           net->alignment & WORD_ALIGNMENT) {
          AddWordLinkRecord(node->exitToken, node, -1, net->time);
        } else if (node->mType & NT_Model && net->alignment & MODEL_ALIGNMENT) {
          AddWordLinkRecord(node->exitToken, node, -1, net->time);
        }

        nlinks = FORWARD_PASS ? node->nlinks : node->nbacklinks;
        links  = FORWARD_PASS ? node->links  : node->backlinks;

        for (i = 0; i < nlinks; i++) {
          FLOAT lmLike = links[i].like * net->lmScale;
          if (node->exitToken->like + lmLike > net->beamThresh 
          && (/*links[i].node->start == UNDEF_TIME ||*/
                links[i].node->start <= net->time) 
    && (  links[i].node->stop  == UNDEF_TIME ||
                links[i].node->stop  >= net->time)
    && (net->SearchPaths != Network::SP_TrueOnly      || 
       (node->mType & NT_True)                || 
       !(links[i].node->mType & NT_Model))) {
#ifdef TRACE_TOKENS
          printf("Node %d -> Node %d ", node->aux, links[i].node->aux);
#endif
            net->PassTokenInNetwork(node->exitToken,
                                    &links[i].node->tokens[0], lmLike);
            if (links[i].node->mType & NT_Model) {
              ActivateModel(net, links[i].node);
            } else {
              assert(links[i].node->isActiveNode ||
                     links[i].node == (FORWARD_PASS ? net->last : net->first));
            }
          }
        }
      }
      if (!(node->mType & NT_Sticky)) KillToken(node->exitToken);
    }

  }


  if (net->collectAlphaBeta) {
    if (FORWARD_PASS) {
      for (node = net->activeModels; node != NULL; node = node->nextActiveModel) {
        if (/*!(node->mType & NT_Model) ||*/ node->tokens[0].like < LOG_MIN) continue;
        WriteAlpha(net->time, node, 0, &node->tokens[0]);
      }
    } else {
      for (node = net->activeModels; node != NULL; node = node->nextActiveModel) {
        if (/*!(node->mType & NT_Model) ||*/ node->tokens[0].like < LOG_MIN
          || BackwardPruning(net->time, node, node->hmm->mNStates-1)) continue;
        WriteBeta(net->time, node, node->hmm->mNStates-1, &node->tokens[0]);
      }
    }
  }

//  Go through newly activeted models and costruct list of active nodes
//  for (node=net->activeModels; node&&!node->isActive; node=node->nextActiveModel) {
//    ActivateNode(net, node);
//  }
  assert(!HasCycle(net));
}

void TokenPropagationInModels(Network *net, FLOAT *observation)
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

  net->bestToken = NULL;
/*  if (net->threshState) {
    threshOutProb = net->OutputProbability(net->threshState, observation, net);
  } */

  for (node = net->activeModels; node != NULL; node = node->nextActiveModel) {
//  for (node = net->first; node != NULL; node = node->next) {
//    if (!(node->mType & NT_Model)) continue;

    hmm = node->hmm;

    if ((/*node->start != UNDEF_TIME &&*/node->start >= net->time)
    || (  node->stop  != UNDEF_TIME &&  node->stop  <  net->time)
    || net->SearchPaths == Network::SP_TrueOnly && !(node->mType & NT_True)) {
      for (i = 0; i < hmm->mNStates-1; i++) {
        KillToken(&node->tokens[i]);
      }
      DeactivateModel(net, node);
      continue;
    }

    if (net->accumType == AT_MPE && FORWARD_PASS && IS_ACTIVE(node->tokens[0])) {
      FloatInLog fil_lmpa =
        {node->tokens[0].like + log(fabs(node->phoneAccuracy)),
         node->phoneAccuracy < 0};

      node->tokens[0].accuracy = FIL_Add(node->tokens[0].accuracy, fil_lmpa);
    }

    for (i = 0; i < hmm->mNStates-1; i++) {
      assert(!IS_ACTIVE(net->auxTokens[i]));
      if (IS_ACTIVE(node->tokens[i])) {
        net->auxTokens[i] = node->tokens[i];
        node->tokens[i].like = LOG_0;
        node->tokens[i].wlr  = NULL;
      }
    }

    int keepModelActive = FALSE;

    assert(!IS_ACTIVE(node->tokens[hmm->mNStates-1]));

    for (j = 1; j < hmm->mNStates-1; j++) {
      state_idx = (FORWARD_PASS ? j : hmm->mNStates-1 - j);

      if (net->collectAlphaBeta &&
        !FORWARD_PASS &&
        BackwardPruning(net->time, node, state_idx)) {
        continue; // backward pruning after forward pass
      }

      for (i = 0; i < hmm->mNStates-1; i++) {
        from = FORWARD_PASS ? i : hmm->mNStates-1 - j;
        to   = FORWARD_PASS ? j : hmm->mNStates-1 - i;

        assert(!IS_ACTIVE(net->auxTokens[i]) || node->isActive);

        if (hmm->mpTransition->matrix[from * hmm->mNStates + to] > LOG_MIN &&
          IS_ACTIVE(net->auxTokens[i])) {
          FLOAT transP = hmm->mpTransition->matrix[from * hmm->mNStates + to];

#ifdef TRACE_TOKENS
          printf("Model %d State %d -> State %d ",  node->aux, i, j);
#endif
          if (net->PassTokenInModel(&net->auxTokens[i], &node->tokens[j],
                                   transP * net->tranScale)) {
            winingToken = i;
          }
        }
      }

      // if (IS_ACTIVE(node->tokens[j])) {
      if (node->tokens[j].like > net->beamThresh) {
        FLOAT outProb = net->OutputProbability(hmm->state[state_idx-1],
                                                   observation, net);
        outProb *= net->outpScale;

        /*if (outProb < threshOutProb) {
          outProb = threshOutProb;
        }*/

        if (net->collectAlphaBeta && !FORWARD_PASS) {
          WriteBeta(net->time, node, state_idx, &node->tokens[j]);
        }

        if (net->accumType == AT_MFE && node->mType & NT_True) {
          FloatInLog fil_like = {node->tokens[j].like, 0};
          node->tokens[j].accuracy = FIL_Add(node->tokens[j].accuracy, fil_like);
        }

        node->tokens[j].accuracy.logvalue += outProb;
        node->tokens[j].like              += outProb;

        if (net->collectAlphaBeta && FORWARD_PASS) {
          WriteAlpha(net->time, node, state_idx, &node->tokens[j]);
        }

        if (net->alignment & STATE_ALIGNMENT && winingToken > 0 &&
           (winingToken != j || net->alignment & FRAME_ALIGNMENT)) {
           AddWordLinkRecord(&node->tokens[j], node,
                             (FORWARD_PASS ? winingToken
                                           : hmm->mNStates-1 - winingToken)-1,
                             net->time-1);
        }

        net->activeTokens++;
        keepModelActive = TRUE;
        assert(node->isActive);
      } else {
        assert(node->isActive || !IS_ACTIVE(node->tokens[j]));
        KillToken(&node->tokens[j]);
      }
    }

    for (i = 0; i < hmm->mNStates-1; i++) {
      KillToken(&net->auxTokens[i]);
    }

    if (!keepModelActive) DeactivateModel(net, node);

    state_idx = (FORWARD_PASS ? hmm->mNStates - 1 : 0);
    assert(!IS_ACTIVE(node->tokens[hmm->mNStates - 1]));

    if (!keepModelActive ||
       (net->collectAlphaBeta && !FORWARD_PASS &&
        BackwardPruning(net->time-1, node, state_idx))) {
      // backward pruning after forward pass
      continue;
      //KillToken(&node->tokens[hmm->mNStates - 1]);
    }

    for (i = 1; i < hmm->mNStates-1; i++) {
      from = FORWARD_PASS ? i : 0;
      to   = FORWARD_PASS ? hmm->mNStates-1 : hmm->mNStates-1 - i;

      if (IS_ACTIVE(node->tokens[i])) {
        if (!net->bestToken || net->bestToken->like < node->tokens[i].like) {
          net->bestToken = &node->tokens[i];
          net->bestNode  = node;
        }

        if (hmm->mpTransition->matrix[from * hmm->mNStates + to] > LOG_MIN) {
          FLOAT transP = hmm->mpTransition->matrix[from * hmm->mNStates + to];
#ifdef TRACE_TOKENS
          printf("Model %d State %d -> Exit State ",  node->aux, i);
#endif
          if (net->PassTokenInModel(&node->tokens[i],
                                  &node->tokens[hmm->mNStates - 1],
                                  transP * net->tranScale)) {
            winingToken = i;
          }
        }
      }
    }

    if (IS_ACTIVE(node->tokens[hmm->mNStates - 1])) {
      if (net->accumType == AT_MPE && !FORWARD_PASS) {
        FloatInLog fil_lmpa =
          {node->tokens[hmm->mNStates - 1].like + log(fabs(node->phoneAccuracy)),
           node->phoneAccuracy < 0};

        node->tokens[hmm->mNStates - 1].accuracy =
          FIL_Add(node->tokens[hmm->mNStates - 1].accuracy, fil_lmpa);
      }

  //    ActivateNode(net, node);
      if (net->alignment & STATE_ALIGNMENT) {
        AddWordLinkRecord(&node->tokens[hmm->mNStates - 1], node,
                          (FORWARD_PASS ? winingToken :
                                          hmm->mNStates-1 - winingToken-1)-1,
                          net->time);
      }
    }
  }
  assert(!HasCycle(net));
}

void TokenPropagationDone(Network *net)
{
  int j;
  Node *node;

  KillToken(FORWARD_PASS ? net->last->exitToken
                         : net->first->exitToken);

//  for (i=0; i < net->nnodes; i++) {
//    node = &net->nodes[i];
    for (node = net->first; node != NULL; node = node->next) {
    int numOfTokens;

    if (!(node->mType & NT_Model)) {
      assert(!IS_ACTIVE(*node->exitToken));
      continue;
    }

    numOfTokens = node->hmm->mNStates;

    for (j=0; j < numOfTokens; j++) {
      KillToken(&node->tokens[j]);
    }
  }
}

int PassTokenMax(Token *from, Token *to, FLOAT like)
{
  int ret = 0;
#ifdef TRACE_TOKENS
  printf("(%.2f + %.2f -> %.2f = ", from->like, like, to->like);
#endif
  if (!IS_ACTIVE(*to) || from->like + like > to->like) {
    KillToken(to);

    ret = 1;
    *to = *from;
    to->like += like;
    if (to->wlr) {
      to->wlr->refs++;
    }
  }
#ifdef TRACE_TOKENS
  printf("%.2f)\n", to->like);
#endif
  return ret;
}

/*int PassTokenSum(Token *from, Token *to, FLOAT like)
{
  double tl;
  int ret = 0;
#ifdef TRACE_TOKENS
   printf("(%.2f + %.2f -> %.2f = ", from->like, like, to->like);
#endif
  if (IS_ACTIVE(*to)) {
    tl = LogAdd(to->like, from->like + like);
  } else {
    tl = from->like + like;
  }

  if (!IS_ACTIVE(*to) || from->like + like > to->bestlike) {
    KillToken(to);

    ret = 1;
    *to = *from;
    to->bestlike = from->like + like;
    if (to->wlr) {
      to->wlr->refs++;
    }
  }

  to->like = tl;
#ifdef TRACE_TOKENS
  printf("%.2f)\n", to->like);
#endif
  return ret;
}*/

int PassTokenSum(Token *from, Token *to, FLOAT like)
{
  double tl;
  FloatInLog fe, fil_from_like = {like, 0};
  int ret = 0;
#ifdef TRACE_TOKENS
   printf("(%.2f + %.2f -> %.2f = ", from->like, like, to->like);
#endif
  if (IS_ACTIVE(*to)) {
    tl = LogAdd(to->like, from->like + like);
    fe = FIL_Add(to->accuracy, FIL_Mul(from->accuracy, fil_from_like));
  } else {
    tl = from->like + like;
    fe = FIL_Mul(from->accuracy, fil_from_like);
  }

  if (!IS_ACTIVE(*to) || from->like + like > to->bestlike) {
    KillToken(to);

    ret = 1;
    *to = *from;
    to->bestlike = from->like + like;
    if (to->wlr) {
      to->wlr->refs++;
    }
  }

  to->like   = tl;
  to->accuracy = fe;
#ifdef TRACE_TOKENS
  printf("%.2f)\n", to->like);
#endif
  return ret;
}

void FreeWordLinkRecords(WLR *wlr)
{
  if (wlr != NULL) {
    --wlr->refs;
    assert(wlr->refs >= 0);

    if (wlr->refs == 0) {
      FreeWordLinkRecords(wlr->next);
#ifdef DEBUG_MSGS
      assert(wlr->freed == 0);
      wlr->freed = 1;
#else
      free(wlr);
#endif
    }
  }
}

void KillToken(Token *token)
{
  token->like = LOG_0;
  FreeWordLinkRecords(token->wlr);
  token->wlr  = NULL;
  token->twlr = NULL;
}

void AddWordLinkRecord(Token *token, Node *node, int state_idx, int time)
{
  WLR *wlr;

  if ((wlr = (WLR *) malloc(sizeof(WLR))) == NULL) {
    Error("Insufficient memory");
  }

  wlr->state_idx = state_idx;
  wlr->like  = token->like;
  wlr->node  = node;
  wlr->time  = time;
  wlr->next  = token->wlr;
  wlr->refs  = 1;
  token->wlr  = wlr;
#ifdef bordel_staff
  if (!token->twlr) token->twlr = wlr;
#endif // bordel_staff
#ifdef DEBUG_MSGS
  wlr->tmpNext = firstWLR;
  firstWLR = wlr;
  wlr->freed = 0;
#endif
}

FLOAT DiagCGaussianDensity(Mixture *mix, FLOAT  *obs, Network *net)
{
  FLOAT like = 0.0;
  int j;

  if (net && net->mixPCache[mix->mID].time == net->time) {
    return net->mixPCache[mix->mID].value;
  }

  for (j = 0; j < mix->mpMean->mVectorSize; j++) {
    like += SQR(obs[j] - mix->mpMean->mVector[j]) * mix->mpVariance->mVector[j];
  }

  like = -0.5 * (mix->mGConst + like);

  if (net) {
    net->mixPCache[mix->mID].time  = net->time;
    net->mixPCache[mix->mID].value = like;
  }

  return like;
}

#ifdef DEBUG_MSGS
int gaus_computaions = 0;
#endif

FLOAT DiagCGaussianMixtureDensity(State *state, FLOAT *obs, Network *net)
{
  size_t  i;
  FLOAT   like = LOG_0;

  assert(state->mOutPdfKind == KID_DiagC);

  if (net && net->outPCache[state->mID].time == net->time) {
    return net->outPCache[state->mID].value;
  }

  for (i = 0; i < state->mNumberOfMixtures; i++) {
    FLOAT glike;
    Mixture *mix = state->mpMixture[i].estimates;

    obs = XFormPass(mix->mpInputXForm, obs,
                    net ? net->time : UNDEF_TIME,
                    net ? net->propagDir : FORWARD);

    assert(obs != NULL);
    glike = DiagCGaussianDensity(mix, obs, net);
    like  = LogAdd(like, glike + state->mpMixture[i].weight);
  }

#ifdef DEBUG_MSGS
  gaus_computaions++;
#endif

  if (net) {
    net->outPCache[state->mID].time = net->time;
    net->outPCache[state->mID].value = like;
  }
  return like;
}

FLOAT FromObservationAtStateId(State *state, FLOAT *obs, Network *net)
{
  obs = XFormPass(net->hmmSet->mpInputXForm, obs,
                  net ? net->time : UNDEF_TIME,
                  net ? net->propagDir : FORWARD);
  assert(obs != NULL);
  return obs[state->PDF_obs_coef];
}

Label *GetLabels(Token *token)
{
  WLR *wlr;
  Label *tmp, *level[3] = {NULL, NULL, NULL};
  int li = 0;

  if (!token || !IS_ACTIVE(*token)) {
    return NULL;
  }

  for (wlr = token->wlr; wlr != NULL; wlr = wlr->next) {
    if ((tmp = (Label *) malloc(sizeof(Label))) == NULL) {
      Error("Insufficient memory");
    }

    tmp->score = wlr->like;
    tmp->stop  = wlr->time;
    tmp->id    = wlr->state_idx;

    if (wlr->node->mType & NT_Model) {
      li = wlr->state_idx >= 0 ? 0 : 1;
      tmp->data = wlr->node->hmm;
      tmp->mpName = wlr->node->hmm->mpMacro->mpName;
      tmp->nextLevel = level[li+1] ? level[li+1] : level[2];
    } else if (wlr->node->pronun->outSymbol) {
      li = 2;
      tmp->data = wlr->node->pronun->word;
      tmp->mpName = wlr->node->pronun->outSymbol;
      tmp->nextLevel = NULL;
    } else {
      free(tmp);
      continue;
    }

    if (level[li]) {

   // if previous label mis its label on lower level, just make it
/*      if (li > 0 && (!level[li-1] || level[li-1]->nextLevel != level[li])) {
        Label *tmp2;
        if ((tmp2 = (Label *) malloc(sizeof(Label))) == NULL) {
          Error("Insufficient memory");
        }
        tmp2->nextLevel = level[li];
        tmp2->mpName  = level[li]->mpName;
        tmp2->score = level[li]->score;
        tmp2->stop  = level[li]->stop;
        tmp2->id    = level[li]->id;

        if (level[li-1]) {
          level[li-1]->score -= tmp2->score;
          level[li-1]->start  = tmp2->stop;
        }

        tmp2->next  = level[li-1];
        level[li-1] = tmp2;
      }*/

      level[li]->score -= tmp->score;
      level[li]->start = tmp->stop;
    }

    tmp->next = level[li];
    level[li] = tmp;
  }
  for (li = 0; li < 3; li++) if (level[li]) level[li]->start = 0;
  return level[0] ? level[0] : level[1] ? level[1] : level[2];
}

void ViterbiInit(Network *net)
{
  net->PassTokenInModel   = &PassTokenMax;
  net->PassTokenInNetwork = &PassTokenMax;
  net->propagDir = FORWARD;
  net->hmmSet->ResetXFormInstances();

  net->time = 0; // Must not be set to -net->hmmSet->mTotalDelay yet
                 // otherwise token cannot enter first model node
                 // with start set to 0
  TokenPropagationInit(net);
  net->time = -net->hmmSet->mTotalDelay;
}

#ifdef DEBUG_MSGS
void PrintNumOfWLR() {
  WLR *wlr = firstWLR;
  int nwlrs = 0, nwlrsnf = 0;
  while (wlr) {
    nwlrs++;
    if (!wlr->freed) {
      nwlrsnf++;
    }
    wlr = wlr->tmpNext;
  }
  printf("%d Released: %d\n", nwlrsnf, nwlrs - nwlrsnf);
}
#endif

void ViterbiStep(Network *net, FLOAT *observation)
{
  net->time++;
  net->hmmSet->UpdateStacks(observation, net->time, net->propagDir);

  if (net->time <= 0) {
    return;
  }

  TokenPropagationInModels(net, observation);
  TokenPropagationInNetwork(net);
#ifdef DEBUG_MSGS
  printf("Frame: %ld Nuberm of WLR Active: ", net->time); PrintNumOfWLR();
#endif

}

FLOAT ViterbiDone(Network *net, Label **labels)
{
  FLOAT totLike = LOG_0;
  if (labels) {
    if (IS_ACTIVE(*net->last->exitToken)) {
      totLike = net->last->exitToken->like;
      *labels =  GetLabels(net->last->exitToken);
    } else {
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


/*void BaumWelchInit(ModelSet *hmms) {
  //ResetAccumsForHMMSet(hmms);
  hmms.ResetAccums();
}*/

struct FWBWRet {
  double totLike;
  FloatInLog avgAccuracy;
};

struct FWBWRet ForwardBackward(Network *net, FLOAT *obsMx, int nFrames)
{
  int i;
  Cache *outPCache;
  struct FWBWRet ret;
  ModelSet *hmms = net->hmmSet;

  outPCache = (Cache *) malloc(nFrames * hmms->mNStates * sizeof(Cache));

  if (outPCache == NULL) {
    Error("Insufficient memory");
  }

  for (i = 0; i < nFrames * hmms->mNStates; i++) {
    outPCache[i].time  = UNDEF_TIME;
    outPCache[i].value = LOG_0;
  }

  free(net->outPCache);
  net->outPCache = NULL;

  net->PassTokenInModel   = &PassTokenSum;
  net->PassTokenInNetwork = &PassTokenSum;
  net->alignment          = NO_ALIGNMENT;


  //Forward Pass
  net->propagDir = FORWARD;
  net->collectAlphaBeta = 1;
  net->time = -hmms->mTotalDelay;

  net->hmmSet->ResetXFormInstances();
  for (i = 0; i < hmms->mTotalDelay; i++) {
    net->time++;
    hmms->UpdateStacks(obsMx + hmms->mInputVectorSize * i, net->time, FORWARD);
  }

  TokenPropagationInit(net);

  for (i = hmms->mTotalDelay; i < nFrames+hmms->mTotalDelay; i++) {
    net->time++;
    net->outPCache = outPCache + hmms->mNStates * (net->time-1);

    net->hmmSet->UpdateStacks(obsMx + hmms->mInputVectorSize * i,
                              net->time, net->propagDir);

    TokenPropagationInModels(net,  obsMx + hmms->mInputVectorSize * i);
    TokenPropagationInNetwork(net);
  }

  if (!IS_ACTIVE(*net->last->exitToken)) { // No token survivered
    TokenPropagationDone(net);
    FreeFWBWRecords(net);
    net->outPCache = outPCache;
    ret.totLike = LOG_0;
    return ret;
  }

  ret.totLike = net->last->exitToken->like; //  totalLikelihood;
  TokenPropagationDone(net);

  //Backward Pass
  net->propagDir = BACKWARD;
  net->time = nFrames+hmms->mTotalDelay;

  for (i = nFrames + hmms->mTotalDelay - 1; i >= nFrames; i--) {
//  We do not need any features, in backward prop. All output probab. are cached.
//  UpdateStacks(hmms, obsMx + hmms->mInputVectorSize * i,
//                 net->time, BACKWARD);
    net->time--;
  }

  net->outPCache = NULL;  // TokenPropagationInit would reset valid 1st frm cache
  TokenPropagationInit(net);

  for (i = nFrames-1; i >= 0; i--) {
    net->outPCache = outPCache + hmms->mNStates * (net->time-1);

//  We do not need any features, in backward prop. All output probab. are cached.
//  UpdateStacks(net->hmmSet, obsMx + hmms->mInputVectorSize * i,     |
//                 net->time, net->propagDir);                   |
//                                                               V
    TokenPropagationInModels(net,  NULL); //obsMx + hmms->mInputVectorSize * i);
    net->time--;
    TokenPropagationInNetwork(net);
  }

  net->outPCache = outPCache;

  if (!IS_ACTIVE(*net->first->exitToken)) { // No token survivered
    TokenPropagationDone(net);
    FreeFWBWRecords(net);
    ret.totLike = LOG_0;
    return ret;
  }

  ret.totLike = HIGHER_OF(ret.totLike, net->first->exitToken->like); //  totalLikelihood;
  // Backward pass P can differ from forward pass P because of the precision
  // problems. Take the higher one to decrease the possibility of getting
  // an occupation probability (when normalizing by P) higher that one.

  FloatInLog fil_ret_totLike = {ret.totLike, 0};
  ret.avgAccuracy  = FIL_Div(net->first->exitToken->accuracy, fil_ret_totLike);
  TokenPropagationDone(net);

  // There may be remaining records in alphaBetaListReverse unused in
  // backward pass. Free them and set alphaBetaListReverse's to NULLs;

  Node *node;
  for (node = net->first; node != NULL; node = node->next) {
    if (!(node->mType & NT_Model)) continue;

    while (node->alphaBetaListReverse) {
      FWBWR *fwbwr = node->alphaBetaListReverse;
      node->alphaBetaListReverse = fwbwr->next;
      free(fwbwr);
    }
  }

  return ret;
}

FLOAT MCEReest(Network *net, FLOAT *obsMx, FLOAT *obsMx2, int nFrames, FLOAT weight, FLOAT sigSlope)
{
  struct FWBWRet fwbw;
  FLOAT TP, P, F;

//  FLOAT updateDir;
  int i, j, k;
  int t;
  ModelSet *hmmsAlig = net->hmmSet;
  ModelSet *hmmsUpdt = net->hmmSetToUpdate;
  Node *node;

  net->accumType = AT_ML;
  net->PassTokenInModel   = &PassTokenSum;
  net->PassTokenInNetwork = &PassTokenSum;
  net->propagDir          = FORWARD;
  net->alignment          = NO_ALIGNMENT;

  net->SearchPaths        = Network::SP_TrueOnly;
  net->hmmSet->ResetXFormInstances();

  net->time = 0; // Must not be set to -net->hmmSet->totalDelay yet
                 // otherwise token cannot enter first model node
                 // with start set to 0
  TokenPropagationInit(net);
  net->time = -net->hmmSet->mTotalDelay;

  for (t = 0; t < nFrames+hmmsAlig->mTotalDelay; t++) {
    ViterbiStep(net, obsMx + hmmsAlig->mInputVectorSize * t);
  }

  TP = net->last->exitToken->like;
  ViterbiDone(net, NULL);

  if (TP <= LOG_MIN) return LOG_0;


  ////////////////// Denominator accumulation //////////////////
  net->SearchPaths = Network::SP_All;

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

  for (node = net->first; node != NULL; node = node->next) {
    if (node->mType & NT_Model && node->alphaBetaList != NULL &&
       node->alphaBetaList->time == 0) {
      node->alphaBetaListReverse = node->alphaBetaList;
      node->alphaBetaList = node->alphaBetaList->next;
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


  // net->mixPCache might be used to cache likelihoods of mixtures of target
  // models. Reallocate the cache to fit mixtures of both models and reset it.
  k = HIGHER_OF(hmmsUpdt->mNMixtures, hmmsAlig->mNMixtures);
  net->mixPCache = (Cache *) realloc(net->mixPCache, k * sizeof(Cache));
  if (net->mixPCache == NULL) Error("Insufficient memory");

  for (i = 0; i < k; i++) net->mixPCache[i].time = UNDEF_TIME;

// Update accumulators
  for (net->time = 0; net->time < nFrames; net->time++) {//for every frame
    FLOAT *obs  =obsMx +hmmsAlig->mInputVectorSize*(net->time+hmmsAlig->mTotalDelay);
    FLOAT *obs2 =obsMx2+hmmsUpdt->mInputVectorSize*(net->time+hmmsUpdt->mTotalDelay);
    hmmsAlig->UpdateStacks(obs, net->time, FORWARD);
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->UpdateStacks(obs2, net->time, FORWARD);
    }

    for (node = net->first; node != NULL; node = node->next) { //for every model
//    for (k=0; k < net->nnodes; k++) {
//      Node *node = &net->nodes[k];
      if (node->mType & NT_Model &&
         node->alphaBetaList != NULL &&
         node->alphaBetaList->time == net->time+1) {

        struct AlphaBeta *st;
        int Nq       = node->hmm->mNStates;
        st = node->alphaBetaList->state;

        for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
          if (st[j].alpha + st[j].beta - P > MIN_LOG_WEGIHT) {
            assert(node->alphaBetaListReverse->time == net->time);

            ReestState(net, node, j-1,
                       (st[j].alpha + st[j].beta - P)  * net->ocpScale,
                        -weight, obs, obs2);
          }
        }
        if (node->alphaBetaListReverse) free(node->alphaBetaListReverse);
        node->alphaBetaListReverse = node->alphaBetaList;
        node->alphaBetaList = node->alphaBetaList->next;
      }
    }
  }

  for (node = net->first; node != NULL; node = node->next) {
    if (node->alphaBetaListReverse != NULL)
      free(node->alphaBetaListReverse);
  }


  ////////////////// Numerator accumulation //////////////////
  net->SearchPaths = Network::SP_TrueOnly;


  ForwardBackward(net, obsMx, nFrames);

  for (node = net->first; node != NULL; node = node->next) {
    if (node->mType & NT_Model && node->alphaBetaList != NULL &&
       node->alphaBetaList->time == 0) {
      node->alphaBetaListReverse = node->alphaBetaList;
      node->alphaBetaList = node->alphaBetaList->next;
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


  // net->mixPCache might be used to cache likelihoods of mixtures of target
  // models. Reallocate the cache to fit mixtures of both models and reset it.
  k = HIGHER_OF(hmmsUpdt->mNMixtures, hmmsAlig->mNMixtures);
  net->mixPCache = (Cache *) realloc(net->mixPCache, k * sizeof(Cache));
  if (net->mixPCache == NULL) Error("Insufficient memory");

  for (i = 0; i < k; i++) net->mixPCache[i].time = UNDEF_TIME;

// Update accumulators
  for (net->time = 0; net->time < nFrames; net->time++) {//for every frame
    FLOAT *obs  =obsMx +hmmsAlig->mInputVectorSize*(net->time+hmmsAlig->mTotalDelay);
    FLOAT *obs2 =obsMx2+hmmsUpdt->mInputVectorSize*(net->time+hmmsUpdt->mTotalDelay);
    hmmsAlig->UpdateStacks(obs, net->time, FORWARD);
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->UpdateStacks(obs2, net->time, FORWARD);
    }

    for (node = net->first; node != NULL; node = node->next) 
    { //for every model
      if (node->mType & NT_Model &&
         node->alphaBetaList != NULL &&
         node->alphaBetaList->time == net->time+1) 
      {

        struct AlphaBeta *st;
        int Nq       = node->hmm->mNStates;
        FLOAT *aq    = node->hmm->        mpTransition->matrix;
        FLOAT *aqacc = node->hmmToUpdate->mpTransition->matrix + SQR(Nq);
//        int qt_1 = (net->nNetStates * net->time) + node->estate_id;
//        int qt = qt_1 + net->nNetStates;

        st = node->alphaBetaList->state;

        if (//!net->mmi_den_pass &&
           st[Nq-1].alpha + st[Nq-1].beta - TP > MIN_LOG_WEGIHT) 
        {
          for (i = 0; i < Nq - 1; i++) 
          {
            LOG_INC(aqacc[i * Nq + Nq-1], aq[i * Nq + Nq-1]  * net->tranScale +
                                         (st[i].alpha                         +
                                          st[Nq-1].beta - TP) * net->ocpScale);
          }
        }

        for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
          if (st[j].alpha + st[j].beta - TP > MIN_LOG_WEGIHT) {
            FLOAT bjtO =net->outPCache[hmmsAlig->mNStates * net->time +
                                       node->hmm->state[j-1]->mID].value;
            // ForwardBackward() set net->outPCache to contain out prob. for all frames

            assert(node->alphaBetaListReverse->time == net->time);

//            if (!net->mmi_den_pass) {
            for (i = 0; i < Nq - 1; i++) {
              LOG_INC(aqacc[i * Nq + j],
                      aq[i * Nq + j]    * net->tranScale +
                      (node->alphaBetaListReverse->state[i].alpha +
                      bjtO              * net->outpScale +
                      st[j].beta - TP)   * net->ocpScale);
            }
//            }

            ReestState(net, node, j-1,
                       (st[j].alpha + st[j].beta - TP)  * net->ocpScale,
                        weight, obs, obs2);

// For True MCE
//            ReestState(net, node, j-1,
//                       (st[j].alpha + st[j].beta - TP + LogAdd(TP,P) - P)  * net->ocpScale,
//                        weight, obs, obs2);

          }
        }
        
        if (node->alphaBetaListReverse) 
          free(node->alphaBetaListReverse);
          
        node->alphaBetaListReverse = node->alphaBetaList;
        node->alphaBetaList = node->alphaBetaList->next;
      }
    }
  }

  for (node = net->first; node != NULL; node = node->next) {
    if (node->alphaBetaListReverse != NULL)
      free(node->alphaBetaListReverse);
  }

  net->accumType = AT_MCE;
  return TP;
}

FLOAT BaumWelchReest(Network *net, FLOAT *obsMx, FLOAT *obsMx2, int nFrames, FLOAT weight)
{
  struct FWBWRet fwbw;
  FLOAT P, updateDir;
  int i, j, k;
  ModelSet *hmmsAlig = net->hmmSet;
  ModelSet *hmmsUpdt = net->hmmSetToUpdate;
  Node *node;

  fwbw = ForwardBackward(net, obsMx, nFrames);
  P = fwbw.totLike;
  if (P < LOG_MIN) return LOG_0;

#ifdef MOTIF
  FLOAT *ocprob = (FLOAT *) malloc(net->nNetStates * (nFrames+1) * sizeof(FLOAT));
  for (i=0; i<net->nNetStates * (nFrames+1); i++) ocprob[i] = 0;
#endif

  for (node = net->first; node != NULL; node = node->next) {
    if (node->mType & NT_Model && node->alphaBetaList != NULL &&
       node->alphaBetaList->time == 0) {
      node->alphaBetaListReverse = node->alphaBetaList;
      node->alphaBetaList = node->alphaBetaList->next;
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


  // net->mixPCache might be used to cache likelihoods of mixtures of target
  // models. Reallocate the cache to fit mixtures of both models and reset it.
  k = HIGHER_OF(hmmsUpdt->mNStates, hmmsAlig->mNStates);
  net->mixPCache = (Cache *) realloc(net->mixPCache, k * sizeof(Cache));
  if (net->mixPCache == NULL) Error("Insufficient memory");

  for (i = 0; i < k; i++) net->mixPCache[i].time = UNDEF_TIME;

// Update accumulators
  for (net->time = 0; net->time < nFrames; net->time++) {//for every frame
    FLOAT *obs  =obsMx +hmmsAlig->mInputVectorSize*(net->time+hmmsAlig->mTotalDelay);
    FLOAT *obs2 =obsMx2+hmmsUpdt->mInputVectorSize*(net->time+hmmsUpdt->mTotalDelay);
    hmmsAlig->UpdateStacks(obs, net->time, FORWARD);
    if (hmmsAlig != hmmsUpdt) {
      hmmsUpdt->UpdateStacks(obs2, net->time, FORWARD);
    }

    for (node = net->first; node != NULL; node = node->next) { //for every model
//    for (k=0; k < net->nnodes; k++) {
//      Node *node = &net->nodes[k];
      if (node->mType & NT_Model &&
         node->alphaBetaList != NULL &&
         node->alphaBetaList->time == net->time+1) {

        struct AlphaBeta *st;
        int Nq       = node->hmm->mNStates;
        FLOAT *aq    = node->hmm->        mpTransition->matrix;
        FLOAT *aqacc = node->hmmToUpdate->mpTransition->matrix + SQR(Nq);
//        int qt_1 = (net->nNetStates * net->time) + node->estate_id;
//        int qt = qt_1 + net->nNetStates;

        st = node->alphaBetaList->state;

        if (//!net->mmi_den_pass &&
           st[Nq-1].alpha + st[Nq-1].beta - P > MIN_LOG_WEGIHT) {
          for (i = 0; i < Nq - 1; i++) {
            LOG_INC(aqacc[i * Nq + Nq-1], aq[i * Nq + Nq-1]  * net->tranScale +
                                         (st[i].alpha                         +
                                          st[Nq-1].beta - P) * net->ocpScale);
          }
        }

        for (j = 1; j < Nq - 1; j++) {                   //for every emitting state
          if (st[j].alpha + st[j].beta - P > MIN_LOG_WEGIHT) {
#ifdef MOTIF
          ocprob[net->nNetStates * (net->time+1) + node->estate_id + j]
//           = node->phoneAccuracy;
            = exp(st[j].alpha+st[j].beta-P) *
              ((1-2*st[j].alphaAccuracy.negative) * exp(st[j].alphaAccuracy.logvalue - st[j].alpha) +
               (1-2*st[j].betaAccuracy.negative)  * exp(st[j].betaAccuracy.logvalue  - st[j].beta)
               - (1-2*fwbw.avgAccuracy.negative) * exp(fwbw.avgAccuracy.logvalue)
              );

#endif
//            int qt_1   = qt - net->nNetStates;
            FLOAT bjtO =net->outPCache[hmmsAlig->mNStates * net->time +
                                       node->hmm->state[j-1]->mID].value;
            // ForwardBackward() set net->outPCache to contain out prob. for all frames

            assert(node->alphaBetaListReverse->time == net->time);

//            if (!net->mmi_den_pass) {
            for (i = 0; i < Nq - 1; i++) {
              LOG_INC(aqacc[i * Nq + j],
                      aq[i * Nq + j]    * net->tranScale +
                      (node->alphaBetaListReverse->state[i].alpha +
                      bjtO              * net->outpScale +
                      st[j].beta - P)   * net->ocpScale);
            }
//            }

            if (net->accumType == AT_MFE || net->accumType == AT_MPE) {
              updateDir = (1-2*st[j].alphaAccuracy.negative) * exp(st[j].alphaAccuracy.logvalue - st[j].alpha) +
                          (1-2*st[j].betaAccuracy.negative)  * exp(st[j].betaAccuracy.logvalue  - st[j].beta)  -
                          (1-2*fwbw.avgAccuracy.negative)    * exp(fwbw.avgAccuracy.logvalue);
            } else {
              updateDir = 1.0;
            }

            ReestState(net, node, j-1,
                       (st[j].alpha + st[j].beta - P)  * net->ocpScale,
                        updateDir*weight, obs, obs2);
          }
        }
        if (node->alphaBetaListReverse) free(node->alphaBetaListReverse);
        node->alphaBetaListReverse = node->alphaBetaList;
        node->alphaBetaList = node->alphaBetaList->next;
      }
    }
  }

  for (node = net->first; node != NULL; node = node->next) {
    if (node->alphaBetaListReverse != NULL)
      free(node->alphaBetaListReverse);
  }

#ifdef MOTIF
  FLOAT max = LOG_0;
  printf("tranScale: %f\noutpScale: %f\n",net->tranScale, net->outpScale);
  for (i = 0; i < net->nNetStates * (nFrames+1); i++) max = HIGHER_OF(max, ocprob[i]);

  imagesc(ocprob, net->nNetStates, (nFrames+1),
          sizeof(FLOAT) == 4 ? "float" : "double",
          NULL,  cm_color, "OccProb");

  for (j=0; j<nFrames+1; j++) {
    for (i=0; i<net->nNetStates; i++) {
      printf("%6.3g ", (ocprob[net->nNetStates * j + i]));
    }
    printf("\n");
  }
//  for (i = 0; i < net->nNetStates * (nFrames+1); i++) {
//    if (ocprob[i] < LOG_MIN) ocprob[i] = max;
//  }
//  imagesc(ocprob, net->nNetStates, (nFrames+1), STR(FLOAT), NULL, cm_color, "Log OccProb");
  free(ocprob);
#endif


  return P;
}


FLOAT ViterbiReest(Network *net, FLOAT *obsMx, FLOAT *obsMx2, int nFrames, FLOAT weight)
{
  int t;
  WLR *wlr;
  Node *prevnode = NULL;
  FLOAT P;
  Cache *outPCache;
  ModelSet *hmmsAlig = net->hmmSet;
  ModelSet *hmmsUpdt = net->hmmSetToUpdate;

  outPCache = (Cache *) malloc(nFrames * hmmsAlig->mNStates * sizeof(Cache));
  if (outPCache == NULL) Error("Insufficient memory");

  for (t = 0; t < nFrames * hmmsAlig->mNStates; t++) {
    outPCache[t].time  = UNDEF_TIME;
    outPCache[t].value = LOG_0;
  }

  free(net->outPCache);
  net->outPCache = NULL;
  net->alignment = STATE_ALIGNMENT;
  ViterbiInit(net);
  nFrames += hmmsAlig->mTotalDelay;

  for (t = 0; t < nFrames; t++) {
    if (t >= hmmsAlig->mTotalDelay) {
      net->outPCache = outPCache + hmmsAlig->mNStates*(t-hmmsAlig->mTotalDelay);
    }
    ViterbiStep(net, obsMx + hmmsAlig->mInputVectorSize * t);
  }

  net->outPCache = outPCache;

  if (!IS_ACTIVE(*net->last->exitToken)) {
    ViterbiDone(net, NULL);
    return LOG_0;
  }

  hmmsAlig->ResetXFormInstances();

  if (hmmsAlig != hmmsUpdt) {
    hmmsUpdt->ResetXFormInstances();
  }

  // invert order of WRLs
  wlr = net->last->exitToken->wlr;
  
  while (wlr->next != NULL) 
  {
    WLR *tmp = wlr->next->next;
    wlr->next->next = net->last->exitToken->wlr;
    net->last->exitToken->wlr = wlr->next;
    wlr->next = tmp;
  }

  for (net->time = -hmmsAlig->mTotalDelay; net->time < 0; net->time++) 
  {
    FLOAT *obs = obsMx+hmmsAlig->mInputVectorSize*(net->time+hmmsAlig->mTotalDelay);
    hmmsAlig->UpdateStacks(obs, net->time, FORWARD);
  }

  if (hmmsAlig != hmmsUpdt) 
  {
    FLOAT *obs2= obsMx2+hmmsUpdt->mInputVectorSize*(net->time+hmmsUpdt->mTotalDelay);
    for (net->time = -hmmsUpdt->mTotalDelay; net->time < 0; net->time++) 
    {
      hmmsUpdt->UpdateStacks(obs2, net->time, FORWARD);
    }
  }

// Update accumulators
  for (wlr = net->last->exitToken->wlr; wlr != NULL; wlr = wlr->next) 
  {
    Node *node   = wlr->node;
    int Nq       = node->hmmToUpdate->mNStates;
    FLOAT *aqacc = node->hmmToUpdate->mpTransition->matrix + SQR(Nq);
    int currstate = wlr->state_idx+1;
    int nextstate = (wlr->next && node == wlr->next->node)
                    ? wlr->next->state_idx+1 : Nq-1;
    int duration  = wlr->time - net->time;



    if (prevnode != node) {
//      if (!net->mmi_den_pass)
      LOG_INC(aqacc[currstate], 0 /*ln(1)*/);
      prevnode = node;
    }

//    if (!net->mmi_den_pass) { // So far we dont do any MMI estimation of trasitions
    LOG_INC(aqacc[currstate * Nq + currstate], log(duration-1));
    LOG_INC(aqacc[currstate * Nq + nextstate], 0 /*ln(1)*/);
//    }

    for (; net->time < wlr->time; net->time++) {
      //for every frame of segment
      FLOAT *obs =obsMx +hmmsAlig->mInputVectorSize*(net->time+hmmsAlig->mTotalDelay);
      FLOAT *obs2=obsMx2+hmmsUpdt->mInputVectorSize*(net->time+hmmsUpdt->mTotalDelay);

      hmmsAlig->UpdateStacks(obs, net->time, FORWARD);
      if (hmmsAlig != hmmsUpdt) {
        hmmsUpdt->UpdateStacks(obs2, net->time, FORWARD);
      }
      ReestState(net, node, currstate-1, 0.0, 1.0*weight, obs, obs2);
    }
  }

  P = net->last->exitToken->wlr->like;
  ViterbiDone(net, NULL);
  return P;
}

void UpdateXFormStatCache(XFormStatCache *xfsc,
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
    if (xfsc->norm > 0) 
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
    xfsc->norm++;
  } 
  
  else if (topXForm->mXFormType == XT_COMPOSITE) 
  {
    CompositeXForm *cxf = (CompositeXForm *) topXForm;
    for (i=0; i<cxf->layer[0].mNBlocks; i++) 
    {
      UpdateXFormStatCache(xfsc, cxf->layer[0].block[i], input);
      input += cxf->layer[0].block[i]->mInSize;
    }
  }
}

void UpdateXFormInstanceStatCaches(XFormInstance *xformInstance,
                                   FLOAT *observation, int time)
{
  int i, j;
  FLOAT *obs;

  if (xformInstance == NULL || xformInstance->statCacheTime == time) return;

  xformInstance->statCacheTime = time;

  if (xformInstance->mNumberOfXFormStatCaches == 0) return;

  if (xformInstance->mpInput) {
    UpdateXFormInstanceStatCaches(xformInstance->mpInput, observation, time);
  }

  for (i = 0; i < xformInstance->mNumberOfXFormStatCaches; i++) {
    XFormStatCache *xfsc = &xformInstance->mpXFormStatCache[i];

    if (xfsc->mpUpperLevelStats != NULL &&
       xfsc->mpUpperLevelStats->mpStats == xfsc->mpStats) { //just link to upper level?
       xfsc->norm = xfsc->mpUpperLevelStats->norm;
      continue;
    }

    obs = XFormPass(xformInstance->mpInput, observation, time, FORWARD);
    xfsc->norm = 0;
    UpdateXFormStatCache(xfsc, xformInstance->mpXForm, obs);
    if (xfsc->mpUpperLevelStats != NULL) {
      int size = xfsc->mpXForm->mInSize;
      for (j = 0; j < size + size*(size+1)/2; j++) {
        xfsc->mpStats[j] += xfsc->mpUpperLevelStats->mpStats[j];
      }
      xfsc->norm += xfsc->mpUpperLevelStats->norm;
    }
  }
}

void ReestState(Network *net, Node *node,
                int state_idx, FLOAT logPriorProb, FLOAT updateDir,
                FLOAT *obs, FLOAT *obs2) {
  int i, j, k, m;
  State *state  = node->hmm->        state[state_idx];
  State *state2 = node->hmmToUpdate->state[state_idx];
  FLOAT bjtO    = LOG_0;
  int nmixtures;

  if (!net->hmmSetToUpdate->mGaussLvl2ModelReest
  && net->hmmSet != net->hmmSetToUpdate) {
    // Occupation probabilities of mixtures are computed using target model
    state = state2; obs = obs2;
  } else if (state->mNumberOfMixtures <= state2->mNumberOfMixtures) {
    bjtO = net->outPCache[net->hmmSet->mNStates * net->time + state->mID].value;
  }

  nmixtures = LOWER_OF(state->mNumberOfMixtures, state2->mNumberOfMixtures);

  if (bjtO < LOG_MIN) {
    // State likelihood was not available in cache because
    // - occupation probabilities of mixtures are computed using target model
    // - not all mixtures of alignment model are used for computation of
    //   occupation probabilities (state2->num_mix < state->num_mix)
    for (m = 0; m < nmixtures; m++) {
      Mixture *mix = state->mpMixture[m].estimates;
      FLOAT cjm    = state->mpMixture[m].weight;
      FLOAT *xobs  = XFormPass(mix->mpInputXForm, obs, net->time, FORWARD);
      FLOAT bjmtO  = DiagCGaussianDensity(mix, xobs, net);
      bjtO  = LogAdd(bjtO, cjm + bjmtO);
    }
  }

  for (m = 0; m < nmixtures; m++) {                  //for every emitting state
    Mixture *mix = state->mpMixture[m].estimates;
    FLOAT cjm    = state->mpMixture[m].weight;
    FLOAT *xobs  = XFormPass(mix->mpInputXForm, obs, net->time, FORWARD);
    FLOAT bjmtO  = DiagCGaussianDensity(mix, xobs, net);
    FLOAT Lqjmt  = logPriorProb - bjtO + cjm + bjmtO;

    if (Lqjmt > MIN_LOG_WEGIHT) {
      Mixture *mix       = state2->mpMixture[m].estimates;
      int vec_size       = mix->mpMean->mVectorSize;
      FLOAT *mnvec       = mix->mpMean->mVector;
      FLOAT *mnacc       = mix->mpMean->mVector     +     vec_size;
      FLOAT *vvacc       = mix->mpVariance->mVector +     vec_size;
      FLOAT *vmacc       = mix->mpVariance->mVector + 2 * vec_size;
      XFormInstance *ixf = mix->mpInputXForm;
      FLOAT *xobs        = XFormPass(ixf, obs2, net->time, FORWARD);

/*      if (net->mmi_den_pass) {
        mnacc += vec_size + 1;
        vvacc += 2 * vec_size + 1;
        vmacc += 2 * vec_size + 1;
      }*/

      Lqjmt = exp(Lqjmt) * updateDir;

      for (i = 0; i < vec_size; i++) {                 // Update
        mnacc[i] += Lqjmt * xobs[i];                  // mean
        if (net->hmmSetToUpdate->mUpdateMask & UM_OLDMEANVAR) {
          vvacc[i] += Lqjmt * SQR(xobs[i]-mnvec[i]);  // var
        } else {
          vvacc[i] += Lqjmt * SQR(xobs[i]);           // scatter
          vmacc[i] += Lqjmt * xobs[i];                //mean for var.
        }
      }

      mnacc[vec_size] += Lqjmt; //norms for mean
      vmacc[vec_size] += Lqjmt; //norms for variance

//      if (net->mmi_den_pass) {
//        state2->mpMixture[m].weight_accum_den += Lqjmt;
//      } else {
      state2->mpMixture[m].weight_accum     += Lqjmt; //  Update weight accum
//      }

      if (Lqjmt < 0) {
        state2->mpMixture[m].weight_accum_den -= Lqjmt;
      }

      if (ixf == NULL || ixf->mNumberOfXFormStatCaches == 0) continue;

      UpdateXFormInstanceStatCaches(ixf, obs2, net->time);

      for (i = 0; i < ixf->mNumberOfXFormStatCaches; i++) {
        XFormStatCache *xfsc = &ixf->mpXFormStatCache[i];
        Variance *var  = state2->mpMixture[m].estimates->mpVariance;
        Mean     *mean = state2->mpMixture[m].estimates->mpMean;

        for (j = 0; j < var->mNumberOfXFormStatAccums; j++) {
          XFormStatAccum *xfsa = &var->mpXFormStatAccum[j];
          if (xfsa->mpXForm == xfsc->mpXForm) {
            int size = xfsc->mpXForm->mInSize;
            for (k = 0; k < size+size*(size+1)/2; k++) {
              xfsa->mpStats[k] += xfsc->mpStats[k] * Lqjmt;
            }
            xfsa->norm += xfsc->norm * Lqjmt;
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
            xfsa->norm += xfsc->norm * Lqjmt;
            break;
          }
        }
      }
    }
  }
}


/*FLOAT *StateOccupationProbability(Network *net, FLOAT *obsMx, ModelSet *hmms,
                                  int nFrames, FLOAT **outProbOrMahDist, int getMahalDist)
{
  int i, j, k;
  int nNetModels;
  int nEmitingStates;
  FLOAT totalLike;
  FLOAT *occupProb;

  FLOAT *beta = (FLOAT *) malloc(net->nNetStates * (nFrames+1) * sizeof(FLOAT));
  FLOAT *alfa = (FLOAT *) malloc(net->nNetStates * (nFrames+1) * sizeof(FLOAT));
  Cache *outPCache = (Cache *) malloc(nFrames * hmms->mNStates * sizeof(Cache));

  if (outPCache == NULL || beta == NULL || alfa == NULL) {
    Error("Insufficient memory");
  }

  for (i = 0; i < net->nNetStates * (nFrames+1); i++) {
    beta[i] = alfa[i] = LOG_0;
  }

  nNetModels = 0;
  for (i=0; i < net->nnodes; i++) {
    if (net->nodes[i].type & NT_Model) nNetModels++;
  }

  nEmitingStates = net->nNetStates - 2 * nNetModels;

  for (i = 0; i < nFrames * hmms->mNStates; i++) {
    outPCache[i].time = UNDEF_TIME;
    outPCache[i].value = LOG_0;
  }

  free(net->outPCache);


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
  net->alignment          = NO_ALIGNMENT;

  //Backward Pass
  net->propagDir = BACKWARD;

  net->time = nFrames;
  TokenPropagationInit(net,
                       beta + net->nNetStates * net->time,
                       NULL);

  for (i = nFrames-1; i >= 0; i--) {

    net->outPCache = outPCache + hmms->mNStates * i;
    TokenPropagationInModels(net,  obsMx + hmms->mInputVectorSize * i,
                             beta + net->nNetStates * net->time,
                             NULL);

    if (outProbOrMahDist != NULL && getMahalDist !=4) {
      int state_counter = 0;

      for (k=0; k < net->nnodes; k++) {
        Node *node = &net->nodes[k];
        if (node->mType & NT_Model) {
          for (j = 0; j < node->hmm->mNStates - 2; j++, state_counter++) {
            FLOAT tmpf = net->OutputProbability(node->hmm->state[j],
                                       obsMx + hmms->mInputVectorSize * i, net);
            switch (getMahalDist) {
              case 0:
                break;
              case 1:
               tmpf = log((tmpf / -0.5) - node->hmm->state[j]->mpMixture[0].estimates->mGConst) * -0.5;
               break;
              case 2:
               tmpf = log((tmpf / -0.5) - node->hmm->state[j]->mpMixture[0].estimates->mGConst) * -1;
               break;
              case 3:
               tmpf += node->hmm->state[j]->mpMixture[0].estimates->mGConst * 0.5;
//               tmpf /= hmms->mInputVectorSize;
               break;
            }

            (*outProbOrMahDist)[i * nEmitingStates + state_counter] = tmpf;
          }
        }
      }
    }

    net->time--;
    TokenPropagationInNetwork(net,
                              beta + net->nNetStates * net->time, NULL);
  }

  if (!IS_ACTIVE(*net->first->exitToken)) {
    TokenPropagationDone(net);
    net->outPCache = outPCache;
    free(alfa);
    for (i=0; i < net->nNetStates * nFrames; i++) {
      beta[i] = LOG_0;
    }
    return beta; // No token survivered
  }

  totalLike = net->first->exitToken->like;
//  totalLike = 0;
  TokenPropagationDone(net);

  //Forward Pass
  net->propagDir = FORWARD;
  net->time = 0;
  TokenPropagationInit(net,
                       beta + net->nNetStates * net->time,
                       alfa + net->nNetStates * net->time);

  for (i = 0; i < nFrames; i++) {
    net->outPCache = outPCache + hmms->mNStates * i;
    net->time++;
    TokenPropagationInModels(net,  obsMx + hmms->mInputVectorSize * i,
                             beta + net->nNetStates * net->time,
                             alfa + net->nNetStates * net->time);

    TokenPropagationInNetwork(net,
                              beta + net->nNetStates * net->time,
                              alfa + net->nNetStates * net->time);

    if (outProbOrMahDist != NULL && getMahalDist == 4) {
      for (j=0; j < nEmitingStates; j++) {
        (*outProbOrMahDist)[i * nEmitingStates + j] = totalLike;
      }
    }
  }

  net->outPCache = outPCache;

  TokenPropagationDone(net);

  for (i = 1; i <= nFrames; i++) {
    int state_counter = 0;

    for (k=0; k < net->nnodes; k++) {
      Node *node = &net->nodes[k];
      if (node->mType & NT_Model) {
        for (j = 0; j < node->hmm->mNStates - 2; j++, state_counter++) {
          int idx = (net->nNetStates * i) + node->estate_id + j + 1;
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

WLR *TimePruning(Network *net, int frame_delay)
{
  size_t    i;
  Node *    node;
  Token *   token = net->bestToken;
  WLR *     twlr;
  WLR *     rwlr=NULL;

  if (frame_delay > net->time-1 || !token) 
    return NULL;

  if (token->twlr != NULL &&
     token->twlr->time == net->time-1 - frame_delay) 
  {
    rwlr = token->twlr;
  }

  for (node = net->first; node != NULL; node = node->next) 
  {
    if (!(node->mType & NT_Model)) 
      continue;

    for (i = 0; i < node->hmm->mNStates-1; i++) 
    {
      if (IS_ACTIVE(node->tokens[i])) 
      {
        Token *token = &node->tokens[i];

        if (token->twlr != NULL &&
           token->twlr->time == net->time-1 - frame_delay) 
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
              for (twlr = token->wlr; twlr->next != token->twlr; twlr = twlr->next)
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
