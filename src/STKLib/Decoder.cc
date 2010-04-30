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
#include "mymath.h"

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
#ifndef NDEBUG
  int test_for_cycle = 0;
  int HasCycleCounter = 100000;
#endif

  
#ifdef DEBUG_MSGS
  WordLinkRecord * firstWLR;
#endif
  
  
int nbest_lattices = -1;
std::map<std::string,FLOAT> state_posteriors;


#if 0
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
#endif 

      
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
      FLOAT     glike = Decoder<DecoderNetwork>::DiagCGaussianDensity(mix, l_obs, NULL) + pState->mpMixture[i].mWeight;

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
          NetworkType::Node *node = &net->mpNodes[k];
          if (node->mType & NT_MODEL) {
            for (j = 0; j < node->mpHmm->mNStates - 2; j++, state_counter++) {
              FLOAT tmpf = net->OutputProbability(node->mpHmm->mpState[j],
                                        obsMx + hmms->mInputVectorSize * i, net);
              switch (getMahalDist) {
                case 0:
                  break;
                case 1:
                tmpf = my_log((tmpf / -0.5) - node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->GConst()) * -0.5;
                break;
                case 2:
                tmpf = my_log((tmpf / -0.5) - node->mpHmm->mpState[j]->mpMixture[0].mpEstimates->GConst()) * -1;
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
        NetworkType::Node *node = &net->mpNodes[k];
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
  /*void 
  Decoder::
  AddLinkToLattice(NetworkType::Node *from, NetworkType::Node *to, FLOAT lmLike)
  {
    assert(from->mpAnr && from->mpAnr && from->mpAnr->mpExitToken->mpWlr);
    long long from_time = from->mpAnr->mpExitToken->mpWlr->mpNext 
                          ? from->mpAnr->mpExitToken->mpWlr->mpNext->mTime : 0;
    
    std::ostringstream ss;
    ss << ios::hex
       << from
       << from_time;
               
    NetworkType::Node * lattice_from = find_or_create_node(&mLatticeNodeHash, ss.str().c_str(), &mpLatticeLastNode);
    lattice_from->mType  = from->mType;
    lattice_from->mpName = from->mpName; // !!! Relys on union, stays for all mpName, mpHmm and mpPronun

    ss.str("");
    ss << ios::hex 
       << to
       << from->mpAnr->mpExitToken->mpWlr->mTime;

    NetworkType::Node * lattice_to = find_or_create_node(&mLatticeNodeHash, ss.str().c_str(), &mpLatticeLastNode);
    lattice_to->mType  = to->mType;
    lattice_to->mpName = to->mpName; // !!! Relys on union, stays for all mpName, mpHmm and mpPronun
    
    int nl = ++lattice_from->rNLinks();

    lattice_from->rpLinks() = (NetworkType::LinkType *) realloc(lattice_from->rpLinks(), nl * sizeof(NetworkType::LinkType));
    if (lattice_from->rpLinks() == NULL) Error("Insufficient memory");
          
    lattice_from->rpLinks()[nl-1].mpNode = lattice_to;
    lattice_from->rpLinks()[nl-1].mLike = lmLike;
    lattice_from->rpLinks()[nl-1].mpNode->rNBackLinks()++;
    
  }*/
  // AddLinkToLattice(NetworkType::Node *from, NetworkType::Node *to, FLOAT lmLike)
  //***************************************************************************
  

  
    

  //***************************************************************************
  //***************************************************************************
  // void
  // Decoder:: 
  // SortNodes()
  // {
  //   int       i;
  //   int       j;
  //   NetworkType::Node *    chain;
  //   NetworkType::Node *    last;
  //   NetworkType::Node *    node;
  // 
  //   // Sort nodes for forward (Viterbi) propagation
  // 
  //   for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
  //   {
  //     node->mAux = node->rNBackLinks();
  //   }
  // 
  //   for (i = 0; i < rNetwork().pFirst()->rNLinks(); i++) 
  //   {
  //     rNetwork().pFirst()->rpLinks()[i].pNode()->mAux--;
  //   }
  // 
  //   last = rNetwork().pFirst();
  //   chain = rNetwork().pFirst()->mpNext;
  // 
  //   while (chain) 
  //   {
  //     bool    short_curcuit = true;
  //     NetworkType::Node ** curPtr = &chain;
  //     i = 0;
  // 
  //     while (*curPtr) 
  //     {
  //       if ((((*curPtr)->mType & NT_MODEL) && !((*curPtr)->mType & NT_TEE))
  //         || (*curPtr)->mAux == 0) 
  //       {
  //         for (j = 0; j < (*curPtr)->rNLinks(); j++) 
  //         {
  //           (*curPtr)->rpLinks()[j].pNode()->mAux--;
  //         }
  // 
  //         last = (last->mpNext = *curPtr);
  //         last->mAux = i++;
  //         *curPtr = (*curPtr)->mpNext;
  //         short_curcuit = false;
  //       } 
  //       else 
  //       {
  //         curPtr = &(*curPtr)->mpNext;
  //       }
  //     }
  // 
  //     if (short_curcuit) 
  //     {
  // //      fprintf(stderr, "Nodes in loop: ");
  // //      for (curPtr = &chain; *curPtr; curPtr = &(*curPtr)->next)
  // //        fprintf(stderr, "%d %d", *curPtr - mpNodes, (*curPtr)->mType);
  // //      fprintf(stderr, "\n");
  //       Error("Loop of non-emiting nodes found in network");
  //     }
  //   }
  // 
  //   last->mpNext = NULL;
  // 
  //   /// !!! What is this sorting links good for ???
  //   for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
  //   {
  //     if (node->rNLinks() > 1)
  //       qsort(node->rpLinks(), node->rNLinks(), sizeof(NetworkType::LinkType), cmplnk);
  //   }
  // 
  // // Sort nodes for backward propagation
  // 
  //   for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext)
  //     node->mAux = node->rNLinks();
  // 
  //   for (i = 0; i < rNetwork().pLast()->rNBackLinks(); i++)
  //     rNetwork().pLast()->rpBackLinks()[i].pNode()->mAux--;
  // 
  //   last = rNetwork().pLast();
  //   chain = rNetwork().pLast()->mpBackNext;
  //   i = 0;
  // 
  //   while (chain) 
  //   {
  //     bool short_curcuit = true;
  //     NetworkType::Node **curPtr = &chain;
  // 
  //     while (*curPtr) 
  //     {
  //       if ((((*curPtr)->mType & NT_MODEL) && !((*curPtr)->mType & NT_TEE))
  //         || (*curPtr)->mAux == 0) 
  //       {
  //         for (j = 0; j < (*curPtr)->rNBackLinks(); j++) 
  //         {
  //           (*curPtr)->rpBackLinks()[j].pNode()->mAux--;
  //         }
  // 
  //         last = (last->mpBackNext = *curPtr);
  //         last->mAux = i++;
  //         *curPtr = (*curPtr)->mpBackNext;
  //         short_curcuit = false;
  //       } 
  //       else 
  //       {
  //         curPtr = &(*curPtr)->mpBackNext;
  //       }
  //     }
  // 
  //     assert(!short_curcuit); // Shouldn't happen, since it didnot happen before
  //   }
  // 
  //   last->mpBackNext = NULL;
  // 
  //   /// !!! What is this sorting links good for ???
  //   for (node = rNetwork().pFirst(); node != NULL; node = node->mpNext) 
  //   {
  //     if (node->rNBackLinks() > 1)
  //       qsort(node->rpBackLinks(), node->rNBackLinks(), sizeof(NetworkType::LinkType), cmplnk);
  //   }
  // } // Decoder::SortNodes();
  

}; // namespace STK



