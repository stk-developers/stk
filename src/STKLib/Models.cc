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

#include "Models.h"
#include "Decoder.h"
#include "stkstream.h"
#include "mymath.h"
#include "common.h"
#include "mymath.h"

// MATLAB Engine
#ifdef MATLAB_ENGINE
#  include "engine.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include <functional>
#include <map>



namespace STK
{

//  const char*     gpHListFilter;
//  bool            gHmmsIgnoreMacroRedefinition = true;
  FLOAT           gWeightAccumDen;
  
  
  FunctionTable gFuncTable[] = 
  {
    {sigmoid_vec, KID_Sigmoid}, {log_vec,     KID_Log},
    {exp_vec,     KID_Exp},     {sqrt_vec,    KID_Sqrt},
    {softmax_vec, KID_SoftMax}
  };
  
  size_t gFuncTableSize = sizeof(gFuncTable)/sizeof(*gFuncTable);
  enum StatType {MEAN_STATS, COV_STATS, CMLLR_STATS};
  
  
  //***************************************************************************
  //***************************************************************************
  void
  CatMVectors(BasicVector<FLOAT>& rV, Matrix<FLOAT>* pM, size_t nM)
  {
    size_t j(0);
    
    for (size_t n = 0; n < nM; n++)
    {
      for (size_t i = 0; i < pM[n].Cols(); i++)
      {
        rV[j++] = pM[n][0][i];
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void
  CatClusterWeightPartialVectors(BasicVector<FLOAT>& rV, BiasXform** pM, size_t nM)
  {
    size_t j(0);
    
    for (size_t n = 0; n < nM; n++)
    {
      for (size_t i = 0; i < pM[n]->mInSize; i++)
      {
        rV[j++] = pM[n]->mVector[0][i];
      }
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Mixture&
  Mixture::
  AddToClusterWeightVectorsAccums(Matrix<FLOAT>* pGw, BasicVector<FLOAT>* pKw)
  {
    FLOAT          occ_counts;                                                  // occupation counts for current mixture
    Matrix<FLOAT>  aux_mat1(mpMean->mClusterMatrixT);                           // auxiliary matrix
    
    // old vector
    aux_mat1.DiagScale(mpVariance->mVector);
    
    // we go through each accumulator set
    for (size_t sub=0; sub < mpMean->mNClusterWeightVectors; sub++)
    {
      occ_counts = mpMean->mpOccProbAccums[sub];                
      
      pGw[sub].AddCMMtMul(occ_counts, aux_mat1, mpMean->mClusterMatrixT);
      pKw[sub].AddCMVMul(1, aux_mat1, mpMean->mCwvAccum[sub]);
    } // for (size_t sub=0; sub < mpMean->mCwvAccum.Rows())
        
    return *this;
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  ComputeClusterWeightVectorAccums(
    int               macro_type, 
    HMMSetNodeName    nodeName, 
    MacroData*        pData, 
    void*             pUserData)
  {
    ClusterWeightAccumUserData* cwa = reinterpret_cast<ClusterWeightAccumUserData*>(pUserData);
    Mixture*             mix = reinterpret_cast<Mixture*>(pData);
    
    if (mix->mpMean->mNClusterWeightVectors > 0)
    {
      mix->AddToClusterWeightVectorsAccums(cwa->mpGw, cwa->mpKw);
    }
  }  
  
  //***************************************************************************
  //***************************************************************************
  Mixture&
  Mixture::
  ResetClusterWeightVectorsAccums()
  {
    for (size_t i = 0; i < mpMean->mNClusterWeightVectors; i++)
    {
      mpMean->mpOccProbAccums[i] = 0.0;
    }
    
    mpMean->mCwvAccum.Clear();
    return *this;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ReplaceItem(
    int               macro_type, 
    HMMSetNodeName    nodeName, 
    MacroData *       pData, 
    void *            pUserData)
  {
    ReplaceItemUserData*   ud = (ReplaceItemUserData*) pUserData;
    size_t                  i;
    size_t                  j;
  
    if (macro_type == 'h') 
    {
      Hmm *hmm = (Hmm *) pData;
      
      if (ud->mType == 's') 
      {
        for (i = 0; i < hmm->mNStates-2; i++) 
        {
          if (hmm->mpState[i] == ud->mpOldData) 
          {
            hmm->mpState[i] = (State*) ud->mpNewData;
          }
        }
      } 
      else if (hmm->mpTransition == ud->mpOldData) 
      {
        hmm->mpTransition = (Transition*) ud->mpNewData;
      }
    } 
    
    else if (macro_type == 's') 
    {
      State *state = (State *) pData;
      if (state->mOutPdfKind != KID_PDFObsVec) 
      {
//!!! Test for macro_type == 'm'
        for (i = 0; i < state->mNMixtures; i++) 
        {
          if (state->mpMixture[i].mpEstimates == ud->mpOldData) 
          {
            state->mpMixture[i].mpEstimates = (Mixture*) ud->mpNewData;
          }
        }
      }
    } 
    
    else if (macro_type == 'm') 
    {
      Mixture* mixture = static_cast<Mixture*>(pData);
      if (mixture->mpMean      == ud->mpOldData) mixture->mpMean       = (Mean*)          ud->mpNewData;
      if (mixture->mpVariance  == ud->mpOldData) mixture->mpVariance   = (Variance*)      ud->mpNewData;
      if (mixture->mpInputXform== ud->mpOldData) mixture->mpInputXform = (XformInstance*) ud->mpNewData;
      
      // For cluster adaptive training, BiasXform holds the weights vector. If
      // new bias is defined, we want to check all means and potentially update
      // the weights refference and recalculate the values
      if ('x' == ud->mType                                               
      && (XT_BIAS == static_cast<Xform*>(ud->mpNewData)->mXformType))
      {
        Mean* mean(mixture->mpMean);
        
        // find if any of the mpWeights points to the old data
        i=0;
        while (i < mean->mNClusterWeightVectors 
        &&    (mean->mpClusterWeightVectors[i] != static_cast<BiasXform*>(ud->mpOldData)))
        {i++;}
        
        // if we found something... 
        if (i < mean->mNClusterWeightVectors)
        {
          // the cluster parameters accumulators need to be updated as speaker 
          // has changed
          mixture->UpdateClusterParametersAccums();
          
          BiasXform* new_weights = static_cast<BiasXform*>(ud->mpNewData);
          
          mean->mpClusterWeightVectors[i] = new_weights;
          mean->RecalculateCAT();
          
          if (mean->mCwvAccum.Rows() > 0)
          {
            mean->ResetClusterWeightVectorsAccums(i);
          }
        }        
      }
    } 
    
    else if (macro_type == 'x') 
    {
      CompositeXform *cxf = (CompositeXform *) pData;
      if (cxf->mXformType == XT_COMPOSITE) 
      {
        for (i = 0; i < cxf->mNLayers; i++) 
        {
          for (j = 0; j < cxf->mpLayer[i].mNBlocks; j++) 
          {
            if (cxf->mpLayer[i].mpBlock[j] == ud->mpOldData) 
            {
              cxf->mpLayer[i].mpBlock[j] = (Xform*) ud->mpNewData;
            }
          }
        }
      }
    } 
    
    else if (macro_type == 'j') 
    {
      XformInstance *xformInstance = (XformInstance *) pData;
      if (xformInstance->mpInput == ud->mpOldData) xformInstance->mpInput = (XformInstance*) ud->mpNewData;
      if (xformInstance->mpXform == ud->mpOldData) xformInstance->mpXform = (Xform*)         ud->mpNewData;
    }
  }
    
  //***************************************************************************
  //***************************************************************************
  bool 
  IsXformIn1stLayer(Xform *xform, Xform *topXform)
  {
    size_t      i;
    if (topXform == NULL)  return false;
    if (topXform == xform) return true;
  
    if (topXform->mXformType == XT_COMPOSITE) 
    {
      CompositeXform *cxf = static_cast<CompositeXform *>(topXform);
      
      for (i=0; i < cxf->mpLayer[0].mNBlocks; i++) 
      {
        if (IsXformIn1stLayer(xform, cxf->mpLayer[0].mpBlock[i])) 
          return true;
      }
    }
    
    return false;
  }
  
  //***************************************************************************
  //***************************************************************************
  bool 
  Is1Layer1BlockLinearXform(Xform * pXform)
  {
    CompositeXform *cxf = static_cast<CompositeXform *>(pXform);
    
    if (cxf == NULL)                                        return false;
    if (cxf->mXformType == XT_LINEAR)                       return true;
    if (cxf->mXformType != XT_COMPOSITE)                    return false;
    if (cxf->mNLayers > 1 || cxf->mpLayer[0].mNBlocks > 1)  return false;
    
    return Is1Layer1BlockLinearXform(cxf->mpLayer[0].mpBlock[0]);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  Macro*
  FindMacro(MyHSearchData *macro_hash, const char *name) 
  {
    ENTRY   e={0}; // {0} is just to make compiler happy
    ENTRY*  ep; 
    e.key = (char *) name;
    my_hsearch_r(e, FIND, &ep, macro_hash);
    return (Macro *) (ep ? ep->data : NULL);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ReleaseMacroHash(MyHSearchData *macro_hash) 
  {
    unsigned int i;
    for (i = 0; i < macro_hash->mNEntries; i++) 
    {
      Macro *macro = (Macro *) macro_hash->mpEntry[i]->data;
      free(macro->mpName);
      free(macro->mpFileName);
      delete macro;
      macro_hash->mpEntry[i]->data = NULL;
    }
    my_hdestroy_r(macro_hash, 0);
  }
  
  
  //***************************************************************************
  //***************************************************************************
  void 
  ReleaseItem(int macro_type,HMMSetNodeName nodeName,MacroData * pData, void * pUserData)
  {
    delete pData;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  qsmacrocmp(const void *a, const void *b) 
  {
    return strcmp(((Macro *)a)->mpName, ((Macro *)b)->mpName);
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  ResetAccum(int macro_type, HMMSetNodeName nodeName,
             MacroData* pData, void* pUserData) 
  {
    size_t    i;
    size_t    size    = 0;
    FLOAT*    vector  = NULL;
    ModelSet  *p_mode_set = reinterpret_cast<ModelSet *>(pUserData);
  
    if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      if (macro_type == mt_mean) {
        size   = reinterpret_cast<Mean *>(pData)->VectorSize();
        vector = reinterpret_cast<Mean *>(pData)->mpAccums;
        size   = (size + 1) * p_mode_set->mAllocAccums;     // mAllocAccums == 2 for MMI update
      } else if (macro_type == mt_variance) {
        size   = reinterpret_cast<Variance *>(pData)->VectorSize();
        vector = reinterpret_cast<Variance *>(pData)->mpAccums;
        size   = (size * 2 + 1) * p_mode_set->mAllocAccums; // mAllocAccums == 2 for MMI update
      }
  
      //for (i = 0; i < size; i++) vector[i] = 0;
      memset(vector, 0, size*sizeof(FLOAT));
  
    } 
    else if (macro_type == mt_state) 
    {
      State *state = (State *) pData;
      if (state->mOutPdfKind == KID_DiagC) 
      {
        for (i = 0; i < state->mNMixtures; i++) 
        {
          state->mpMixture[i].mWeightAccum     = 0;
          state->mpMixture[i].mWeightAccumDen  = 0;
        }
      }
    } 
    else if (macro_type == mt_transition) 
    {
      size   = SQR(reinterpret_cast<Transition *>(pData)->mNStates);
      vector =     reinterpret_cast<Transition *>(pData)->mpMatrixO + size;
  
      for (i = 0; i < size; i++) vector[i] = LOG_0;
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  GlobalStats(
    int                 macro_type, 
    HMMSetNodeName      nn, 
    MacroData*          pData, 
    void*               pUserData)
  {
    size_t        i;
    
    GlobalStatsUserData *ud = (GlobalStatsUserData *) pUserData;
  
    if (macro_type == mt_state) 
    {
      State * state = (State *) pData;
      
      if (state->mOutPdfKind == KID_DiagC) 
      {
        for (i = 0; i < state->mNMixtures; i++) 
        {
          state->mpMixture[i].mWeightAccum += 1;
        }
      }
    } 
    else 
    { // macro_type == mt_mixture
      Mixture * mixture   = (Mixture *) pData;
      size_t    vec_size  = mixture->mpMean->VectorSize();
      FLOAT *   obs       = XformPass(mixture->mpInputXform,ud->observation,ud->mTime,FORWARD);
      
      for (i = 0; i < vec_size; i++) 
      {
        //mixture->mpMean->mpVectorO[vec_size + i] += obs[i];
        mixture->mpMean->mpAccums[i] += obs[i];
      }
      
      //mixture->mpMean->mpVectorO[2 * vec_size] += 1;
      mixture->mpMean->mpAccums[vec_size] += 1;
  
      for (i = 0; i < vec_size; i++) 
      {
        //mixture->mpVariance->mpVectorO[vec_size  +i] += SQR(obs[i]);
        //mixture->mpVariance->mpVectorO[2*vec_size+i] += obs[i];
        mixture->mpVariance->mpAccums[i] += SQR(obs[i]);
        mixture->mpVariance->mpAccums[vec_size+i] += obs[i];
      }
      
      //mixture->mpVariance->mpVectorO[3 * vec_size] += 1;
      mixture->mpVariance->mpAccums[2 * vec_size] += 1;
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  FLOAT*
  XformPass(
    XformInstance*        pXformInst, 
    FLOAT*                pInputVector, 
    int                   time, 
    PropagDirectionType   dir)
  {
    if (pXformInst == NULL) {
      return pInputVector;
    }
  
    if (time != UNDEF_TIME && pXformInst->mTime == time) {
      return pXformInst->pOutputData();
    }
  
    pXformInst->mTime = time;
  
    // recursively pass previous transformations
    if (pXformInst->mpInput) {
      pInputVector = XformPass(pXformInst->mpInput, pInputVector, time, dir);
    }
  
    // evaluate this transformation
    pXformInst->mpXform->Evaluate(
        pInputVector,
        pXformInst->pOutputData(),
        pXformInst->mpMemory,
        dir);
  
    return pXformInst->pOutputData();
  }

    
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  AllocXformStatAccums(
    XformStatAccum **     xformStatAccum,
    size_t         *      nxformStatAccums,
    XformInstance  *      xformInstance,
    enum StatType         stat_type) 
  {
    size_t    i;
    size_t    j;
    
    if (xformInstance == NULL) 
      return;
  
    for (i = 0; i < xformInstance->mNumberOfXformStatCaches; i++) 
    {
      XformStatCache *xfsc = &xformInstance->mpXformStatCache[i];
      XformStatAccum *xfsa = *xformStatAccum;
  
      for (j = 0; j < *nxformStatAccums; j++, xfsa++) 
      {
        if (xfsa->mpXform == xfsc->mpXform) 
          break;
      }
  
      if (j == *nxformStatAccums) 
      {
        size_t size = xfsc->mpXform->mInSize; //mean : mean+covariance
        size = (stat_type == MEAN_STATS) ? size : 
               (stat_type == COV_STATS)  ? size+size*(size+1)/2
                                         : 0; // CMLLR_STATS
  
        *xformStatAccum =
          (XformStatAccum *) realloc(*xformStatAccum,
                                    sizeof(XformStatAccum) * ++*nxformStatAccums);
  
        if (*xformStatAccum == NULL)
          Error("Insufficient memory");
  
        xfsa = *xformStatAccum + *nxformStatAccums - 1;
        
        if(size == 0) {
          xfsa->mpStats = NULL;
        } 
        else if ((xfsa->mpStats = (FLOAT *) malloc(sizeof(FLOAT) * size)) == NULL) {
          Error("Insufficient memory");
        }
  
        xfsa->mpXform    = xfsc->mpXform;
        xfsa->mNorm     = 0.0;
        
        for (j = 0; j < size; j++) 
          xfsa->mpStats[j] = 0.0;
      }
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  AllocateXformStatCachesAndAccums(
    int                 macro_type, 
    HMMSetNodeName      nodeName,
    MacroData       *   pData, 
    void            *   hmm_set)
  {
  
    ModelSet *p_hmm_set = reinterpret_cast<ModelSet *>(hmm_set);
    if (macro_type == mt_XformInstance) 
    {
      //Allocate Xform stat caches for XformInstance
  
      XformInstance *   xfi =(XformInstance *) pData;
      size_t            i;
      size_t            j;
  
      for (i=0; i < p_hmm_set->mNumberOfXformsToUpdate; i++) 
      {
        Xform *xform = p_hmm_set->mpXformToUpdate[i].mpXform;
        int instanceContainXfrom = IsXformIn1stLayer(xform, xfi->mpXform);
  
        //Does instance one level up contain cache for this xform
        XformStatCache *upperLevelStats = NULL;
        if (xfi->mpInput != NULL) {
          for (j=0; j < xfi->mpInput->mNumberOfXformStatCaches; j++) {
            if (xfi->mpInput->mpXformStatCache[j].mpXform == xform) {
              upperLevelStats = &xfi->mpInput->mpXformStatCache[j];
              break;
            }
          }
        }
  
        if (instanceContainXfrom || upperLevelStats != NULL) 
        {
          XformStatCache *xfsc;
  
          xfi->mpXformStatCache = static_cast<XformStatCache*>(
              realloc(xfi->mpXformStatCache, sizeof(XformStatCache) * 
                ++xfi->mNumberOfXformStatCaches));
  
          if (xfi->mpXformStatCache == NULL) {
            Error("Insufficient memory");
          }
  
          xfsc = &xfi->mpXformStatCache[xfi->mNumberOfXformStatCaches-1];
  
          if (instanceContainXfrom) {
            int size = xform->mInSize;
            size = size+size*(size+1)/2;
  
            if ((xfsc->mpStats = (FLOAT *) malloc(sizeof(FLOAT) * size))==NULL) {
              Error("Insufficient memory");
            }
          } else {
            xfsc->mpStats = upperLevelStats->mpStats;
          }
  
          xfsc->mNorm = 0;
          xfsc->mpXform = xform;
          xfsc->mpUpperLevelStats = upperLevelStats;
        }
      }
    } else if (macro_type == mt_mixture) {
      //Allocate Xform stat accumulators for mean and covariance
  
      Mixture *mix = (Mixture *) pData;
      AllocXformStatAccums(&mix->mpMean->mpXformStatAccum,
                          &mix->mpMean->mNumberOfXformStatAccums,
                          mix->mpInputXform, p_hmm_set->mCmllrStats ? CMLLR_STATS : MEAN_STATS);
  
      AllocXformStatAccums(&mix->mpVariance->mpXformStatAccum,
                          &mix->mpVariance->mNumberOfXformStatAccums,
                          mix->mpInputXform, p_hmm_set->mCmllrStats ? CMLLR_STATS : COV_STATS);
  
      if (mix->mpInputXform == NULL || mix->mpInputXform->mNumberOfXformStatCaches == 0)
        return;
        
      if (mix->mpInputXform->mNumberOfXformStatCaches != 1 ||
        !Is1Layer1BlockLinearXform(mix->mpInputXform->mpXform) ||
        mix->mpInputXform->mpXformStatCache[0].mpUpperLevelStats != NULL) 
      {
        mix->mpVariance->mUpdatableFromStatAccums = false;
        mix->mpMean    ->mUpdatableFromStatAccums = false;
        p_hmm_set->mAllMixuresUpdatableFromStatAccums = false;  
      } 
      
      else if (mix->mpMean->mNumberOfXformStatAccums != 1) 
      {
        assert(mix->mpMean->mNumberOfXformStatAccums > 1);
        mix->mpMean->mUpdatableFromStatAccums = false;
        p_hmm_set->mAllMixuresUpdatableFromStatAccums = false;
      } 
      
      else if (mix->mpVariance->mNumberOfXformStatAccums != 1) 
      {
        assert(mix->mpVariance->mNumberOfXformStatAccums > 1);
        mix->mpVariance->mUpdatableFromStatAccums = false;
        p_hmm_set->mAllMixuresUpdatableFromStatAccums = false;
      }
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  NormalizeStatsForXform(int macro_type, HMMSetNodeName nodeName,
                         MacroData * pData, void * pUserData) 
  {
    XformStatAccum *xfsa = NULL;
    int i, j, k, nxfsa = 0, size;
    FLOAT *mean, *cov, inorm;
  
    if (macro_type == mt_mean) 
    {
      xfsa  = ((Mean *)pData)->mpXformStatAccum;
      nxfsa = ((Mean *)pData)->mNumberOfXformStatAccums;
    } 
    
    else if (macro_type == mt_variance) 
    {
      xfsa  = ((Variance *)pData)->mpXformStatAccum;
      nxfsa = ((Variance *)pData)->mNumberOfXformStatAccums;
    }
  
    for (i = 0; i < nxfsa; i++) 
    {
      size = xfsa[i].mpXform->mInSize;
      mean = xfsa[i].mpStats;
      cov  = xfsa[i].mpStats + size;
      inorm = 1.0 / xfsa[i].mNorm;
  
      for (j = 0; j < size; j++) 
        mean[j] *= inorm; //normalize means
  
      if (macro_type == mt_variance) 
      {
        for (k=0; k < size; k++) 
        {
          for (j=0; j <= k; j++) 
          {                 //normalize covariances
            cov[k*(k+1)/2+j] = cov[k*(k+1)/2+j] * inorm - mean[k] * mean[j];
          }
        }
      }
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  WriteStatsForXform(int macro_type, HMMSetNodeName nodeName,
                     MacroData * pData, void * pUserData) 
  {
    XformStatsFileNames*     file = NULL;
    XformStatAccum*               xfsa = NULL;
    size_t                        i;
    size_t                        j;
    size_t                        k;
    size_t                        nxfsa = 0;
    int                           cc = 0;
    size_t                        size;
    FLOAT*                        mean;
    FLOAT*                        cov;
    WriteStatsForXformUserData*   ud = (WriteStatsForXformUserData *) pUserData;
  
    if (macro_type == mt_mean) 
    {
      file  = &ud->mMeanFile;
      xfsa  = ((Mean *)pData)->mpXformStatAccum;
      nxfsa = ((Mean *)pData)->mNumberOfXformStatAccums;
    } 
    else if (macro_type == mt_variance) 
    {
      file  = &ud->mCovFile;
      xfsa  = ((Variance *)pData)->mpXformStatAccum;
      nxfsa = ((Variance *)pData)->mNumberOfXformStatAccums;
    }
  
    for (i = 0; i < nxfsa && xfsa[i].mpXform != (Xform *) ud->mpXform; i++)
      ;
    
    if (i == nxfsa) 
      return;
  
    if (fprintf(file->mpOccupP, "%s ", nodeName) < 0 ||
        WriteNumber(file->mpOccupP, xfsa[i].mNorm) < 0 ||
        fputs("\n", file->mpOccupP) == EOF)
    {
      Error("Cannot write to file: %s", file->mpOccupN);
    }
  
    size = xfsa[i].mpXform->mInSize;
    mean = xfsa[i].mpStats;
    cov  = xfsa[i].mpStats + size;
  
    if (macro_type == mt_mean) 
    {
      if (ud->mBinary) 
      {
        if (!isBigEndian()) 
          for (i = 0; i < size; i++) swapFLOAT(mean[i]);
          
        cc |= (fwrite(mean, sizeof(FLOAT), size, file->mpStatsP) != size);
        
        if (!isBigEndian()) 
          for (i = 0; i < size; i++) swapFLOAT(mean[i]);
      } 
      else 
      {
        for (j=0;j<size;j++) 
        {
          cc |= WriteNumber(file->mpStatsP, mean[j]) < 0;
          cc |= fputc(' ', file->mpStatsP) == EOF;
        }
        
        cc |= fputs("\n", file->mpStatsP) == EOF;
      }
    } 
    else 
    {
      if (ud->mBinary) 
      {
        size = size*(size+1)/2;
        if (!isBigEndian()) for (i = 0; i < size; i++) swapFLOAT(cov[i]);
        cc |= fwrite(cov, sizeof(FLOAT), size, file->mpStatsP) != size;
        if (!isBigEndian()) for (i = 0; i < size; i++) swapFLOAT(cov[i]);
      } 
      else
      {
        for (k=0; k < size; k++) 
        {
          for (j=0;j<=k;j++) {
            cc |= WriteNumber(file->mpStatsP, cov[k*(k+1)/2+j]) < 0;
            cc |= fputc(' ', file->mpStatsP) == EOF;
          }
  
          for (;j<size; j++) {
            cc |= WriteNumber(file->mpStatsP, cov[j*(j+1)/2+k]) < 0;
            cc |= fputc(' ', file->mpStatsP) == EOF;
          }
  
          cc |= fputs("\n", file->mpStatsP) == EOF;
        }
        cc |= fputs("\n", file->mpStatsP) == EOF;
      }
    }
  
    if (cc) {
      Error("Cannot write to file %s", file->mpStatsN);
    }
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  ReadStatsForXform(int macro_type, HMMSetNodeName nodeName,
                    MacroData * pData, void *pUserData) 
  {
    char                        buff[128];
    XformStatsFileNames *  file = NULL;
    XformStatAccum *            xfsa = NULL;
    size_t                      i;
    size_t                      j;
    size_t                      k;
    size_t                      nxfsa = 0;
    int                         cc = 0;
    size_t                      size;
    FLOAT *                     mean;
    FLOAT *                     cov;
    FLOAT                       f;
    WriteStatsForXformUserData * ud = (WriteStatsForXformUserData *) pUserData;
  
    if (macro_type == mt_mean) 
    {
      file  = &ud->mMeanFile;
      xfsa  = ((Mean *)pData)->mpXformStatAccum;
      nxfsa = ((Mean *)pData)->mNumberOfXformStatAccums;
    } 
    else if (macro_type == mt_variance) 
    {
      file  = &ud->mCovFile;
      xfsa  = ((Variance *)pData)->mpXformStatAccum;
      nxfsa = ((Variance *)pData)->mNumberOfXformStatAccums;
    }
  
    for (i = 0; i < nxfsa && xfsa[i].mpXform != (Xform *) ud->mpXform; i++)
      ;
      
    if (i == nxfsa) 
      return;
  
    j = fscanf(file->mpOccupP, "%128s", buff);
    j += ReadNumber(file->mpOccupP, &xfsa[i].mNorm);
    fscanf(file->mpOccupP, "\n");
    
    if (j < 1) 
      Error("Unexpected end of file: %s", file->mpOccupN);
    
    else if (strcmp(buff, nodeName)) 
      Error("'%s' expected but '%s' found in file: %s",nodeName,buff,file->mpOccupN);
    
    else if (j < 2) 
      Error("Decimal number expected after '%s'in file: %s", buff, file->mpOccupN);
    
    size = xfsa[i].mpXform->mInSize;
    mean = xfsa[i].mpStats;
    cov  = xfsa[i].mpStats + size;
  
    if (macro_type == mt_mean) 
    {
      if (ud->mBinary) 
      {
        j = fread(mean, sizeof(FLOAT), size, file->mpStatsP);
        cc |= j != size;
        if (!isBigEndian()) for (i = 0; i < size; i++) swapFLOAT(mean[i]);
      } 
      else 
      {
        for (j=0;j<size;j++)
        {
          cc |= !ReadNumber(file->mpStatsP, &mean[j]);
          fscanf(file->mpStatsP, " ");
        }
      }
    } 
    else 
    {
      if (ud->mBinary) 
      {
        size = size*(size+1)/2;
        cc |= (fread(cov, sizeof(FLOAT), size, file->mpStatsP) != size);
        
        if (!isBigEndian()) 
          for (i = 0; i < size; i++) swapFLOAT(cov[i]);
      } 
      else
      {
        for (k=0; k < size; k++) 
        {
          for (j=0;j<k;j++) 
          {
            cc |= !ReadNumber(file->mpStatsP, &f);
            fscanf(file->mpStatsP, " ");
            if (f != cov[k*(k+1)/2+j]) 
            {
              Error("Covariance matrix '%s' in file '%s' must be symetric",
                    nodeName, file->mpStatsN);
            }
          }
  
          for (;j<size; j++) {
            cc |= !ReadNumber(file->mpStatsP, &cov[j*(j+1)/2+k]);
            fscanf(file->mpStatsP, " ");
          }
        }
      }
    }
  
    if (ferror(file->mpStatsP)) 
    {
      Error("Cannot read file '%s'", file->mpStatsN);
    } 
    else if (cc) 
    {
      Error("Invalid file with Xform statistics '%s'", file->mpStatsN);
    }
  }
  

  //*****************************************************************************
  //*****************************************************************************
  void 
  NormalizeAccum(int macro_type, HMMSetNodeName nodeName,
                 MacroData * pData, void *pUserData) 
  {
    size_t      i;
    size_t      j;
    size_t      size   = 0;
    FLOAT *     vector = NULL;
  
    if (macro_type == mt_mean || macro_type == mt_variance) 
    {
      XformStatAccum *  xfsa = NULL;
      size_t            nxfsa = 0;
  
      if (macro_type == mt_mean) 
      {
        xfsa   = ((Mean *)pData)->mpXformStatAccum;
        nxfsa  = ((Mean *)pData)->mNumberOfXformStatAccums;
        size   = ((Mean *)pData)->VectorSize();
        //vector = ((Mean *)pData)->mpVectorO+size;
        vector = ((Mean *)pData)->mpAccums;
        size   = size + 1;
      } 
      else if (macro_type == mt_variance) 
      {
        xfsa   = ((Variance *)pData)->mpXformStatAccum;
        nxfsa  = ((Variance *)pData)->mNumberOfXformStatAccums;
        size   = ((Variance *)pData)->VectorSize();
        //vector = ((Variance *)pData)->mpVectorO+size;
        vector = ((Variance *)pData)->mpAccums;
        size   = size * 2 + 1;
      }
  
      for (i=0; i < size; i++) 
        vector[i] /= vector[size-1];
  
      for (i = 0; i < nxfsa; i++) 
      {
        size = xfsa[i].mpXform->mInSize;
        size = (macro_type == mt_mean) ? size : size+size*(size+1)/2;
  
        for (j=0; j < size; j++) 
          xfsa[i].mpStats[j] /= xfsa[i].mNorm;
        
        xfsa[i].mNorm = 1.0;
      }
    } 
    else if (macro_type == mt_state) 
    {
      State *state = (State *) pData;
      
      if (state->mOutPdfKind == KID_DiagC) 
      {
        FLOAT accum_sum = 0.0;
  
        for (i = 0; i < state->mNMixtures; i++)
          accum_sum += state->mpMixture[i].mWeightAccum;
  
        if (accum_sum > 0.0) 
        {
          for (i = 0; i < state->mNMixtures; i++)
            state->mpMixture[i].mWeightAccum /= accum_sum;
        }
      }
    } 
    else if (macro_type == mt_transition) 
    {
      size_t nstates = ((Transition *) pData)->mNStates;
      vector = ((Transition *) pData)->mpMatrixO + SQR(nstates);
  
      for (i=0; i < nstates; i++) 
      {
        FLOAT nrm = LOG_0;
        
        for (j=0; j < nstates; j++) 
        {
          LOG_INC(nrm, vector[i * nstates + j]);
        }
        
        if (nrm < LOG_MIN) nrm = 0.0;
        
        for (j=0; j < nstates; j++) 
        {
          vector[i * nstates + j] -= nrm;
        }
      }
    }
  }
  

  //**************************************************************************
  //**************************************************************************
  //   Hmm class section
  //**************************************************************************
  //**************************************************************************
  Hmm::
  Hmm(size_t nStates) : 
    mpTransition(NULL)
  {
    // we allocate pointers for states. The -2 is for the non-emmiting states
    mpState   = new State*[nStates - 2];
    mNStates  = nStates;
  }
  
  //**************************************************************************
  //**************************************************************************
  Hmm::
  ~Hmm()
  {
    delete [] mpState;
  }
  
  //**************************************************************************  
  //**************************************************************************
  void
  Hmm::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    size_t i;
  
    for (i = 0; i < mNStates - 2; i++) 
    {
      if (!this->mpState[i]->mpMacro) 
        mpState[i]->UpdateFromAccums(pModelSet, this);
    }
  
    if (!mpTransition->mpMacro) 
      mpTransition->UpdateFromAccums(pModelSet);
  } // UpdateFromAccums(const ModelSet * pModelSet)

  //*****************************************************************************
  //*****************************************************************************
  //virtual 
  void 
  Hmm::
  Scan(int mask, HMMSetNodeName nodeName,
       ScanAction action, void *pUserData)
  {
    size_t    i;
    size_t    n = 0;
    char *    chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    
    if (mask & MTM_HMM && mask & MTM_PRESCAN) 
      action(mt_hmm, nodeName, this, pUserData);
  
    if (mask & (MTM_ALL & ~(MTM_HMM | MTM_TRANSITION))) 
    {
      for (i=0; i < mNStates-2; i++) 
      {
        if (!mpState[i]->mpMacro) 
        {
          if (n > 0 ) snprintf(chptr, n, ".state[%d]", (int) i+2);
          mpState[i]->Scan(mask, nodeName, action, pUserData);
        }
      }
    }
  
    if (mask & MTM_TRANSITION && !mpTransition->mpMacro) 
    {
      if (n > 0) strncpy(chptr, ".transP", n);
      action(mt_transition, nodeName, mpTransition, pUserData);
    }
  
    if (mask & MTM_HMM && !(mask & MTM_PRESCAN)) 
    {
      if (n > 0) chptr = '\0';
      action(mt_hmm, nodeName, this, pUserData);
    }
  }
  
  
  //**************************************************************************
  //**************************************************************************
  //   Mixture section
  //**************************************************************************
  //**************************************************************************
  void 
  Mixture::
  ComputeGConst()
  {
    FLOAT cov_det = 0;
    size_t i;
  
    for (i = 0; i < mpVariance->VectorSize(); i++) 
    {
      // old vector
      //cov_det -= log(mpVariance->mpVectorO[i]);
      cov_det -= my_log(mpVariance->mVector[i]);
    }
    mGConst = cov_det + M_LOG_2PI * mpVariance->VectorSize();
  }  
  
  //***************************************************************************
  //***************************************************************************
  Variance *
  Mixture::
  FloorVariance(const ModelSet * pModelSet)
  {
    size_t   i;                                   // element index
    FLOAT    g_floor = pModelSet->mMinVariance;   // global varfloor
    Variance *f_floor = NULL;                     // flooring vector
    
    
    f_floor = (mpInputXform != NULL) ? mpInputXform->mpVarFloor 
                                     : pModelSet->mpVarFloor;
    
    // none of the above needs to be set, so 
    if (f_floor != NULL)
    {
      assert(f_floor->VectorSize() == mpVariance->VectorSize());
      
      // go through each element and update
      for (i = 0; i < mpVariance->VectorSize(); i++)
      {
        // old vector
        //mpVariance->mpVectorO[i] = LOWER_OF(mpVariance->mpVectorO[i],
        //                                    f_floor->mpVectorO[i]);
        mpVariance->mVector[i] = LOWER_OF(mpVariance->mVector[i],
                                          f_floor->mVector[i]);
      }
    }
    
    // we still may have global varfloor constant set
    if (g_floor > 0.0)
    {
      // go through each element and update
      g_floor =  1.0 / g_floor;
      for (i = 0; i < mpVariance->VectorSize(); i++)
      {
        // old vector
        //mpVariance->mpVectorO[i] = LOWER_OF(mpVariance->mpVectorO[i], g_floor);
        mpVariance->mVector[i] = LOWER_OF(mpVariance->mVector[i], g_floor);
      }
    }
    
    return mpVariance;
  } // FloorVariance(...)
  
  void
  Mixture::
  ISmoothing(FLOAT numPlusDenNorm,
             const ModelSet * pModelSet)
  {
    if (!pModelSet->mUpdateMask & UM_MAP)
      return;
      
    if(pModelSet->mISmoothingMaxOccup > 0.0
    && numPlusDenNorm > pModelSet->mISmoothingMaxOccup)
    {
      return;
    }

    int i;      
    int vec_size    = mpVariance->VectorSize();
    FLOAT *var_vec = mpVariance->mpAccums;
    FLOAT *mean_vec  = mpVariance->mpAccums + 1 * vec_size;
    FLOAT *norm     = mpVariance->mpAccums + 2 * vec_size;

    FLOAT *var_pri  = mpVariance->mpPrior->mVector.pData();
    FLOAT *mean_pri = mpMean    ->mpPrior->mVector.pData();

    FLOAT tau = pModelSet->mMapTau;

    if(pModelSet->JSmoothing && numPlusDenNorm > *norm)
    {
      tau = tau * fabs(*norm) / numPlusDenNorm;
    }
    for (i = 0; i < vec_size; i++) 
    {
      mean_vec[i] += tau * mean_pri[i];
      var_vec[i]  += tau * (1.0/var_pri[i] + SQR(mean_pri[i]));
    }
    *norm += tau;
  }
  
  //**************************************************************************
  //**************************************************************************
  void 
  Mixture::
  UpdateFromAccums(const ModelSet * pModelSet) 
  {
    double  Djm;
    int     i;

    if (UT_TwoAccumSetEBW == pModelSet->mUpdateType)
    {
      int vec_size    = mpVariance->VectorSize();
      // old vector
      //FLOAT *mean_vec = mpMean->mpVectorO;
      //FLOAT *var_vec  = mpVariance->mpVectorO;
      FLOAT* mean_vec = mpMean->mVector.pData();
      FLOAT* var_vec  = mpVariance->mVector.pData();
  
      //FLOAT *vac_num  = var_vec + 1 * vec_size;
      //FLOAT *mac_num  = var_vec + 2 * vec_size;
      //FLOAT *nrm_num  = var_vec + 3 * vec_size;
      FLOAT *vac_num  = mpVariance->mpAccums;
      FLOAT *mac_num  = mpVariance->mpAccums + 1 * vec_size;
      FLOAT *nrm_num  = mpVariance->mpAccums + 2 * vec_size;
  
      FLOAT *vac_den  = vac_num + 2 * vec_size + 1;
      FLOAT *mac_den  = mac_num + 2 * vec_size + 1;
      FLOAT *nrm_den  = nrm_num + 2 * vec_size + 1;
  
//      FLOAT *var_pri  = mpVariance->mpPrior->mVector.pData();
//      FLOAT *mean_pri = mpMean    ->mpPrior->mVector.pData();

      if (!pModelSet->mISmoothAfterD) 
      {
        // Obsolete way of I-smoothing making use of numearator statistics
        for (i = 0; i < vec_size; i++) 
        {
          mac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
          vac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
        }
        *nrm_num   += pModelSet->MMI_tauI;

        // New way of I-smoothing or general MAP update making use of prior model set
        ISmoothing(*nrm_num + *nrm_den, pModelSet);
/*        if (pModelSet->mUpdateMask & UM_MAP)
        {
          for (i = 0; i < vec_size; i++) {
            mac_num[i] += pModelSet->mMapTau * mean_pri[i];
            vac_num[i] += pModelSet->mMapTau * (1.0/var_pri[i] + SQR(mean_pri[i]));
          }
          *nrm_num += pModelSet->mMapTau;
        }*/
      }
  
      Djm = 0.0;
      
      if (pModelSet->mUpdateMask & UM_VARIANCE) 
      {
        // Find minimum Djm leading to positive update of variances
        for (i = 0; i < vec_size; i++) {
          double macn_macd = mac_num[i]-mac_den[i];
          double vacn_vacd = vac_num[i]-vac_den[i];
          double nrmn_nrmd = *nrm_num - *nrm_den;
          double a  = 1/var_vec[i];
          double b  = vacn_vacd + nrmn_nrmd * (1/var_vec[i] + SQR(mean_vec[i])) -
                      2 * macn_macd * mean_vec[i];
          double c  = nrmn_nrmd * vacn_vacd - SQR(macn_macd);
          double Dd = (- b + sqrt(SQR(b) - 4 * a * c)) / (2 * a);
    
          Djm = HIGHER_OF(Djm, Dd);
        }
      }
  
      Djm = HIGHER_OF(pModelSet->MMI_h * Djm, pModelSet->MMI_E * *nrm_den);
  
      if (pModelSet->mISmoothAfterD) 
      {
        // Obsolete way of I-smoothing making use of numearator statistics
        for (i = 0; i < vec_size; i++) 
        {
          mac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
          vac_num[i] *= (*nrm_num + pModelSet->MMI_tauI) / *nrm_num;
        }
        *nrm_num   += pModelSet->MMI_tauI;

        // New way of I-smoothing or general MAP update making use of prior model set
        ISmoothing(*nrm_num + *nrm_den, pModelSet);
        /*if (pModelSet->mUpdateMask & UM_MAP)
        {
          for (i = 0; i < vec_size; i++) {
            mac_num[i] += pModelSet->mMapTau * mean_pri[i];
            vac_num[i] += pModelSet->mMapTau * (1.0/var_pri[i] + SQR(mean_pri[i]));
          }
          *nrm_num += pModelSet->mMapTau;
        }*/
      }
      
      for (i = 0; i < vec_size; i++) 
      {
        double macn_macd = mac_num[i]-mac_den[i];
        double vacn_vacd = vac_num[i]-vac_den[i];
        double nrmn_nrmd = *nrm_num - *nrm_den;
  
        double new_mean = (macn_macd + Djm * mean_vec[i]) / (nrmn_nrmd + Djm);
        
        if (pModelSet->mUpdateMask & UM_VARIANCE) 
        {
          var_vec[i]   = 1/((vacn_vacd + Djm * (1/var_vec[i] + SQR(mean_vec[i]))) /
                            (nrmn_nrmd + Djm) - SQR(new_mean));
        }
        
        if (pModelSet->mUpdateMask & UM_MEAN) 
        {
          mean_vec[i]  = new_mean;
        }
      }
    } 
    
    // //////////
    // MPE update
    else if ((UT_EBW == pModelSet->mUpdateType) && 0 == mAccumK.Rows()) 
    { 
    

    
      int    vec_size = mpVariance->VectorSize();
      FLOAT* mean_vec = mpMean->mVector.pData();
      FLOAT* var_vec  = mpVariance->mVector.pData();
  
      //FLOAT *vac_mpe  = var_vec + 1 * vec_size;
      //FLOAT *mac_mpe  = var_vec + 2 * vec_size;
      //FLOAT *nrm_mpe  = var_vec + 3 * vec_size;
      FLOAT *vac_mpe  = mpVariance->mpAccums;
      FLOAT *mac_mpe  = mpVariance->mpAccums + 1 * vec_size;
      FLOAT *nrm_mpe  = mpVariance->mpAccums + 2 * vec_size;
  
//      FLOAT *var_pri  = mpVariance->mpPrior->mVector.pData();
//      FLOAT *mean_pri = mpMean    ->mpPrior->mVector.pData();

      if (*nrm_mpe + 2 * gWeightAccumDen <= pModelSet->mMinOccupation) 
      {
//!!!Should not be better to move this after I-smoothing???
        if (mpMacro)
          Warning("Low occupation of '%s', mixture is not updated", mpMacro->mpName);
        else 
          Warning("Low occupation of mixture, mixture is not updated");
        return;
      }
  
      if (!pModelSet->mISmoothAfterD) 
      {
        ISmoothing(*nrm_mpe + 2 * gWeightAccumDen, pModelSet);
        // I-smoothing or general MAP update
/*        if (pModelSet->mUpdateMask & UM_MAP)
        {
          FLOAT tau;
          if(!pModelSet->JSmoothing
          || (*nrm_mpe + 2 * gWeightAccumDen) <= *nrm_mpe) {
            tau = pModelSet->mMapTau;
          } else {
            tau = pModelSet->mMapTau * *nrm_mpe / (*nrm_mpe + 2 * gWeightAccumDen);
          }
          for (i = 0; i < vec_size; i++) 
          {
            mac_mpe[i] += tau * mean_pri[i];
            vac_mpe[i] += tau * (1.0/var_pri[i] + SQR(mean_pri[i]));
          }
          *nrm_mpe += tau;
        }*/
      }
      
      Djm = 0.0;
      
      if (pModelSet->mUpdateMask & UM_VARIANCE) 
      {
        // Find minimum Djm leading to positive update of variances
        for (i = 0; i < vec_size; i++) 
        {
          double macn_macd = mac_mpe[i];
          double vacn_vacd = vac_mpe[i];
          double nrmn_nrmd = *nrm_mpe;
          double a  = 1/var_vec[i];
          double b  = vacn_vacd + nrmn_nrmd * (1/var_vec[i] + SQR(mean_vec[i])) -
                    2 * macn_macd * mean_vec[i];
          double c  = nrmn_nrmd * vacn_vacd - SQR(macn_macd);
          double Dd = (- b + sqrt(SQR(b) - 4 * a * c)) / (2 * a);
  
          Djm = HIGHER_OF(Djm, Dd);
        }
      }
  
//printf("Djm=%f, gWeightAccumDen=%f\n", Djm, gWeightAccumDen);
      // !!! gWeightAccumDen is passed using quite ugly hack that work
      // !!! only if mixtures are not shared by more states - MUST BE REWRITEN
      Djm = HIGHER_OF(pModelSet->MMI_h * Djm, pModelSet->MMI_E * gWeightAccumDen);
      //(*nrm_mpe + 2 * gWeightAccumDen)); SERest_altE
  
      if (pModelSet->mISmoothAfterD) 
      {
        ISmoothing(*nrm_mpe + 2 * gWeightAccumDen, pModelSet);
        // I-smoothing
/*
        if (pModelSet->mUpdateMask & UM_MAP)
        {
          for (i = 0; i < vec_size; i++) {
            mac_mpe[i] += pModelSet->mMapTau * mean_pri[i];
            vac_mpe[i] += pModelSet->mMapTau * (1.0/var_pri[i] + SQR(mean_pri[i]));
          }
          *nrm_mpe += pModelSet->mMapTau;
        }*/
      }
  
      for (i = 0; i < vec_size; i++) 
      {
        double macn_macd = mac_mpe[i];
        double vacn_vacd = vac_mpe[i];
        double nrmn_nrmd = *nrm_mpe;
        double new_mean  = (macn_macd + Djm * mean_vec[i]) / (nrmn_nrmd + Djm);
        if (pModelSet->mUpdateMask & UM_VARIANCE) 
        {
          var_vec[i]     = 1 / ((vacn_vacd + Djm * (1/var_vec[i] + SQR(mean_vec[i]))) /
                                (nrmn_nrmd + Djm) - SQR(new_mean));
        }

        if (pModelSet->mUpdateMask & UM_MEAN)  
        {
          mean_vec[i]    = new_mean;
        }
      }
    }
    
    // Cluster parameters update
    else if (0 < mAccumK.Rows())
    {

      // update accumulators from partial simplified accumulators
      UpdateClusterParametersAccums();

      // if discriminative training:
      if (UT_EBW == pModelSet->mUpdateType)
      {
        FLOAT D       = pModelSet->MMI_E * gWeightAccumDen;
        FLOAT gamma_n = mpVariance->mpAccums[mpVariance->VectorSize()*2] + gWeightAccumDen;
        
        // G_D = G / gamma_n
        Matrix<FLOAT>&      G_D(mAccumGd);
        G_D.DivC(gamma_n);
        
        // K_D = G_D * M^T
        Matrix<FLOAT>      K_D;
        K_D.AddMMMul(G_D, mpMean->mClusterMatrixT);
                
        // L_D = 1/invVar_old + M_old * G_D * M_old' =
        //       1/invVar_old + MG_D * M_old'
        BasicVector<FLOAT> L_D(mpVariance->mVector);
        Matrix<FLOAT>      MG_D(mpMean->mClusterMatrixT.Cols(),
                                mpMean->mClusterMatrixT.Rows());
        
        // :TODO:
        // Test these procedures, they are new...
        MG_D.AddCMtMMul(1.0, mpMean->mClusterMatrixT, G_D);
        L_D.AddDiagCMMMul(1.0, MG_D, mpMean->mClusterMatrixT);
                
        //G += D * G_D
        //K += D * K_D
        //L += D * L_D        
        mAccumG.AddCMMul(D, G_D);
        mAccumK.AddCMMul(D, K_D);
        mAccumL.AddCVMul(D, L_D);
        
        //:KLUDGE^2:
        // this is a total mess
        mpVariance->mpAccums[mpVariance->VectorSize()*2] += D;
      }
      
      // perform update of parameters from accumulators
      UpdateClusterParametersFromAccums(pModelSet);
    }
    
    // ///////
    // ordinary update    
    else 
    {
      // !!! MAP update should not be here !!!
      int vec_size    = mpVariance->VectorSize();
//      FLOAT *mean_vec = mpMean->mVector.pData();
//      FLOAT *var_vec  = mpVariance->mVector.pData();
  
      FLOAT* mmac  = mpMean->mpAccums;
      FLOAT* mnrm  = mpMean->mpAccums + 1 * vec_size;
      FLOAT* vvac  = mpVariance->mpAccums;
      FLOAT* vmac  = mpVariance->mpAccums + 1 * vec_size;
      FLOAT* vnrm  = mpVariance->mpAccums + 2 * vec_size;
  
      FLOAT* var_pri  = mpVariance->mpPrior->mVector.pData();
      FLOAT* mean_pri = mpMean    ->mpPrior->mVector.pData();
      
      if (pModelSet->mUpdateMask & UM_MAP)
      {
        for (i = 0; i < vec_size; i++) {
          mmac[i] += pModelSet->mMapTau * mean_pri[i];
          vmac[i] += pModelSet->mMapTau * mean_pri[i];
          vvac[i] += pModelSet->mMapTau * (1.0/var_pri[i] + SQR(mean_pri[i]));
        }
        *mnrm += pModelSet->mMapTau;
        *vnrm += pModelSet->mMapTau;
      }
      // !!! MAP update should not be here !!!
      
      if (!mpVariance->mpMacro)
        mpVariance->UpdateFromAccums(pModelSet);
  
      if (!mpMean->mpMacro)
        mpMean->UpdateFromAccums(pModelSet);
    }
    
    // Perform some variance flooring
    FloorVariance(pModelSet);
    
    // recompute the GConst for the Gaussian mixture
    ComputeGConst();
  }; // UpdateFromAccums(const ModelSet * pModelSet) 

  
  //***************************************************************************
  //***************************************************************************
  //virtual 
  void 
  Mixture::
  Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *pUserData)
  {
    int n = 0;
    char *chptr = NULL;
    
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_MIXTURE && mask & MTM_PRESCAN) {
      action(mt_mixture, nodeName, this, pUserData);
    }
  
    if (mask & MTM_MEAN && !mpMean->mpMacro) {
      if (n > 0) strncpy(chptr, ".mean", n);
      action(mt_mean, nodeName, mpMean, pUserData);
    }
  
    if (mask & MTM_VARIANCE && !mpVariance->mpMacro) {
      if (n > 0) strncpy(chptr, ".cov", n);
      action(mt_variance, nodeName, mpVariance, pUserData);
    }
  
    if (mask & MTM_XFORM_INSTANCE && mpInputXform &&
      !mpInputXform->mpMacro) {
      if (n > 0) strncpy(chptr, ".input", n);
      mpInputXform->Scan(mask, nodeName, action, pUserData);
    }
  
    if (mask & MTM_MIXTURE && !(mask & MTM_PRESCAN)) {
      if (n > 0) chptr = '\0';
      action(mt_mixture, nodeName, this, pUserData);
    }
  }

  
  //***************************************************************************
  //***************************************************************************
  void 
  Mixture::
  UpdateClusterParametersAccums()
  {
    // create a helping weight vector which is a concatenation of all elementary
    // weight vectors
    BasicVector<FLOAT> aux_vec(mAccumG.Rows());
    CatClusterWeightPartialVectors(aux_vec, mpMean->mpClusterWeightVectors, 
      mpMean->mNClusterWeightVectors);
    
    // add const * vec * vec^T
    mAccumG.AddCVVtMul(mPartialAccumG, aux_vec, aux_vec);
    mAccumK.AddCVVtMul(1, aux_vec, mPartialAccumK);
    
    // clear the partial accumulator as they are strictly speaker dependent
    mPartialAccumG = 0.0;
    mPartialAccumK.Clear();
    
    // for discriminative training:
    if (mAccumGd.IsInitialized())
    {
      mAccumGd.AddCVVtMul(mPartialAccumGd, aux_vec, aux_vec);      
      mPartialAccumGd = 0.0;
    }
  }    
          

  //****************************************************************************
  //****************************************************************************
  void 
  Mixture::
  UpdateClusterParametersFromAccums(const ModelSet * pModelSet)
  {
    if (pModelSet->mUpdateMask & UM_MEAN)
    {
      // update the mean matrix
      Matrix<FLOAT>   g_inv(mAccumG);
      
      g_inv.Invert();
      mpMean->mClusterMatrixT.RepMMMul(g_inv, mAccumK);
    }
      
    if (pModelSet->mUpdateMask & UM_VARIANCE)
    {
      FLOAT           tmp_val;
      //: KLUDGE:
      // norms for variances are stored in vector[vec_zize*3]. Not nice...
      for (size_t i = 0; i < mpVariance->VectorSize(); i++)
      {
        tmp_val = 0;
         
        for (size_t j = 0; j < mAccumK.Rows(); j++)
        {
          tmp_val += mpMean->mClusterMatrixT[j][i] * mAccumK[j][i]; 
        }
        
        // old vector
        //mpVariance->mpVectorO[i] = (mpVariance->mpVectorO[mpVariance->VectorSize()*3]) / 
        //  (mAccumL[i] - tmp_val);
        mpVariance->mVector[i] = (mpVariance->mpAccums[mpVariance->VectorSize()*2]) / 
          (mAccumL[i] - tmp_val);        
      }
    }
  }    
  
    
  //**************************************************************************  
  //**************************************************************************  
  // Mean section
  //**************************************************************************  
  //**************************************************************************  
  Mean::
  Mean(size_t vectorSize, int allocateAccums) :
    mVector(vectorSize)
  {
    size_t accum_size = 0;
    void*  free_vec;
    
    if (0 < allocateAccums)
    {
      // allocateAccums == 2 for MMI update
      accum_size = align<16>(((vectorSize + 1) * allocateAccums) * sizeof(FLOAT)); 
      
      mpAccums = static_cast<FLOAT*>
        (stk_memalign(16, accum_size, &free_vec));
#ifdef STK_MEMALIGN_MANUAL
      mpAccumsFree = static_cast<FLOAT*>(free_vec);
#endif
    }
    else
    {
      mpAccums = NULL;
#ifdef STK_MEMALIGN_MANUAL
      mpAccumsFree = NULL;
#endif
    }
    
    mpXformStatAccum          = NULL;
    mNumberOfXformStatAccums  = 0;
    mUpdatableFromStatAccums  = true;
    mpPrior                   = NULL;
    mpClusterWeightVectors    = NULL;
    mNClusterWeightVectors    = 0;
    mpOccProbAccums           = NULL;
  }  

  //**************************************************************************  
  //**************************************************************************  
  Mean::
  ~Mean()
  {
    if (mpXformStatAccum != NULL)
    {
      free(mpXformStatAccum->mpStats);
      free(mpXformStatAccum);
    }
    
    // delete refferences to weights Bias xforms
    if (NULL != mpClusterWeightVectors)
    { 
      delete [] mpClusterWeightVectors;
    }
    
    if (NULL != mpOccProbAccums)
    {
      delete [] mpOccProbAccums;
    }
    
#ifdef STK_MEMALIGN_MANUAL
    free(mpAccumsFree);
#else
    free(mpAccums);
#endif
  }  
  
  //**************************************************************************  
  //**************************************************************************  
  void
  Mean::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    size_t   i;
    //FLOAT*   vec  = mpVectorO;
    //FLOAT *  acc  = mpVectorO + 1 * VectorSize();
    //FLOAT    nrm  = mpVectorO  [2 * VectorSize()];
    FLOAT*   vec  = mVector.pData();
    FLOAT*   acc  = mpAccums;
    FLOAT    nrm  = !pModelSet->mModelUpdateDoesNotNormalize ? mpAccums[VectorSize()] : 1;
  
    if (pModelSet->mUpdateMask & UM_MEAN) 
    {
      if (mNumberOfXformStatAccums == 0) 
      {              
        for (i = 0; i < VectorSize(); i++) 
        {
          vec[i] = acc[i] / nrm;
        }
      } 
      else 
      { // Updating from xform statistics
        size_t r;
        size_t c;
        size_t pos = 0;
        
        if (!mUpdatableFromStatAccums) 
          return;
        
        // If 'mUpdatableFromStatAccums' is true then 'mNumberOfXformStatAccums' is equal to 1
        for (i=0; i < mNumberOfXformStatAccums; i++) 
        {
          LinearXform*  xform   = (LinearXform*) mpXformStatAccum[i].mpXform;
          size_t        in_size = xform->mInSize;
          FLOAT*        mnv     = mpXformStatAccum[i].mpStats;
  
          assert(xform->mXformType == XT_LINEAR);
          assert(pos + xform->mOutSize <= VectorSize());
          
          for (r = 0; r < xform->mOutSize; r++) 
          {
            //mpVectorO[pos + r] = 0;
            mVector[pos + r] = 0;
            
            for (c = 0; c < in_size; c++) 
            {
              //mpVectorO[pos + r] += mnv[c] * xform->mpMatrixO[in_size * r + c];
              //mpVectorO[pos + r] += mnv[c] * xform->mMatrix[r][c];
              mVector[pos + r] += mnv[c] * xform->mMatrix[r][c];
            }
          }
          pos += xform->mOutSize;
        }
        
        assert(pos == VectorSize());
      }
    }
  } // UpdateFromAccums(const ModelSet * pModelSet)
  
  
  //**************************************************************************  
  //**************************************************************************  
  void 
  Mean::
  RecalculateCAT()
  {
    if (NULL != mpClusterWeightVectors)
    {
      //:KLUDGE: optimize this
      //memset(mpVectorO, 0, sizeof(FLOAT) * VectorSize());
      mVector.Clear();
      
      int v = 0;
      
      // go through weights vectors
      for (size_t w = 0; w < mNClusterWeightVectors; w++)
      {
        // go through rows in the cluster matrix
        for (size_t i = 0; i < mpClusterWeightVectors[w]->mInSize; i++)
        {
          // go through cols in the cluster matrix
          for (size_t j = 0; j < VectorSize(); j++)
          {
            //mpVectorO[j] += mClusterMatrixT(v, j) * mpClusterWeightVectors[w]->mVector[0][i];
            mVector[j] += mClusterMatrixT[v][j] * mpClusterWeightVectors[w]->mVector[0][i];
          }
          // move to next cluster mean vector
          v++;
        }
      }
      
    }
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  // Variance section
  //**************************************************************************  
  //**************************************************************************  
  Variance::
  Variance(size_t vectorSize, int allocateAccums) :
    mVector(vectorSize)
  {
    void* free_vec;
    
    if (0 < allocateAccums) 
    {
      // allocateAccums == 2 for MMI update
      size_t accum_size = align<16>(((2*vectorSize + 1) * allocateAccums)*sizeof(FLOAT)); 
          
      mpAccums = static_cast<FLOAT*>
        (stk_memalign(16, accum_size, &free_vec));
#ifdef STK_MEMALIGN_MANUAL
      mpAccumsFree = static_cast<FLOAT*>(free_vec);
#endif
    }
    else
    {
      mpAccums = NULL;
#ifdef STK_MEMALIGN_MANUAL
      mpAccumsFree = NULL;
#endif
    }
    
    //mVectorSize               = vectorSize;    
    mpXformStatAccum          = NULL;
    mNumberOfXformStatAccums  = 0;
    mUpdatableFromStatAccums  = true;
    mpPrior                   = NULL;
  }  
  
  //**************************************************************************  
  //**************************************************************************  
  Variance::
  ~Variance()
  {
    if (mpXformStatAccum != NULL)
    {
      free(mpXformStatAccum->mpStats);
      free(mpXformStatAccum);
    }
    
#ifdef STK_MEMALIGN_MANUAL
//    free(mpVectorOFree);
#else
//    free(mpVectorO);
#endif

#ifdef STK_MEMALIGN_MANUAL
    free(mpAccumsFree);
#else
    free(mpAccums);
#endif
  }  
  
  //**************************************************************************  
  //**************************************************************************  
  void
  Variance::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    size_t    i;
  
    if (pModelSet->mUpdateMask & UM_VARIANCE) 
    {
      if (mNumberOfXformStatAccums == 0) 
      {
        // old vector
        //FLOAT *vec  = mpVectorO;
        FLOAT *vec  = mVector.pData();
        //FLOAT *vac  = mpVectorO + 1 * VectorSize(); // varriance accum 
        //FLOAT *mac  = mpVectorO + 2 * VectorSize(); // mean accum 
        //FLOAT *nrm  = mpVectorO + 3 * VectorSize(); // norm - occupation count
        FLOAT *vac  = mpAccums; // varriance accum 
        FLOAT *mac  = mpAccums + 1 * VectorSize(); // mean accum 
        FLOAT *nrm  = mpAccums + 2 * VectorSize(); // norm - occupation count
        
        for (i = 0; i < VectorSize(); i++) 
        {
          if (pModelSet->mUpdateMask & UM_OLDMEANVAR) 
          {
            vec[i] = *nrm / vac[i];
          } 
          else 
          {
            vec[i] = 1 / HIGHER_OF(0.0, vac[i] / *nrm - SQR(mac[i] / *nrm)); 
                      // HIGHER_OF is to avoid negative variances
          }
          
          /* This block is moved one level higher
          
          // !!! Need for transformation dependent varFloors
          if (pModelSet->mpVarFloor && 
              pModelSet->mpVarFloor->VectorSize() == VectorSize()) 
          {
            vec[i] = LOWER_OF(vec[i], pModelSet->mpVarFloor->mpVectorO[i]);
          }
          */
        }
      } 
      
      else // Updating from xform statistics
      { 
        size_t r;
        size_t c;
        size_t t;
        size_t pos = 0;
        
        if (!mUpdatableFromStatAccums) 
          return;
  
  //    If 'mUpdatableFromStatAccums' is true then 'mNumberOfXformStatAccums' is equal to 1
        for (i=0; i<mNumberOfXformStatAccums; i++) 
        {
          LinearXform*  xform = (LinearXform*) mpXformStatAccum[i].mpXform;
          size_t        in_size = xform->mInSize;
          FLOAT*        cov  = mpXformStatAccum[i].mpStats + in_size;
  
          assert(xform->mXformType == XT_LINEAR);
          assert(pos + xform->mOutSize <= VectorSize());
          
          for (r = 0; r < xform->mOutSize; r++) 
          {
            // old vector
            //mpVectorO[pos + r] = 0.0;
            mVector[pos + r] = 0.0;
  
            for (c = 0; c < in_size; c++) 
            {
              FLOAT aux = 0;
              for (t = 0; t <= c; t++) 
              {
                //aux += cov[c * (c+1)/2 + t] * xform->mpMatrixO[in_size * r + t];
                aux += cov[c * (c+1)/2 + t] * xform->mMatrix[r][t];
              }
  
              for (; t < in_size; t++) 
              {
                //aux += cov[t * (t+1)/2 + c] * xform->mpMatrixO[in_size * r + t];
                aux += cov[t * (t+1)/2 + c] * xform->mMatrix[r][t];
              }
              
              //mpVectorO[pos + r] += aux     * xform->mpMatrixO[in_size * r + c];
              // old vector
              //mpVectorO[pos + r] += aux     * xform->mMatrix[r][c];
              mVector[pos + r] += aux     * xform->mMatrix[r][c];              
            }
            
            // old vector
            //mpVectorO[pos + r] = 1 / HIGHER_OF(0.0, mpVectorO[pos + r]);
            mVector[pos + r] = 1 / HIGHER_OF(0.0, mVector[pos + r]);

            /* This block is moved one level higher
            if (pModelSet->mpVarFloor && 
                pModelSet->mpVarFloor->VectorSize() == VectorSize()) // !!! Need for transformation dependent varFloors
            {
                 mpVectorO[pos + r] = LOWER_OF(mpVectorO[pos + r], pModelSet->mpVarFloor->mpVectorO[pos + r]);
            }
            */
          } // for (r = 0; r < xform->mOutSize; r++) 
          
          pos += xform->mOutSize;
        } // for (i=0; i<mNumberOfXformStatAccums; i++) 
        
        assert(pos == VectorSize());
      }
    }
  }

  
  //**************************************************************************  
  //**************************************************************************  
  // Transition section
  //**************************************************************************  
  //**************************************************************************  
  Transition::
  Transition(size_t nStates, int allocateAccums)
  {
    size_t  alloc_size;
    alloc_size = 0 < allocateAccums ? 2 * SQR(nStates) + nStates: SQR(nStates);
    
    mpMatrixO = new FLOAT[alloc_size];
    mNStates = nStates;
    mpPrior  = NULL;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  Transition::
  ~Transition()
  {
    delete [] mpMatrixO;
  }


  //**************************************************************************  
  //**************************************************************************  
  void
  Transition::
  UpdateFromAccums(const ModelSet * pModelSet)
  {
    int       i;
    int       j;
    int       nstates = mNStates;
    FLOAT *   vec = mpMatrixO;
    FLOAT *   acc = mpMatrixO + 1 * SQR(mNStates);
  
    if (pModelSet->mUpdateMask & UM_TRANSITION) 
    {
      for (i=0; i < nstates; i++) 
      {
        FLOAT nrm = LOG_0;
        for (j=0; j < nstates; j++) 
        {
          LOG_INC(nrm, acc[i * nstates + j]);
        }
        
        if (nrm == LOG_0) 
          nrm = 0;
        
        for (j=0; j < nstates; j++) 
        {
          vec[i * nstates + j] = acc[i * nstates + j] - nrm; // it is in log
        }
      }
    }
  }


  //**************************************************************************  
  //**************************************************************************  
  void
  Transition::
  Print(std::ostream& rOstr)
  {
    FLOAT* p_vec = mpMatrixO;
    size_t i;
    size_t j;

    // go over all rows
    for (i=0; i < mNStates; i++) {
      // go over all cols
      for (j=0; j < mNStates; j++) {
        // print it
        rOstr << p_vec[i * mNStates + j];

        // print space for all but the last number
        if (j + 1 < mNStates) {
          rOstr << " ";
        }
      }

      // print EOL after the last number
      rOstr << std::endl;
    }
  }

  //**************************************************************************  
  //**************************************************************************  
  void 
  Transition::
  PrintCout()
  { this->Print(std::cout); }

  
  //**************************************************************************  
  //**************************************************************************  
  // Transition section
  //**************************************************************************  
  //**************************************************************************  
  
  //***************************************************************************
  //***************************************************************************
  State::
  State(size_t numMixtures)
  {
    mpMixture = (numMixtures > 0) ? new MixtureLink[numMixtures] : NULL;
    mpPrior   = NULL;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  State::
  ~State()
  {
    if (NULL != mpMixture)
      delete [] mpMixture;
    
    //std::cerr << "State destructor" << std::endl;      
  }

  
  //***************************************************************************
  //***************************************************************************
  void 
  State::
  UpdateFromAccums(const ModelSet * pModelSet, const Hmm * pHmm) 
  {
    size_t i, k;
  
    if (mOutPdfKind == KID_DiagC) 
    {
      FLOAT accum_sum = 0;
  
//    if (hmm_set->mUpdateMask & UM_WEIGHT) {

      for (i = 0; i < mNMixtures; i++)
      {
        if(UT_EBW == pModelSet->mUpdateType || (pModelSet->mMapTau > 0.0 && pModelSet->mUpdateMask & UM_MAP))
        {
          // For MMI update from single combined num-den accumulator, we can use numerator counts.
          ///accum_sum += mpMixture[i].mWeightAccumDen;

          //... but it does not make sense for MPE (we would have to distinguish MMI and MPE update),
          //!!! So we do not have mWeightAccum and prior weights are taken instead.
          //!!! Anyway, this is something not very nice.
          if(mpPrior->mpMixture[i].mWeight < LOG_MIN)
          {
            mpMixture[i].mWeightAccum = 0.0;
          }
          else
          {
            mpMixture[i].mWeightAccum = my_exp(mpPrior->mpMixture[i].mWeight);
          }
        }
/*
        // It can be useful to update prior model with -w 0.0, which avoid discarding Gaussians
        // with with zero (or very low) occupation. This way, we can make sure that the Gaussians
        // in the prior and the model that is being updated here will match. Now it is the time
        // to discard such Gaussians that have have LOG_0 (or very low) weight in the prior model.        
        if(pModelSet->mMapTau > 0.0
        && (mpPrior->mpMixture[i].mWeight < LOG_MIN 
           || exp(mpPrior->mpMixture[i].mWeight) < pModelSet->mMinMixWeight)) {
          mpMixture[i].mWeightAccum = 0.0;
        }
*/
        accum_sum += mpMixture[i].mWeightAccum;
      }

      if (accum_sum <= 0.0 && pModelSet->mMinMixWeight > 0.0) 
      {
        if (mpMacro) 
        {
          Warning("No occupation of '%s', state is not updated",
                  mpMacro->mpName);
        } 
        else 
        {
          size_t j; // find the state number
          for (j=0; j < pHmm->mNStates && pHmm->mpState[j] != this; j++)
            ;
          Warning("No occupation of '%s[%d]', state is not updated",
                  pHmm->mpMacro->mpName, (int) j + 1);
        }
        return;
      }
  
        // Remove mixtures with low weight
//      if (pModelSet->mUpdateMask & UM_WEIGHT) 
      {

//!!! In MPE update code, we do not update mixture with low counts
//!!! Here, we remove mixture with low weight, but even if weights are not to be re-estimated
//!!! We should behave more consistently

        for (i = 0, k = 0; i < mNMixtures; i++, k++) 
        {
          if (mpMixture[i].mWeightAccum / accum_sum < pModelSet->mMinMixWeight) 
          {
            if (pHmm)
            {
              size_t j;
              for (j=0; j < pHmm->mNStates && pHmm->mpState[j] != this; j++)
                ; // find the state number
              Warning("Discarding mixture %d of Hmm %s state %d because of too low mixture weight",
                      (int) k, pHmm->mpMacro->mpName, (int) j + 1);
            } 
            else 
            {
              assert(mpMacro);
              Warning("Discarding mixture %d of state %s because of too low mixture weight",
                      (int) k, mpMacro->mpName);
            }
            accum_sum -= mpMixture[i].mWeightAccum;
  
            if (!mpMixture[i].mpEstimates->mpMacro) 
            {
              mpMixture[i].mpEstimates->Scan(MTM_ALL,NULL,ReleaseItem,NULL);
            }
  
            mpMixture[i--] = mpMixture[--mNMixtures];
            continue;
          }
        }
      }
  
      // in case we don't want to normalize, reset the normalizing constant
      // this is probably not the optimal way to do this, but the easiest
      if (pModelSet->mModelUpdateDoesNotNormalize) {
        accum_sum = 1;
      }

      for (i = 0; i < mNMixtures; i++) 
      {
  //      printf("Weight Acc: %f\n", (float) state->mpMixture[i].mWeightAccum);
        if (pModelSet->mUpdateMask & UM_WEIGHT) 
        {
          mpMixture[i].mWeight = my_log(mpMixture[i].mWeightAccum / accum_sum);
        }
        
        if (!mpMixture[i].mpEstimates->mpMacro) 
        {
        //!!! This is just too ugly hack
          gWeightAccumDen = mpMixture[i].mWeightAccumDen;
          mpMixture[i].mpEstimates->UpdateFromAccums(pModelSet);
        }
      }
    }
  }
  

  //***************************************************************************
  //***************************************************************************
  //virtual 
  void 
  State::
  Scan(int mask, HMMSetNodeName nodeName,  ScanAction action, void *pUserData)
  {
    size_t    i;
    size_t    n = 0;
    char *    chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_STATE && mask & MTM_PRESCAN) 
    {
      action(mt_state, nodeName, this, pUserData);
    }
  
    if (mOutPdfKind != KID_PDFObsVec &&
      mask & (MTM_ALL & ~(MTM_STATE | MTM_HMM | MTM_TRANSITION))) 
    {
      for (i=0; i < mNMixtures; i++) 
      {
        if (!mpMixture[i].mpEstimates->mpMacro) 
        {
          if (n > 0 ) snprintf(chptr, n, ".mix[%d]", (int) i+1);
          mpMixture[i].mpEstimates->Scan(mask, nodeName,
                      action, pUserData);
        }
      }
    }
    if (mask & MTM_STATE && !(mask & MTM_PRESCAN)) 
    {
      if (n > 0) chptr = '\0';
      action(mt_state, nodeName, this, pUserData);
    }
  }
  
  
  //***************************************************************************
  //***************************************************************************
  // XformInstance seciton
  //***************************************************************************
  //***************************************************************************
  XformInstance::
  XformInstance(size_t vectorSize) :
    mOutputVector(vectorSize)
  {
    mpVarFloor       = NULL;
    mpXformStatCache = NULL;

    mpInput          = NULL;
    mpXform          = NULL;
    mTime            = 0;
    mpNext           = NULL;
    mNumberOfXformStatCaches = 0;
    mStatCacheTime   = 0;
    mpMemory         = NULL;
    mTotalDelay      = 0;
  }
  
  //***************************************************************************
  //***************************************************************************
  XformInstance::
  XformInstance(Xform* pXform, size_t vectorSize) :
    mOutputVector(vectorSize)
  {
    mpVarFloor          = NULL;
    mpXformStatCache    = NULL;

    mpInput             = NULL;
    mpXform             = pXform;
    mTime               = 0;
    mpNext              = NULL;
    mNumberOfXformStatCaches = 0;
    mStatCacheTime      = 0;
    mpMemory            = NULL;
    mTotalDelay         = 0;
  }
    
  
  //***************************************************************************
  //***************************************************************************
  XformInstance::
  ~XformInstance()
  {
    if (mpXformStatCache != NULL)
    {
      std::map<void*,int>::iterator i = mStatCachePtrs.find(mpXformStatCache);
      if (i == mStatCachePtrs.end()) {
        mStatCachePtrs[mpXformStatCache] = 1;
      }
      else {
        ++(i->second);
      }
        
      free(mpXformStatCache->mpStats);
      free(mpXformStatCache);

      mpXformStatCache = NULL;
    }
    //delete [] mpOutputVector;
    delete [] mpMemory;
  }
  
  
  //*****************************************************************************
  //*****************************************************************************
  //virtual   
  void 
  XformInstance::
  Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *pUserData)
  {
    int n = 0;
    char *chptr = NULL;
    
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_XFORM_INSTANCE && mask & MTM_PRESCAN) {
      action(mt_XformInstance, nodeName, this, pUserData);
    }
    
    if (mpInput != NULL) {
      if (!mpInput->mpMacro) {
        if (n > 0) strncpy(chptr, ".input", n);
          mpInput->Scan(mask, nodeName, action, pUserData);
      }
    }
  
    if (mask & MTM_XFORM && mpXform != NULL && 
      !mpXform->mpMacro) {
      if (n > 0) strncpy(chptr, ".mpXform", n);
      mpXform->Scan(mask, nodeName, action, pUserData);
    }
    
    if (mask & MTM_XFORM_INSTANCE && !(mask & MTM_PRESCAN)) {
      if (n > 0) chptr = '\0';
      action(mt_XformInstance, nodeName, this, pUserData);
    }
  }
  

  //***************************************************************************
  //***************************************************************************
  // Xform seciton
  //***************************************************************************
  //***************************************************************************
  
  //****************************************************************************
  //****************************************************************************
  //virtual 
  void 
  Xform::
  Scan(int mask, HMMSetNodeName nodeName, ScanAction action, void *pUserData)
  {
    size_t      n = 0;
    char *      chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mask & MTM_PRESCAN) 
    {
      action(mt_Xform, nodeName, this, pUserData);
    }
    
    if (mXformType == XT_COMPOSITE) 
    {
      CompositeXform *  cxf = dynamic_cast<CompositeXform *> (this);
      size_t            i;
      size_t            j;
  
      for (i=0; i < cxf->mNLayers; i++) 
      {
        for (j = 0; j < cxf->mpLayer[i].mNBlocks; j++) 
        {
          if (!cxf->mpLayer[i].mpBlock[j]->mpMacro) 
          {
            if (n > 0) 
              snprintf(chptr, n, ".part[%d,%d]", (int) i+1, (int) j+1);

            cxf->mpLayer[i].mpBlock[j]->Scan(mask, nodeName, action, pUserData);
          }
        }
      }
    } 
    else if (mXformType == XT_FEATURE_MAPPING) 
    {
      static_cast<FeatureMappingXform*>(this)->mpStateTo->Scan(mask, nodeName, action, pUserData);
      static_cast<FeatureMappingXform*>(this)->mpStateFrom->Scan(mask, nodeName, action, pUserData);
    }
    else if (mXformType == XT_GMM_POSTERIORS) 
    {
      static_cast<GmmPosteriorsXform*>(this)->mpState->Scan(mask, nodeName, action, pUserData);
    }
    else if (mXformType == XT_REGIONDEPENDENT) {
      RegionDependentXform*  rdxf = dynamic_cast<RegionDependentXform *> (this);
      size_t j;
  
      for (j = 0; j < rdxf->mNBlocks; j++) {
        if (!rdxf->mpBlock[j]->mpMacro) {
          if (n > 0) {
            snprintf(chptr, n, ".part[%d]", (int) j+1);
          }
          rdxf->mpBlock[j]->Scan(mask, nodeName, action, pUserData);
        }
      }
    }
  
    if (!(mask & MTM_PRESCAN)) 
    {
      if (n > 0) 
        chptr = '\0';
        
      action(mt_Xform, nodeName, this, pUserData);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void Xform::ExpandPredef(int macro_type, HMMSetNodeName nodeName,
                 MacroData * pData, void *pUserData) 
  {
    Xform *p_xform = static_cast<Xform *>(pData);
        
    switch(p_xform->mXformType)
    {
      case XT_LINEAR:
      {
        LinearXform *p_lxform = static_cast<LinearXform *>(p_xform);
        p_lxform->mPredefinedID = PLTID_NONE;
        break;
      }
      case XT_WINDOW:
      {
        WindowXform *p_wxform = static_cast<WindowXform *>(p_xform);
        p_wxform->mUsePredefVector = false;
        break;
      }
      case XT_BIAS:
      {
        BiasXform *p_bxform = static_cast<BiasXform *>(p_xform);
        p_bxform->mUsePredefVector = false;
        break;        
      }
    }
  }
      
    
  //**************************************************************************  
  //**************************************************************************  
  // XformLayer section
  //**************************************************************************  
  //**************************************************************************  
  XformLayer::
  XformLayer() : mOutputVector(), mNBlocks(0), mpBlock(NULL) 
  {
    // we don't do anything here, all is done by the consructor call
  }
  
  //**************************************************************************  
  //**************************************************************************  
  XformLayer::
  ~XformLayer()
  {
    //delete the pointers
    //delete [] mpOutputVector;
    delete [] mpBlock;
  }

  //**************************************************************************  
  //**************************************************************************  
  Xform **
  XformLayer::
  InitBlocks(size_t nBlocks)
  {
    mpBlock  = new Xform*[nBlocks];
    mNBlocks = nBlocks;
    
    for (size_t i = 0; i < nBlocks; i++) 
      mpBlock[i] = NULL;
      
    return mpBlock;
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  // CompositeXform section
  //**************************************************************************  
  //**************************************************************************  
  CompositeXform::
  CompositeXform(size_t nLayers)
  {
    mpLayer = new XformLayer[nLayers];
    
    // set some properties
    mMemorySize = 0;
    mDelay      = 0;
    mNLayers    = nLayers;
    mXformType = XT_COMPOSITE;    
  }
  
  //**************************************************************************  
  //**************************************************************************  
  CompositeXform::
  ~CompositeXform()
  {
    delete [] mpLayer;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  CompositeXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    size_t  i;
    size_t  j;
  
    for (i = 0; i < mNLayers; i++) 
    {
      //FLOAT* in  = (i == 0)          ? pInputVector  : mpLayer[i-1].mpOutputVector;
      //FLOAT* out = (i == mNLayers-1) ? pOutputVector : mpLayer[i].mpOutputVector;
      FLOAT* in  = (i == 0)          ? pInputVector  : mpLayer[i-1].mOutputVector.pData();
      FLOAT* out = (i == mNLayers-1) ? pOutputVector : mpLayer[i].mOutputVector.pData();
      
      for (j = 0; j < mpLayer[i].mNBlocks; j++) 
      {
        mpLayer[i].mpBlock[j]->Evaluate(in, out, pMemory, direction);
        in  += mpLayer[i].mpBlock[j]->mInSize;
        out += mpLayer[i].mpBlock[j]->mOutSize;
        pMemory += mpLayer[i].mpBlock[j]->mMemorySize;
      }
    }
    
    return pOutputVector;
  }; //Evaluate(...)



  //**************************************************************************  
  //**************************************************************************  
  // RegionDependentXform section
  //**************************************************************************  
  //**************************************************************************  
  RegionDependentXform::
  RegionDependentXform(size_t nBlocks):
    mBlockOutputVector()
  {
    mpBlock  = new Xform*[nBlocks];
    mNBlocks = nBlocks;
    
    for (size_t i = 0; i < nBlocks; i++) 
      mpBlock[i] = NULL;
      
    // set some properties
    mMemorySize = 0;
    mDelay      = 0;
    mNBlocks    = nBlocks;
    mXformType  = XT_REGIONDEPENDENT;    
  }
  
  //**************************************************************************  
  //**************************************************************************  
  RegionDependentXform::
  ~RegionDependentXform()
  {
    delete [] mpBlock;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  RegionDependentXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    size_t  i;
    size_t  j;
    FLOAT* region_weights = pInputVector;
    FLOAT* input_vector   = pInputVector + mNBlocks;
    
    for(j = 0; j < mBlockOutputVector.Length(); j++) {
      pOutputVector[j] = 0;
    }  
    for (i = 0; i < mNBlocks; i++) 
    {
      if(pInputVector[i] != 0 || mpBlock[i]->mMemorySize) {
        mpBlock[i]->Evaluate(
          input_vector,
          mBlockOutputVector.pData(), 
          pMemory, 
          direction);
          
        for(j = 0; j < mpBlock[i]->mOutSize; j++) {
          pOutputVector[j] = pOutputVector[j] + region_weights[i] * mBlockOutputVector[j];
        }
        pMemory += mpBlock[i]->mMemorySize;
      }
    }
    
    return pOutputVector;
  }; //Evaluate(...)

  
  
  //**************************************************************************  
  //**************************************************************************  
  // LinearXform section
  //**************************************************************************  
  //**************************************************************************  
  LinearXform::
  LinearXform(size_t inSize, size_t outSize):
    mMatrix(outSize, inSize)
  {
    mOutSize      = outSize;
    mInSize       = inSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_LINEAR;
    mPredefinedID = PLTID_NONE;
  }
    
  
  //**************************************************************************  
  //**************************************************************************  
  LinearXform::
  ~LinearXform()
  {
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT*
  LinearXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    // make the 1x1 matrices multiply fast!!!
    if(mInSize == 1  &&  mOutSize == 1) {
      pOutputVector[0] = *(mMatrix[0]) * pInputVector[0];
      return pOutputVector;
    }

#ifdef HAVE_ATLAS
    memset(pOutputVector, 0, mOutSize*sizeof(FLOAT));
# if DOUBLEPRECISION
    cblas_dgemv(CblasRowMajor, CblasNoTrans, mMatrix.Rows(), mMatrix.Cols(), 1.0, 
        mMatrix.pData(), mMatrix.Stride(), pInputVector, 1, 1.0F, pOutputVector, 1);
# else
    cblas_sgemv(CblasRowMajor, CblasNoTrans, mMatrix.Rows(), mMatrix.Cols(), 1.0, 
        mMatrix.pData(), mMatrix.Stride(), pInputVector, 1, 1.0F, pOutputVector, 1);
# endif
#else

    size_t c; // column counter
    size_t r; // row counter
    
    // matrix multiplication
    for (r = 0; r < mOutSize; r++) 
    {
      pOutputVector[r] = 0.0;
      for (c = 0; c < mInSize; c++) 
      {
        pOutputVector[r] += mMatrix[r][c] * pInputVector[c];
      }
    }
#endif
    return pOutputVector;  
  }; //Evaluate(...)
  

  //**************************************************************************  
  //**************************************************************************  
  // BiasXform section
  //**************************************************************************  
  //**************************************************************************  
  BiasXform::
  BiasXform(size_t vectorSize) :
    mVector(1, vectorSize)
  {
    mInSize       = vectorSize;
    mOutSize      = vectorSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_BIAS;
    mUsePredefVector = false;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  BiasXform::
  ~BiasXform()
  {
    //delete [] mpVectorO;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  BiasXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[i] + mVector[0][i];
    }
    return pOutputVector;
  }; //Evaluate(...)
  

  //**************************************************************************  
  //**************************************************************************  
  // FuncXform section
  //**************************************************************************  
  //**************************************************************************  
  FuncXform::
  FuncXform(size_t size, int funcId)
  {
    mFuncId       = funcId;
    mInSize       = size;
    mOutSize      = size;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_FUNC;
  }  
  
  //**************************************************************************  
  //**************************************************************************  
  FuncXform::
  ~FuncXform()
  {}
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  FuncXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    gFuncTable[mFuncId].funcPtr(pInputVector, pOutputVector, mOutSize);
    return pOutputVector;
  }; //Evaluate(...)
  
      
  //**************************************************************************  
  //**************************************************************************  
  // CopyXform section
  //**************************************************************************  
  //**************************************************************************  
  CopyXform::
  CopyXform(size_t inSize, size_t outSize)
  {
    mpIndices     = new int[outSize];
    mOutSize      = outSize;
    mInSize       = inSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_COPY;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  CopyXform::
  ~ CopyXform()
  {
    delete [] mpIndices;
  }

  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  CopyXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[mpIndices[i]];
    }
    return pOutputVector;
  }; //Evaluate(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  // BlockCopyXform section
  //**************************************************************************  
  //**************************************************************************  
  BlockCopyXform::
  BlockCopyXform(size_t nBlocks)
  {
    mpIndices     = 0;
    mpBlocks      = 0;
    mNBlocks      = 0;
    mNRows        = 0;
    mOutSize      = 0;
    mInSize       = 0;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_BLOCKCOPY;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  BlockCopyXform::
  ~ BlockCopyXform()
  {
    if(mpIndices)
    {
      delete [] mpIndices;
    }
    if(mpBlocks)
    {
      delete [] mpBlocks;
    }
  }  
    
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  BlockCopyXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[mpIndices[i]];
    }
    return pOutputVector;
  }; //Evaluate(...)


  //**************************************************************************  
  //**************************************************************************  
  // TransposeXform section
  //**************************************************************************  
  //**************************************************************************  
  TransposeXform::
  TransposeXform(size_t inRows, size_t inCols)
  {
    mInRows       = inRows;
    mInSize       = inRows * inCols;
    mOutSize      = mInSize;
    mpIndices     = new int[mOutSize];
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_TRANSPOSE;
    
    size_t i;
    size_t j;
    size_t k = 0;
    for(i = 0; i < inCols; i++)
    {
      for(j = 0; j < inRows; j++)
      {
         mpIndices[k] = j * inCols + i;
	 k++;      
      }
    }
  }
  
  //**************************************************************************  
  //**************************************************************************  
  // ConstantXform section
  //**************************************************************************  
  //**************************************************************************  
  ConstantXform::
  ConstantXform(size_t outSize) :
    mVector(1, outSize)
  {
    mInSize       = 0;
    mOutSize      = outSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_CONSTANT;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  TransposeXform::
  ~ TransposeXform()
  {
    delete [] mpIndices;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  TransposeXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[mpIndices[i]];
    }
    return pOutputVector;
  }; //Evaluate(...)

  
  //**************************************************************************  
  //**************************************************************************  
  // WindowXform section
  //**************************************************************************  
  //**************************************************************************  
  WindowXform::
  WindowXform(size_t vectorSize) :
     mVector(vectorSize)
  {
    mOutSize      = vectorSize;
    mInSize       = vectorSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_WINDOW;
    mUsePredefVector = false;
  }
  
  //**************************************************************************  
  //**************************************************************************  
  WindowXform::
  ~ WindowXform()
  {}
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  WindowXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = pInputVector[i] * mVector[i];
    }
    return pOutputVector;
  }; //Evaluate(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  ConstantXform::
  ~ConstantXform()
  {}
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  ConstantXform::
  Evaluate(FLOAT *    pInputVector, 
           FLOAT *    pOutputVector,
           char *     pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    
    for (i = 0; i < mOutSize; i++) 
    {
      pOutputVector[i] = mVector[0][i];
    }
    return pOutputVector;
  }; //Evaluate(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  // FeatureMappingXform section
  //**************************************************************************  
  //**************************************************************************  
  FeatureMappingXform::
  FeatureMappingXform(size_t inSize)
  {
    mOutSize      = inSize;
    mInSize       = inSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_FEATURE_MAPPING;
    mpStateFrom   = NULL;
    mpStateTo     = NULL;
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  // virutal
  FeatureMappingXform::
  ~FeatureMappingXform()
  {
    //std::cerr << "FeatureMappingXform destructor" << std::endl;
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  FeatureMappingXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    State*  state_from = mpStateFrom;
    State*  state_to   = mpStateTo;
    size_t  i;
    size_t  max_i      = 0;
    FLOAT   g_like;
    FLOAT   max_g_like = LOG_0;
    
    // find best input mixture
    for (i = 0; i < state_from->mNMixtures; i++)
    {
      g_like = Decoder<DecoderNetwork>::DiagCGaussianDensity(state_from->mpMixture[i].mpEstimates, 
        pInputVector, NULL) + state_from->mpMixture[i].mWeight;
        
      if (g_like > max_g_like)
      {
        max_i      = i;
        max_g_like = g_like;
      }
    }
    
    // retrieve matching output mixture
    const FLOAT*  m_from  = state_from->mpMixture[max_i].mpEstimates->mpMean->mVector.cpData();
    const FLOAT*  v_from  = state_from->mpMixture[max_i].mpEstimates->mpVariance->mVector.cpData();
    const FLOAT*  m_to    = state_to->mpMixture[max_i].mpEstimates->mpMean->mVector.cpData();
    const FLOAT*  v_to    = state_to->mpMixture[max_i].mpEstimates->mpVariance->mVector.cpData();    
    
    // compute new feature vector
    // compute each element separately
#ifndef OPTIMIZE_LOOPS_FOR_SSE    
    for (i = 0; i < mInSize; i++)
    {
      pOutputVector[i] = (pInputVector[i] - m_from[i]) * v_to[i] / v_from[i] + m_to[i];
    }
#else
          f4vector* pov      = reinterpret_cast<f4vector*>(pOutputVector);
    const f4vector* piv      = reinterpret_cast<const f4vector*>(pInputVector);
    const f4vector* pm_from  = reinterpret_cast<const f4vector*>(m_from);
    const f4vector* pv_from  = reinterpret_cast<const f4vector*>(v_from);
    const f4vector* pm_to    = reinterpret_cast<const f4vector*>(m_to);
    const f4vector* pv_to    = reinterpret_cast<const f4vector*>(v_to);
    
    for (i = 0; i < mInputVectorSize; i += 4)
    {
      pov->v = (piv->v - pm_from->v) * pv_to->v / pv_from->v + pm_to->v;

      pov     ++;
      piv     ++;
      pm_from ++;
      pv_from ++;
      pm_to   ++;
      pv_to   ++;      
    }
#endif
    
    return pOutputVector;
  }; //Evaluate(...)
    
  
  //**************************************************************************  
  //**************************************************************************  
  // GmmPosteriorsXform section
  //**************************************************************************  
  //**************************************************************************  
  GmmPosteriorsXform::
  GmmPosteriorsXform(size_t inSize):
    mNBestLikesVector()
  {
    mInSize       = inSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mXformType    = XT_GMM_POSTERIORS;
    mpState       = NULL;
    mScale        = 0.0;
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  // virutal
  GmmPosteriorsXform::
  ~GmmPosteriorsXform()
  {  
    
    //std::cerr << "GmmPosteriorsXform destructor" << std::endl;
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  GmmPosteriorsXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    size_t  i;
    FLOAT min_like;
    FLOAT like_sum = LOG_0;
    
    // find best input mixture
    for (i = 0; i < mpState->mNMixtures; i++) {
      pOutputVector[i] = (Decoder<DecoderNetwork>::DiagCGaussianDensity(mpState->mpMixture[i].mpEstimates, 
        pInputVector, NULL) + mpState->mpMixture[i].mWeight) * mScale;
    }
    partial_sort_copy(pOutputVector, pOutputVector+mpState->mNMixtures, 
                      mNBestLikesVector.begin(), mNBestLikesVector.end(), std::greater<FLOAT>());
    min_like = mNBestLikesVector.back();
    
    for (i = 0; i < mpState->mNMixtures; i++) {
      if(pOutputVector[i] >= min_like) {
        like_sum = LogAdd(like_sum, pOutputVector[i]);
      }
    }
    for (i = 0; i < mpState->mNMixtures; i++) {
      if(pOutputVector[i] >= min_like) {
        pOutputVector[i] = my_exp(pOutputVector[i]-like_sum);
      } else {
        pOutputVector[i] = 0;
      }
    }
    return pOutputVector;
  }; //Evaluate(...)
    
  
  //**************************************************************************  
  //**************************************************************************  
  // FrantaProductXform section
  //**************************************************************************  
  //**************************************************************************  
  FrantaProductXform::
  FrantaProductXform(size_t inSize, size_t n_parts)
  {
    mOutSize      = inSize / n_parts;
    mInSize       = inSize;
    mMemorySize   = 0;
    mDelay        = 0;
    mNParts       = n_parts;
    mXformType    = XT_FRANTA_PRODUCT;
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  FrantaProductXform::
  ~ FrantaProductXform()
  {
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT*
  FrantaProductXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    size_t i;
    size_t j;
    
    for (i = 0; i < mOutSize; i++)
    {
      pOutputVector[i] = 1.0;
      
      for (j = 0; j < mNParts; j++)
      {
        pOutputVector[i] *= pInputVector[j * mOutSize + i];
      }
    }
    
    return pOutputVector;
  }; //Evaluate(...)

  
    

  //**************************************************************************  
  //**************************************************************************  
  // MatlabXform section
  //**************************************************************************  
  //**************************************************************************  
  MatlabXform::
  MatlabXform(size_t inSize, size_t outSize, size_t delay)
  {
    mOutSize      = outSize;
    mInSize       = inSize;
    mMemorySize   = 0;
    mDelay        = delay;
    mXformType    = XT_MATLAB;
#ifdef MATLAB_ENGINE
#  ifdef DOUBLEPRECISION    
    mpInput       = mxCreateNumericMatrix(1, inSize,  mxSINGLE_CLASS, mxREAL);
#  else
    mpInput       = mxCreateNumericMatrix(1, inSize,  mxDOUBLE_CLASS, mxREAL);
#  endif
    mpOutput      = NULL;
    mpEp          = NULL;
#endif
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  MatlabXform::
  ~ MatlabXform()
  {
#ifdef MATLAB_ENGINE
    if (NULL != mpInput)  mxDestroyArray(mpInput);
    if (NULL != mpEp)     engClose(mpEp);
#endif  
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT*
  MatlabXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
#ifdef MATLAB_ENGINE
    // check whether we have a running instance of Matlab. If not, create one
    if (NULL == mpEp)
    {
      if (!(mpEp = engOpen("\0"))) 
      {
        Error("Can't start MATLAB Engine");
      }      
    }
    
    // copy our data to the MATLAB Matrix struct
    //mxSetData(mpInput, static_cast<void*>(pInputVector));
    memcpy(mxGetData(mpInput), static_cast<void*>(pInputVector), 
      mInSize * sizeof(FLOAT));
    
    engPutVariable(mpEp, "STKInput", mpInput);
    
    // run the program
    engEvalString(mpEp, mProgram.c_str());
    
    if (NULL == (mpOutput = engGetVariable(mpEp, "STKOutput")))
      Error("Cannot retrieve STKOutput (perhaps an error in Matlab script)");
    
    // copy result to pOutputVector
    memcpy(static_cast<void*>(pOutputVector), mxGetData(mpOutput), 
      mOutSize * sizeof(FLOAT));
      
    mxDestroyArray(mpOutput);
#endif    

    return pOutputVector;
  }; //Evaluate(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  // StackingXform section
  //**************************************************************************  
  //**************************************************************************  
  StackingXform::
  StackingXform(size_t stackSize, size_t inSize)
  {
    size_t out_size = stackSize * inSize;
  
    mOutSize      = out_size;
    mInSize       = inSize;
    mMemorySize   = out_size * sizeof(FLOAT);
    mDelay        = stackSize - 1;
    mHorizStack   = 0;
    mXformType    = XT_STACKING;      
  }
  
  //**************************************************************************  
  //**************************************************************************  
  StackingXform::
  ~StackingXform()
  { }
  
  //**************************************************************************  
  //**************************************************************************  
  //virtual 
  FLOAT *
  StackingXform::
  Evaluate(FLOAT*     pInputVector, 
           FLOAT*     pOutputVector,
           char*      pMemory,
           PropagDirectionType  direction)
  {
    FLOAT* stack = reinterpret_cast <FLOAT*> (pMemory);
  
    memmove(stack + (direction == BACKWARD ? mInSize : 0),
            stack + (direction ==  FORWARD ? mInSize : 0),
            (mOutSize - mInSize) * sizeof(FLOAT));
    
    memmove(stack + (direction ==  FORWARD ? mOutSize - mInSize : 0),
            pInputVector, mInSize * sizeof(FLOAT));
  
    if (!mHorizStack) 
    {
      memmove(pOutputVector, stack, mOutSize * sizeof(FLOAT));
    } 
    else 
    { // stacking == HORZ_STACK
      size_t t;
      size_t c;
      size_t stack_size = mOutSize / mInSize;
      
      for (t = 0; t < stack_size; t++) 
      {
        for (c = 0; c < mInSize; c++) 
        {
          pOutputVector[c * stack_size + t] = stack[t * mInSize + c];
        }
      }
    }
    
    return pOutputVector;
  }; //Evaluate(...)
    
  
 
  //**************************************************************************
  //**************************************************************************
  //   ModelSet section
  //**************************************************************************
  //**************************************************************************
  
  //**************************************************************************
  //**************************************************************************  
  void
  ModelSet::
  Init(FlagType flags)
  {
    if (!my_hcreate_r(100, &mHmmHash)            ||
        !my_hcreate_r(100, &mStateHash)          ||
        !my_hcreate_r( 10, &mMixtureHash)        ||
        !my_hcreate_r( 10, &mMeanHash)           ||
        !my_hcreate_r( 10, &mVarianceHash)       ||
        !my_hcreate_r( 10, &mTransitionHash)     ||
        !my_hcreate_r( 10, &mXformInstanceHash)  ||
        !my_hcreate_r( 10, &mXformHash)) 
    {
      Error("Insufficient memory");
    }
  
    /// :TODO: 
    /// Get rid of this by making STK use Matrix instead of FLOAT arrays
    mUseNewMatrix               = false;
    
    this->mpXformInstances      = NULL;
    this->mpInputXform          = NULL;
    this->mpFirstMacro          = NULL;
    this->mpLastMacro           = NULL;
    this->mInputVectorSize      = -1;
    this->mInputVectorStride    = 0;
    this->mParamKind            = -1;
    this->mOutPdfKind           = KID_UNSET;
    this->mDurKind              = KID_UNSET;
    this->mNMixtures            = 0;
    this->mNStates              = 0;
    this->mAllocAccums          = flags & MODEL_SET_WITH_ACCUM         ? 1 :
                                  flags & MODEL_SET_WITH_TWO_ACCUM_SET ? 2 : 0;
    this->mTotalDelay           = 0;
    InitLogMath();
    
    this->mSaveGlobOpts         = true;
  
    //Reestimation params
    this->mMinOccurances        = 3;
    this->mMinMixWeight         = MIN_WEGIHT;
    this->mpVarFloor            = NULL;
    this->mMinVariance          = 0.0;
    this->mUpdateMask           = UM_TRANSITION | UM_MEAN | UM_VARIANCE |
                                  UM_WEIGHT | UM_XFSTATS | UM_XFORM;
    this->mpXformToUpdate       = NULL;
    this->mNumberOfXformsToUpdate     = 0;
    this->mGaussLvl2ModelReest  = 0;
    this->mUpdateType           = UT_ML;
    this->mISmoothAfterD        = false;
    this->JSmoothing            = false;
    this->MMI_E                 = 2.0;
    this->MMI_h                 = 2.0;
    this->MMI_tauI              = 100.0;
    this->mMapTau               = 10.0;
    this->mMinOccupation        = 0.0;
    this->mISmoothingMaxOccup   = -1.0;
    this->mModelUpdateDoesNotNormalize = false;
    mpClusterWeightVectors            = NULL;
    mNClusterWeightVectors            = 0;
    mpGw                              = NULL;
    mpKw                              = NULL;
    //mClusterParametersUpdate    = false;
    InitKwdTable();    
  } // Init(...);


  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  Release()
  {
    size_t i;  
    
    Scan(MTM_REVERSE_PASS | MTM_ALL, NULL, ReleaseItem, NULL);
  
    ReleaseMacroHash(&mHmmHash);
    ReleaseMacroHash(&mStateHash);
    ReleaseMacroHash(&mMixtureHash);
    ReleaseMacroHash(&mMeanHash);
    ReleaseMacroHash(&mVarianceHash);
    ReleaseMacroHash(&mTransitionHash);
    ReleaseMacroHash(&mXformHash);
    ReleaseMacroHash(&mXformInstanceHash);
  
    for (i = 0; i < mNumberOfXformsToUpdate; i++) 
    {
      if (NULL != mpXformToUpdate[i].mpShellCommand)
        free(mpXformToUpdate[i].mpShellCommand);
    }
  
    if (NULL != mpXformToUpdate)
      free(mpXformToUpdate);
        
    if (NULL != mpClusterWeightVectors) delete [] mpClusterWeightVectors;
    if (NULL != mpGw)                   delete [] mpGw;
    if (NULL != mpKw)                   delete [] mpKw;
      
  } // Release();
    
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  ModelSet::
  Scan(int mask, HMMSetNodeName nodeName,
       ScanAction action, void *pUserData)
  {
    Macro *macro;
    
    if (nodeName != NULL) 
      strcpy(nodeName+sizeof(HMMSetNodeName)-4, "...");
  
    // walk through the list of macros and decide what to do... 
    // we also decide which direction to walk based on the MTM_REVERSE_PASS flag
    for (macro = mask & MTM_REVERSE_PASS ? mpLastMacro   : mpFirstMacro;
        macro != NULL;
        macro = mask & MTM_REVERSE_PASS ? macro->prevAll : macro->nextAll) 
    {
      assert(macro->mpData != NULL);
      if (macro->mpData->mpMacro != macro) 
        continue;
      
      if (nodeName != NULL) 
        strncpy(nodeName, macro->mpName, sizeof(HMMSetNodeName)-4);
  
      switch (macro->mType) 
      {
        case mt_Xform:
          if (!(mask & MTM_XFORM)) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_XformInstance:
          if (!(mask & (MTM_XFORM | MTM_XFORM_INSTANCE))) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_mean:
          if (!(mask & MTM_MEAN)) break;
          action(mt_mean, nodeName, macro->mpData, pUserData);
          break;
        case mt_variance:
          if (!(mask & MTM_VARIANCE)) break;
          action(mt_variance, nodeName, macro->mpData, pUserData);
          break;
        case mt_transition:
          if (!(mask & MTM_TRANSITION)) break;
          action(mt_transition, nodeName, macro->mpData, pUserData);
          break;
        case mt_mixture:
          if (!(mask & (MTM_ALL & ~(MTM_STATE | MTM_HMM | MTM_TRANSITION)))) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_state:
          if (!(mask & (MTM_ALL & ~(MTM_HMM | MTM_TRANSITION)))) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        case mt_hmm:
          if (!(mask & MTM_ALL)) break;
          macro->mpData->Scan(mask, nodeName, action, pUserData);
          break;
        default: assert(0);
      }
    }
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  NormalizeAccums()
  {
    Scan(MTM_ALL & ~(MTM_XFORM_INSTANCE|MTM_XFORM), NULL,
               NormalizeAccum, NULL);
  }; // NormalizeAccums
  
    
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  ResetAccums()
  {
    if (mAllocAccums)
    {
      Scan(MTM_STATE | MTM_MEAN | MTM_VARIANCE | MTM_TRANSITION,
              NULL, ResetAccum, this);
    }
  }; // ResetAccums()
  
  //***************************************************************************
  //***************************************************************************
  void 
  ModelSet::
  ExpandPredefXforms()
  {
    Scan(MTM_XFORM, NULL, Xform::ExpandPredef, NULL);
  }  
  
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  DistributeMacroOccurances()
  {
    Macro *   macro;
    size_t    i;
    size_t    j;
    size_t    k;
  
    for (k = 0; k < mHmmHash.mNEntries; k++) 
    {
      macro = static_cast <Macro *> (mHmmHash.mpEntry[k]->data);
      
      if (macro->mpData->mpMacro != macro) 
      { 
        macro->mpData->mpMacro->mOccurances += macro->mOccurances;
      }
      
      Hmm *hmm = static_cast <Hmm *> (macro->mpData);
  
      for (i = 0; i < hmm->mNStates - 2; i++) 
      {
        State *state = hmm->mpState[i];
  
        if (state->mpMacro) 
          state->mpMacro->mOccurances += macro->mOccurances;
  
        if (state->mOutPdfKind == KID_DiagC) 
        {
          for (j = 0; j < state->mNMixtures; j++) 
          {
            Mixture *mixture = state->mpMixture[j].mpEstimates;
  
            if (mixture->mpMacro)
              mixture->mpMacro->mOccurances += macro->mOccurances;
            
            if (mixture->mpMean->mpMacro)
              mixture->mpMean->mpMacro->mOccurances += macro->mOccurances;
            
            if (mixture->mpVariance->mpMacro) 
              mixture->mpVariance->mpMacro->mOccurances += macro->mOccurances;
          }
        }
      }
  
      if (hmm->mpTransition->mpMacro) 
        hmm->mpTransition->mpMacro->mOccurances += macro->mOccurances;
    }    
  }; // DistributeMacroOccurances()
  
  
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  ComputeGlobalStats(FLOAT *observation, int time)
  {
    GlobalStatsUserData ud = {observation, time};
    Scan(MTM_STATE | MTM_MIXTURE, NULL, GlobalStats, &ud);
  }; // ComputeGlobalStats(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  UpdateFromAccums()
  {
    Macro * macro;
    size_t  i;
    
    for (i = 0; i < mMeanHash.mNEntries; i++) 
    {
      macro = (Macro *) mMeanHash.mpEntry[i]->data;
  //  for (macro = mean_list; macro != NULL; macro = macro->mpNext) {
      if (macro->mpData->mpMacro != macro) continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Mean vector", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Mean *) macro->mpData)->UpdateFromAccums(this);
      }
    }
  
    for (i = 0; i < mVarianceHash.mNEntries; i++) 
    {
      macro = (Macro *) mVarianceHash.mpEntry[i]->data;
  //  for (macro = variance_list; macro != NULL; macro = macro->mpNext) {
      if (macro->mpData->mpMacro != macro) continue;
      
      if (strncmp(macro->mpName, "varFloor", 8)) 
      {
        if (macro->mOccurances < mMinOccurances) {
          WARN_FEW_EXAMPLES("Variance vector", macro->mpName, macro->mOccurances);
        } 
        else 
        {
          ((Variance *) macro->mpData)->UpdateFromAccums(this);
        }
      }
    }
  
    for (i = 0; i < mTransitionHash.mNEntries; i++) 
    {
      macro = (Macro *) mTransitionHash.mpEntry[i]->data;
  //  for (macro = transition_list; macro != NULL; macro = macro->mpNext) {
      if (macro->mpData->mpMacro != macro) 
        continue;
        
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Transition matrix ", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Transition *) macro->mpData)->UpdateFromAccums(this);
      }
    }               
      
    for (i = 0; i < mMixtureHash.mNEntries; i++) 
    {
      macro = (Macro *) mMixtureHash.mpEntry[i]->data;
  //  for (macro = mixture_list; macro != NULL; macro = macro->mpNext) {
      if (macro->mpData->mpMacro != macro)
        continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Mixture", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Mixture *) macro->mpData)->UpdateFromAccums(this);
      }
    }
  
    for (i = 0; i < mStateHash.mNEntries; i++) 
    {
      macro = (Macro *) mStateHash.mpEntry[i]->data;
  //  for (macro = hmm_set->state_list; macro != NULL; macro = macro->mpNext) {
      if (macro->mpData->mpMacro != macro) 
        continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("State", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((State *) macro->mpData)->UpdateFromAccums(this, NULL);
      }
    }
  
    for (i = 0; i < mHmmHash.mNEntries; i++) 
    {
      macro = (Macro *) mHmmHash.mpEntry[i]->data;
  //  for (macro = hmm_set->hmm_list; macro != NULL; macro = macro->mpNext) {
      if (macro->mpData->mpMacro != macro) 
        continue;
      
      if (macro->mOccurances < mMinOccurances) 
      {
        WARN_FEW_EXAMPLES("Model", macro->mpName, macro->mOccurances);
      } 
      else 
      {
        ((Hmm *) macro->mpData)->UpdateFromAccums(this);
      }
    }    
  }; // UpdateHMMSetFromAccums(...)
  
  
  //**************************************************************************  
  //**************************************************************************  
  Macro *
  ModelSet::
  pAddMacro(const char type, const std::string & rNewName)
  {
    Macro *                  macro;
    MyHSearchData * hash;
    ENTRY                    e;
    ENTRY *                  ep;
  
    switch (type)
    {
      case 'h': hash = &mHmmHash; break;
      case 's': hash = &mStateHash; break;
      case 'm': hash = &mMixtureHash; break;
      case 'u': hash = &mMeanHash; break;
      case 'v': hash = &mVarianceHash; break;
      case 't': hash = &mTransitionHash; break;
      case 'j': hash = &mXformInstanceHash; break;
      case 'x': hash = &mXformHash; break;
      default:  hash = NULL; break;
    }
    
    if (hash == NULL) 
    {
      return NULL;
    }
  
    if ((macro = FindMacro(hash, rNewName.c_str())) != NULL)
    {
      return macro;
    }
    
    // create new macro
    macro = new Macro;
    
    //if ((macro = (Macro *) malloc(sizeof(Macro))) == NULL ||
    if ((macro->mpName = strdup(rNewName.c_str())) == NULL              
    ||  (macro->mpFileName = NULL, mpCurrentMmfName
    &&  (macro->mpFileName = strdup(mpCurrentMmfName)) == NULL)) 
    {
      Error("Insufficient memory");
    }    
    
    e.key  = macro->mpName;
    e.data = macro;
  
    if (!my_hsearch_r(e, ENTER, &ep, hash)) 
    {
      Error("Insufficient memory");
    }
  
    macro->mpData = NULL;
    macro->mOccurances = 0;
    macro->mType = type;
    // List of all macros is made to be able to save macros in proper order
    macro->nextAll = NULL;
    macro->prevAll = mpLastMacro;
    
    if (!mpFirstMacro /* => !hmm_set->mpLastMacro  */) 
    {
      mpFirstMacro = macro;
    } 
    else 
    {
      mpLastMacro->nextAll = macro;
    }
    
    mpLastMacro = macro;
    
    return macro;
  }; // pAddMacro(...)


  //**************************************************************************  
  //**************************************************************************  
  void 
  ModelSet::
  AllocateAccumulatorsForXformStats() 
  {
    mAllMixuresUpdatableFromStatAccums = true;
    
    if(mCmllrStats) {
      for (size_t i=0; i < mNumberOfXformsToUpdate; i++) {
        if (NULL == mpXformToUpdate[i].mpXform->mpCmllrStats) {
          int size = mpXformToUpdate[i].mpXform->mInSize;
          mpXformToUpdate[i].mpXform->mpCmllrStats = (FLOAT *) calloc(sizeof(FLOAT), ((size+size*(size+1)/2) * (size-1) + 1));
        }
      }
    }
    
    Scan(MTM_XFORM_INSTANCE | MTM_MIXTURE, NULL,
         AllocateXformStatCachesAndAccums, this);
  }


  //***************************************************************************
  //**************************************************************************  
  void 
  ModelSet::
  WriteXformStatsAndRunCommands(const char * pOutDir, bool binary)
  {
    HMMSetNodeName                nodeNameBuffer;
    WriteStatsForXformUserData    userData;
    char                          fileName[1024];
    size_t                        i;
    size_t                        j;
    size_t                        k;
    FILE *                        fp;
  
    struct XfStatsHeader 
    {
      #define PRECISION_FLOAT    0
      #define PRECISION_DOUBLE   1
      #define PRECISION_UNKNOWN 15
      char precision;
      #define STATS_MEAN         0
      #define STATS_COV_LOW_TRI  1
      char stats_type;
      short size;
    } header = {sizeof(FLOAT) == sizeof(float)  ? PRECISION_FLOAT  :
                sizeof(FLOAT) == sizeof(double) ? PRECISION_DOUBLE :
                                                  PRECISION_UNKNOWN};
    userData.mBinary = binary;
    for (k = 0; k < mNumberOfXformsToUpdate; k++) 
    {
      char *    ext;
      char *    shellCommand = mpXformToUpdate[k].mpShellCommand;
      
      userData.mpXform = (LinearXform *) mpXformToUpdate[k].mpXform;
      assert(userData.mpXform->mXformType == XT_LINEAR);
  
      MakeFileName(fileName, userData.mpXform->mpMacro->mpName, pOutDir, NULL);
      ext = fileName + strlen(fileName);
      
      strcpy(ext, ".xms"); userData.mMeanFile.mpStatsN = strdup(fileName);
      strcpy(ext, ".xmo"); userData.mMeanFile.mpOccupN = strdup(fileName);
      strcpy(ext, ".xcs"); userData.mCovFile.mpStatsN = strdup(fileName);
      strcpy(ext, ".xco"); userData.mCovFile.mpOccupN = strdup(fileName);
  
      if (userData.mMeanFile.mpStatsN == NULL || userData.mMeanFile.mpOccupN == NULL||
        userData.mCovFile.mpStatsN  == NULL || userData.mCovFile.mpOccupN  == NULL) 
      {
        Error("Insufficient memory");
      }
  
      userData.mMeanFile.mpStatsP = fopen(userData.mMeanFile.mpStatsN, binary?"w":"wt");
      if (userData.mMeanFile.mpStatsP == NULL) 
      {
        Error("Cannot open output file %s",
              userData.mMeanFile.mpStatsN);
      }
  
      if (binary) 
      {
        header.stats_type = STATS_MEAN;
        header.size = userData.mpXform->mInSize;
        
        if (!isBigEndian()) 
          swap2(header.size);
        
        if (fwrite(&header, sizeof(header), 1, userData.mMeanFile.mpStatsP) != 1) 
        {
          Error("Cannot write to file: %s", userData.mMeanFile.mpStatsN);
        }
      }
  
      userData.mMeanFile.mpOccupP = fopen(userData.mMeanFile.mpOccupN, "wt");
      if (userData.mMeanFile.mpOccupP == NULL) 
      {
        Error("Cannot open output file %s", userData.mMeanFile.mpOccupN);
      }
  
      userData.mCovFile.mpStatsP = fopen(userData.mCovFile.mpStatsN, binary?"w":"wt");
      if (userData.mCovFile.mpStatsP == NULL) 
      {
        Error("Cannot open output file %s", userData.mCovFile.mpStatsN);
      }
  
      if (binary) 
      {
        header.stats_type = STATS_COV_LOW_TRI;
        header.size = userData.mpXform->mInSize;
        if (!isBigEndian()) 
          swap2(header.size);
        
        if (fwrite(&header, sizeof(header), 1, userData.mCovFile.mpStatsP) != 1) 
        {
          Error("Cannot write to file: %s",
                userData.mCovFile.mpStatsN);
        }
      }
  
      userData.mCovFile.mpOccupP = fopen(userData.mCovFile.mpOccupN, "wt");
      if (userData.mCovFile.mpOccupP == NULL) 
      {
        Error("Cannot open output file %s", userData.mCovFile.mpOccupN);
      }
  
      if (!mCmllrStats) {
        Scan(MTM_MEAN | MTM_VARIANCE, nodeNameBuffer, WriteStatsForXform, &userData);
      } else {
        size_t dim;
        size_t out_size= userData.mpXform->mOutSize;
        size_t k_size  = userData.mpXform->mInSize;
        size_t G_size  = k_size*(k_size+1)/2;
        size_t stat_size = k_size+G_size;
        assert(k_size - 1 == out_size);
        
        FLOAT norm = userData.mpXform->mpCmllrStats[stat_size * out_size];
        
        if (fprintf(userData.mMeanFile.mpOccupP, "%s ", userData.mpXform->mpMacro->mpName) < 0 ||
            WriteNumber(userData.mMeanFile.mpOccupP, norm) < 0 ||
            fputs("\n", userData.mMeanFile.mpOccupP) == EOF) {
          Error("Cannot write to file: %s", userData.mMeanFile.mpOccupN);
        }
        if (fprintf(userData.mMeanFile.mpOccupP, "%s ", userData.mpXform->mpMacro->mpName) < 0 ||
            WriteNumber(userData.mMeanFile.mpOccupP, norm) < 0 ||
            fputs("\n", userData.mMeanFile.mpOccupP) == EOF) {
          Error("Cannot write to file: %s", userData.mCovFile.mpOccupN);
        }

        for(dim = 0; dim < out_size; dim++) {
          int cc = 0;
          FLOAT *out_vec = &userData.mpXform->mpCmllrStats[dim * stat_size];

          if (binary) {
            size_t i;
            if (!isBigEndian()) {
              for (i = 0; i < k_size; i++) swapFLOAT(out_vec[i]);
            }          
            cc |= (fwrite(out_vec, sizeof(FLOAT), k_size, userData.mMeanFile.mpStatsP) != k_size);
        
            if (!isBigEndian()) 
              for (i = 0; i < k_size; i++) swapFLOAT(out_vec[i]);
          } else {
            size_t j;
            for (j=0;j<k_size;j++) {
              cc |= WriteNumber(userData.mMeanFile.mpStatsP, out_vec[j]) < 0;
              cc |= fputc(' ', userData.mMeanFile.mpStatsP) == EOF;
            }
            cc |= fputs("\n", userData.mMeanFile.mpStatsP) == EOF;
          }

          if (cc) {
            Error("Cannot write to file %s", userData.mMeanFile.mpStatsN);
          }      

          out_vec = &userData.mpXform->mpCmllrStats[dim * stat_size + k_size];
          if (binary) {
            size_t i;
            if (!isBigEndian()) for (i = 0; i < G_size; i++) swapFLOAT(out_vec[i]);
            cc |= fwrite(out_vec, sizeof(FLOAT), G_size, userData.mCovFile.mpStatsP) != G_size;
            if (!isBigEndian()) for (i = 0; i < G_size; i++) swapFLOAT(out_vec[i]);
          } else {
            size_t k, j;
            for (k=0; k < k_size; k++) {
              for (j=0;j<=k;j++) {
                cc |= WriteNumber(userData.mCovFile.mpStatsP, out_vec[k*(k+1)/2+j]) < 0;
                cc |= fputc(' ', userData.mCovFile.mpStatsP) == EOF;
              }
              for (;j<k_size; j++) {
                cc |= WriteNumber(userData.mCovFile.mpStatsP, out_vec[j*(j+1)/2+k]) < 0;
                cc |= fputc(' ', userData.mCovFile.mpStatsP) == EOF;
              }
              cc |= fputs("\n", userData.mCovFile.mpStatsP) == EOF;
            }
            cc |= fputs("\n", userData.mCovFile.mpStatsP) == EOF;
          }
  
          if (cc) {
            Error("Cannot write to file %s", userData.mCovFile.mpStatsN);
          }      
        }
      }

      fclose(userData.mMeanFile.mpStatsP); fclose(userData.mMeanFile.mpOccupP);
      fclose(userData.mCovFile.mpStatsP);  fclose(userData.mCovFile.mpOccupP);
      free(userData.mMeanFile.mpStatsN);   free(userData.mMeanFile.mpOccupN);
      free(userData.mCovFile.mpStatsN);    free(userData.mCovFile.mpOccupN);
  
      strcpy(ext, ".xfm");
      if ((fp = fopen(fileName, "wt")) == NULL) 
      {
        Error("Cannot open output file %s", fileName);
      }
  
      for (i=0; i < userData.mpXform->mOutSize; i++) 
      {
        for (j=0; j < userData.mpXform->mInSize; j++) 
        {
          fputc(' ', fp);
          WriteNumber(fp, userData.mpXform->mMatrix[i][j]);
                       //userData.mpXform->mpMatrixO[i*userData.mpXform->mInSize+j]);
        }
        fputs("\n", fp);
      }
  
      fclose(fp);
  
      if (shellCommand != NULL && *shellCommand) 
      {
        TraceLog("Executing command: %s", shellCommand);
        system(shellCommand);
      }
    }
  } // WriteXformStatsAndRunCommands(const std::string & rOutDir, bool binary)

  
  //****************************************************************************
  //****************************************************************************
  void 
  ModelSet::
  ReadXformStats(const char * pOutDir, bool binary) 
  {
    HMMSetNodeName                nodeNameBuffer;
    WriteStatsForXformUserData   userData;
    char                          fileName[1024];
    size_t                        k;
    FILE *                        fp;
  
    struct XfStatsHeader 
    {
      char    precision;
      char    stats_type;
      short   size;
    } header;
  
    userData.mBinary = binary;
    
    for (k = 0; k < mNumberOfXformsToUpdate; k++) 
    {
      char *  ext;
      userData.mpXform = (LinearXform *) mpXformToUpdate[k].mpXform;
      assert(userData.mpXform->mXformType == XT_LINEAR);
  
      MakeFileName(fileName, userData.mpXform->mpMacro->mpName, pOutDir, NULL);
      ext = fileName + strlen(fileName);
      
      if(!mCmllrStats) {
        strcpy(ext, ".xms"); userData.mMeanFile.mpStatsN = strdup(fileName);
        strcpy(ext, ".xmo"); userData.mMeanFile.mpOccupN = strdup(fileName);
        strcpy(ext, ".xcs"); userData.mCovFile.mpStatsN = strdup(fileName);
        strcpy(ext, ".xco"); userData.mCovFile.mpOccupN = strdup(fileName);
    
        if (userData.mMeanFile.mpStatsN == NULL || userData.mMeanFile.mpOccupN == NULL||
          userData.mCovFile.mpStatsN  == NULL || userData.mCovFile.mpOccupN  == NULL) 
        {
          Error("Insufficient memory");
        }
  
        userData.mMeanFile.mpStatsP = fopen(userData.mMeanFile.mpStatsN, "r");
        
        if (userData.mMeanFile.mpStatsP == NULL) 
        {
          Error("Cannot open input file '%s'", userData.mMeanFile.mpStatsN);
        }
    
        if (binary) 
        {
          header.precision = -1;
          fread(&header, sizeof(header), 1, userData.mMeanFile.mpStatsP);
          
          if (!isBigEndian()) 
            swap2(header.size);
          
          if (ferror(userData.mMeanFile.mpStatsP))
          {
            Error("Cannot read input file '%s'", userData.mMeanFile.mpStatsN);
          }        
          else if (header.stats_type                != STATS_MEAN
               || static_cast<size_t>(header.size)  != userData.mpXform->mInSize
               || header.precision                  != (sizeof(FLOAT) == sizeof(FLOAT_32)
                                                       ? PRECISION_FLOAT : PRECISION_DOUBLE))
          {
            Error("Invalid header in file '%s'", userData.mMeanFile.mpStatsN);
          }
        }
    
        userData.mMeanFile.mpOccupP = fopen(userData.mMeanFile.mpOccupN, "r");
        
        if (userData.mMeanFile.mpOccupP == NULL) 
        {
          Error("Cannot open input file '%s'", userData.mMeanFile.mpOccupN);
        }
    
        userData.mCovFile.mpStatsP = fopen(userData.mCovFile.mpStatsN, "r");
        if (userData.mCovFile.mpStatsP == NULL) 
        {
          Error("Cannot open input file '%s'", userData.mCovFile.mpStatsN);
        }
    
        if (binary) 
        {
          header.precision = -1;
          fread(&header, sizeof(header), 1, userData.mCovFile.mpStatsP);
          if (!isBigEndian()) swap2(header.size);
          if (ferror(userData.mCovFile.mpStatsP)) {
            Error("Cannot read input file '%s'", userData.mCovFile.mpStatsN);
          } else if (header.stats_type               != STATS_COV_LOW_TRI
                 || static_cast<size_t>(header.size) != userData.mpXform->mInSize
                 || header.precision                 != (sizeof(FLOAT) == sizeof(float)
                                                       ? PRECISION_FLOAT : PRECISION_DOUBLE))
          {
            Error("Invalid header in file '%s'", userData.mCovFile.mpStatsN);
          }
        }
    
        userData.mCovFile.mpOccupP = fopen(userData.mCovFile.mpOccupN, "r");
        if (userData.mCovFile.mpOccupP == NULL) {
          Error("Cannot open output file '%s'", userData.mCovFile.mpOccupN);
        }
    
        Scan(MTM_MEAN | MTM_VARIANCE, nodeNameBuffer,
             ReadStatsForXform, &userData);
    
    
        fclose(userData.mMeanFile.mpStatsP); fclose(userData.mMeanFile.mpOccupP);
        fclose(userData.mCovFile.mpStatsP);  fclose(userData.mCovFile.mpOccupP);
        free(userData.mMeanFile.mpStatsN);   free(userData.mMeanFile.mpOccupN);
        free(userData.mCovFile.mpStatsN);    free(userData.mCovFile.mpOccupN);
      }
  
      strcpy(ext, ".xfm");
      if ((fp = fopen(fileName, "r")) == NULL) {
        Error("Cannot open input xforn file '%s'", fileName);
      }
  
      int c =0;
      
      //for (i=0; i < userData.mpXform->mOutSize * userData.mpXform->mInSize; i++) 
      //{
      //  c |= !ReadNumber(fp, &userData.mpXform->mpMatrixO[i]);
      //  fscanf(fp, " ");
      //}
      
      for (size_t i=0; i < userData.mpXform->mOutSize; i++)
      {
        for (size_t j=0; j < userData.mpXform->mInSize; j++) 
        {
          c |= !ReadNumber(fp, &userData.mpXform->mMatrix[i][j]);
        }
      }
      
      if (ferror(fp)) 
      {
        Error("Cannot read xform file '%s'", fileName);
      } else if (c) {
        Error("Invalid xform file '%s'", fileName);
      }
      fclose(fp);
    }
  }


  //*****************************************************************************  
  //*****************************************************************************  
  void 
  ModelSet::
  WriteHMMStats(const char * pFileName)
  {
    Macro *     macro;
    size_t      i = 0;
    size_t      j;
    size_t      k;
    FILE *      fp;
  
    if ((fp = fopen(pFileName, "wt")) == NULL) {
      Error("Cannot open output file: '%s'", pFileName);
    }
  
    for (macro = mpFirstMacro; macro != NULL; macro = macro->nextAll) 
    {
      if (macro->mpData->mpMacro != macro) 
        continue;
        
      if (macro->mType != mt_hmm) 
        continue;
        
      Hmm *hmm = (Hmm *) macro->mpData;
  
      fprintf(fp, "%4d%*c\"%s\" %4ld ", (int) ++i,
                  HIGHER_OF(0,13 - strlen(macro->mpName)), ' ',
                  macro->mpName, macro->mOccurances);
  
      for (j = 0; j < hmm->mNStates-2; j++) 
      {
        State *state = hmm->mpState[j];
        FLOAT stOccP = 0;
        for (k = 0; k < state->mNMixtures; k++) 
        {
          stOccP += state->mpMixture[k].mWeightAccum;
        }
        fputc(' ', fp);
        WriteNumber(fp, stOccP, 10, 6);
      }
      fputs("\n", fp);
    }
  
    fclose(fp);
  }
  
  
  //**************************************************************************  
  //**************************************************************************  
  void 
  ModelSet::
  ResetXformInstances()
  {
    XformInstance* inst;
    for (inst = mpXformInstances; inst != NULL; inst = inst->mpNext) 
    {
      inst->mStatCacheTime = UNDEF_TIME;
      inst->mTime          = UNDEF_TIME;
    }
  }

    
  //**************************************************************************  
  //**************************************************************************  
  void 
  ModelSet::
  UpdateStacks(FLOAT* obs, int time,  PropagDirectionType dir) 
  {
    XformInstance *inst;
    
    for (inst = mpXformInstances; inst != NULL; inst = inst->mpNext) 
    {
      if (inst->mpXform->mDelay > 0) 
      {
        XformPass(inst, obs, time, dir);
      }
    }
  }
      
    
  //**************************************************************************  
  //**************************************************************************  
  MyHSearchData 
  ModelSet::
  MakeCIPhoneHash()
  {
    unsigned int              i;
    int                       nCIphns = 0;
    MyHSearchData             tmpHash;
    MyHSearchData             retHash;
    ENTRY                     e={0}; // {0} is just to make compiler happy
    ENTRY*                    ep;
  
    if (!my_hcreate_r(100, &tmpHash)) 
      Error("Insufficient memory");
  
    // Find CI HMMs and put them into hash
    for (i = 0; i < mHmmHash.mNEntries; i++) 
    {
      Macro* macro = (Macro*) mHmmHash.mpEntry[i]->data;
      if (strpbrk(macro->mpName, "+-")) continue;
  
      e.key  = macro->mpName;
      e.data = macro->mpData;
      
      if (!my_hsearch_r(e, ENTER, &ep, &tmpHash)) 
        Error("Insufficient memory");
      
      nCIphns++;
    }
  
    // Find CD HMMs and mark corresponding CI HMMs in the hash
    for (i = 0; i < mHmmHash.mNEntries; i++) 
    {
      char *    ciname;
      char      chr;
      int       cinlen;
      Macro *   macro = (Macro *) mHmmHash.mpEntry[i]->data;
      
      if (!strpbrk(macro->mpName, "+-")) 
        continue;
  
      ciname = strrchr(macro->mpName, '-');
      
      if (ciname == NULL) ciname = macro->mpName;
      else ciname++;
      
      cinlen = strcspn(ciname, "+");
      chr = ciname[cinlen];
      ciname[cinlen] = '\0';
      e.key  = ciname;
      my_hsearch_r(e, FIND, &ep, &tmpHash);
      ciname[cinlen] = chr;
  
      if (ep != NULL && ep->data != 0) 
      {
        ep->data = NULL;
        nCIphns--;
      }
    }
    
    if (!my_hcreate_r(nCIphns, &retHash)) 
      Error("Insufficient memory");
  
    // To the new hash, collect only those CI HMMs from tmpHash not having CD versions
    for (i = 0; i < tmpHash.mNEntries; i++) 
    {
      if (tmpHash.mpEntry[i]->data != NULL) 
      {
        Hmm *hmm  = (Hmm *) tmpHash.mpEntry[i]->data;
        size_t isTee = hmm->mpTransition->mpMatrixO[hmm->mNStates - 1] > LOG_MIN;
        e.key  = tmpHash.mpEntry[i]->key;
        e.data = reinterpret_cast<void *>(isTee);
        
        if (!my_hsearch_r(e, ENTER, &ep, &retHash)) 
          Error("Insufficient memory");
      }
    }
    
    my_hdestroy_r(&tmpHash, 0);
    return retHash;
  }
  
    
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  ComputeClusterWeightsVector(size_t w)
  {
    // offset to the accumulator matrix row
    size_t start = 0;
    size_t end;
    
    for (size_t i = 0; i < w; i++)
    {
      start += mpClusterWeightVectors[i]->mInSize;
    }
    
    end = start + mpClusterWeightVectors[w]->mInSize;
    
    mpGw[w].Invert();
    mpClusterWeightVectors[w]->mVector.Clear();
    
    for (size_t i = start; i < end; i++)
    {
      for (size_t j=0; j < mpGw[w].Cols(); j++)
      {
        mpClusterWeightVectors[w]->mVector[0][i-start] += mpGw[w][i][j] * mpKw[w][j];
      }
    }
  }

  //**************************************************************************  
  //**************************************************************************  
  void 
  ModelSet::
  ResetClusterWeightVectorsAccums(size_t i)
  {
    mpGw[i].Clear();
    mpKw[i].Clear();
  }

  
  //**************************************************************************  
  //**************************************************************************  
  void
  ModelSet::
  AttachPriors(ModelSet *pPriorModelSet)
  {
    Macro *macro;
    Macro *pPriorMacro;
    HMMSetNodeName nodeName;
    
    strcpy(nodeName+sizeof(HMMSetNodeName)-4, "...");
  
    for (macro = mpFirstMacro; macro != NULL; macro = macro->nextAll) 
    {
      if (macro->mpData->mpMacro != macro) 
        continue;
        
      MyHSearchData * hash;
      hash = macro->mType == 'h' ? &pPriorModelSet->mHmmHash :
             macro->mType == 's' ? &pPriorModelSet->mStateHash :
             macro->mType == 'm' ? &pPriorModelSet->mMixtureHash :
             macro->mType == 'u' ? &pPriorModelSet->mMeanHash :
             macro->mType == 'v' ? &pPriorModelSet->mVarianceHash :
             macro->mType == 't' ? &pPriorModelSet->mTransitionHash : NULL;
             
      if (hash == NULL)
        continue;
        
      if (NULL == (pPriorMacro = FindMacro(hash, macro->mpName))
      ||  macro->mType != pPriorMacro->mType) 
      {
        Error("Macro ~%c \"%s\" is not defined in prior model set", macro->mType, macro->mpName);
      }
      
      if (nodeName != NULL) 
        strncpy(nodeName, macro->mpName, sizeof(HMMSetNodeName)-4);

      switch (macro->mType) 
      {
        case mt_mean:       
          reinterpret_cast<Mean       *>(macro->mpData)->AttachPriors(nodeName, reinterpret_cast<Mean       *>(pPriorMacro->mpData));
          break;
        case mt_variance:   
          reinterpret_cast<Variance   *>(macro->mpData)->AttachPriors(nodeName, reinterpret_cast<Variance   *>(pPriorMacro->mpData));
          break;
        case mt_transition: 
          reinterpret_cast<Transition *>(macro->mpData)->AttachPriors(nodeName, reinterpret_cast<Transition *>(pPriorMacro->mpData));
          break;
        case mt_mixture:    
          reinterpret_cast<Mixture    *>(macro->mpData)->AttachPriors(nodeName, reinterpret_cast<Mixture    *>(pPriorMacro->mpData));
          break;
        case mt_state:      
          reinterpret_cast<State      *>(macro->mpData)->AttachPriors(nodeName, reinterpret_cast<State      *>(pPriorMacro->mpData));
          break;
        case mt_hmm:        
          reinterpret_cast<Hmm        *>(macro->mpData)->AttachPriors(nodeName, reinterpret_cast<Hmm        *>(pPriorMacro->mpData));
          break;
        case mt_Xform:
        case mt_XformInstance: 
          break;
        default: assert(0);
      }   
    }
  }
  
  //*****************************************************************************
  //*****************************************************************************
  void 
  Hmm::
  AttachPriors(HMMSetNodeName nodeName, Hmm * pPriorHmm)
  {
    size_t    i;
    size_t    n = 0;
    char *    chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mNStates != pPriorHmm->mNStates)
    {
      Error("Mismatch in number of states in target and prior HMM's '%s'", nodeName);
    }
      
    for (i=0; i < mNStates-2; i++) 
    {
      if (!mpState[i]->mpMacro) 
      {
        if (n > 0 ) snprintf(chptr, n, ".state[%d]", (int) i+2);
        mpState[i]->AttachPriors(nodeName, pPriorHmm->mpState[i]);
      }
    }
  
    if (!mpTransition->mpMacro) 
    {
      if (n > 0) strncpy(chptr, ".transP", n);
      mpTransition->AttachPriors(nodeName, pPriorHmm->mpTransition);
    }
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  State::
  AttachPriors(HMMSetNodeName nodeName, State * pPriorState)
  {
    size_t    i;
    size_t    n = 0;
    char *    chptr = NULL;
  
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }
  
    if (mOutPdfKind == KID_PDFObsVec) 
    {
      return;
    }
    
    if (mOutPdfKind != pPriorState->mOutPdfKind)
    {
      Error("Mismatch in OutPdfKind of target and prior states '%s'", nodeName);
    }
    
    if (mNMixtures != pPriorState->mNMixtures)
    {
      Error("Mismatch in number of mixtures in target and prior states '%s'", nodeName);
    }
    
    for (i=0; i < mNMixtures; i++)
    {
      if (!mpMixture[i].mpEstimates->mpMacro)
      {
        if (n > 0 ) snprintf(chptr, n, ".mix[%d]", (int) i+1);
        mpMixture[i].mpEstimates->AttachPriors(nodeName, pPriorState->mpMixture[i].mpEstimates);
      }
    }
    
    mpPrior = pPriorState;
  }
    
  //***************************************************************************
  //***************************************************************************
  void 
  Mixture::
  AttachPriors(HMMSetNodeName nodeName, Mixture * pPriorMixture)
  {
    int n = 0;
    char *chptr = NULL;
    
    if (nodeName != NULL) 
    {
      n = strlen(nodeName);
      chptr = nodeName + n;
      n = sizeof(HMMSetNodeName) - 4 - n;
    }  
  
    if (!mpMean->mpMacro) {
      if (n > 0) strncpy(chptr, ".mean", n);
      mpMean->AttachPriors(nodeName, pPriorMixture->mpMean);
    }
  
    if (!mpVariance->mpMacro) {
      if (n > 0) strncpy(chptr, ".cov", n);
      mpVariance->AttachPriors(nodeName, pPriorMixture->mpVariance);
    }  
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  Mean::
  AttachPriors(HMMSetNodeName nodeName, Mean * pPriorMean)
  {
    mpPrior = pPriorMean;
  }
  
  //***************************************************************************
  //***************************************************************************
  void 
  Variance::
  AttachPriors(HMMSetNodeName nodeName, Variance * pPriorVariance)
  {
    mpPrior = pPriorVariance;
  }

  //***************************************************************************
  //***************************************************************************
  void 
  Transition::
  AttachPriors(HMMSetNodeName nodeName, Transition * pPriorTransition)
  {
    mpPrior = pPriorTransition;
  }  

  //**************************************************************************  
  //**************************************************************************  
  void
  Mean::
  ResetClusterWeightVectorsAccums(size_t i)
  {
    mpOccProbAccums[i] = 0;
    
    for (size_t j = 0; j < mCwvAccum.Cols(); j++)
      mCwvAccum[i][j] = 0;
  }
}; //namespace STK  
