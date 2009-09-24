/***************************************************************************
 *   copyright           : (C) 2004-2005 by Lukas Burget,UPGM,FIT,VUT,Brno *
 *   email               : burget@fit.vutbr.cz                             *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#include <math.h>
#include <assert.h>
#include "SigP.h" 

namespace STK
{
  void GenerateDCTMatrix(Matrix<FLOAT> &rMat, const size_t nBasis, const size_t vctLen, const size_t nRepetitions, bool includeC0)
  {
    // can not generate 
    assert(nBasis < vctLen || (nBasis == vctLen && includeC0));
  
    if(rMat.IsInitialized())
    {
      rMat.Destroy();
    }
    
    size_t rows = nBasis * nRepetitions;    
    rMat.Init(rows, vctLen);
    
    FLOAT norm = static_cast<FLOAT>(sqrt(2.0f / static_cast<FLOAT>(vctLen)));
    FLOAT pi_by_n = static_cast<FLOAT>(M_PI) / static_cast<FLOAT>(vctLen);

    // zeroth base component
    size_t r_idx = 0;
    size_t r;
    size_t c;
    if(includeC0)
    {       
      for(c = 0; c < vctLen; c++)
      {
        rMat[r_idx][c] = norm;
      }
      r_idx++;
    }
    
    // other base components
    size_t b = 1;    
    for(; r_idx < nBasis; r_idx++)
    {
      FLOAT v = pi_by_n * static_cast<FLOAT>(b);
      b++;            
      for(c = 0; c < vctLen; c++)
      {
        rMat[r_idx][c] = norm * static_cast<FLOAT>(cos(v * (static_cast<FLOAT>(c) + 0.5f)));
      }
    }    
    
    // repetition
    size_t i;
    for(i = 1; i < nRepetitions; i++)
    {
      for(r = 0; r < nBasis; r++, r_idx++)
      {
        for(c = 0; c < vctLen; c++)
        {
          rMat[r_idx][c] = rMat[r][c];
        }
      }
    }    
    
  }

  void GenerateConstantMatrix(Matrix<FLOAT> &rMat, const size_t nRows, const size_t nCols, const FLOAT constant)
  {
    if(rMat.IsInitialized())
    {
      rMat.Destroy();
    }
    
    rMat.Init(nRows, nCols);

    size_t i;
    size_t j;

    for(i = 0; i < nRows; i++)
    {
      for(j = 0; j < nCols; j++)
      {
        rMat[i][j] = constant;
      }
    }
  }

  void GenerateDiagMatrix(Matrix<FLOAT> &rMat, const size_t nRows, const size_t nCols, const FLOAT val)
  {
    if(rMat.IsInitialized())
    {
      rMat.Destroy();
    }
    
    rMat.Init(nRows, nCols);

    size_t i;
    size_t j;

    for(i = 0; i < nRows; i++)
    {
      for(j = 0; j < nCols; j++)
      {
        if(i != j)
        {
          rMat[i][j] = static_cast<FLOAT>(0);
        }
        else
        {
          rMat[i][j] = val;
        }
      }
    }  
  }
  
  void GenerateRandomMatrix(Matrix<FLOAT> &rMat, const size_t nRows, const size_t nCols, const FLOAT minVal, const FLOAT maxVal, const unsigned int seed)
  {
    if(rMat.IsInitialized())
    {
      rMat.Destroy();
    }
    
    rMat.Init(nRows, nCols);

    srand(seed);
    FLOAT range = maxVal - minVal; 

    size_t i;
    size_t j;

    for(i = 0; i < nRows; i++)
    {
      for(j = 0; j < nCols; j++)
      {         
        rMat[i][j] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * range + minVal;
      }
    }
  }

  void GenerateHammingWindow(BasicVector<FLOAT> &rVct, const size_t hammingLen, const size_t nRepetitions)
  {
    if(rVct.IsInitialized())
    {
      rVct.Destroy();
    }

    size_t vct_len = hammingLen * nRepetitions;      
    rVct.Init(vct_len);
    
    size_t j;
    for(j = 0; j < hammingLen; j++)
    {
      rVct[j] = 0.54f - 0.46f * cos(2.0f * (float)M_PI * j / (hammingLen - 1));
    }
    
    size_t i;
    for(i = 1; i < nRepetitions; i++)
    {
      for(j = 0; j < hammingLen; j++)
      {
        rVct[i * hammingLen + j] = rVct[j]; 
      }
    }    
  }  

  void GenerateRandomWindow(BasicVector<FLOAT> &rVct, const size_t vctLen, const FLOAT minVal, const FLOAT maxVal, unsigned int seed)
  {
    if(rVct.IsInitialized())
    {
      rVct.Destroy();
    }
    
    rVct.Init(vctLen);

    srand(seed);
    FLOAT range = maxVal - minVal; 
    
    size_t j;
    for(j = 0; j < vctLen; j++)
    {
      rVct[j] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * range + minVal;
    }
  }    

  void GenerateConstantWindow(BasicVector<FLOAT> &rVct, const size_t vctLen, const FLOAT constant)
  {
    if(rVct.IsInitialized())
    {
      rVct.Destroy();
    }
    
    rVct.Init(vctLen);

    size_t j;
    for(j = 0; j < vctLen; j++)
    {
      rVct[j] = constant;
    }
  }    

  void GenerateLinSpaceWindow(BasicVector<FLOAT> &rVct, const size_t windowLen, const FLOAT startVal, const FLOAT endVal, const size_t nRepetitions)
  {
    if(rVct.IsInitialized())
    {
      rVct.Destroy();
    }

    size_t vct_len = windowLen * nRepetitions;          
    rVct.Init(vct_len);

    FLOAT k = (endVal - startVal) / static_cast<FLOAT>(windowLen - 1);

    size_t j;
    for(j = 0; j < windowLen; j++)
    {
      rVct[j] = k * static_cast<FLOAT>(j) + startVal;
    }

    size_t i;
    for(i = 1; i < nRepetitions; i++)
    {
      for(j = 0; j < windowLen; j++)
      {
        rVct[i * windowLen + j] = rVct[j]; 
      }
    }    
  }    

  void GenerateTriangWindow(BasicVector<FLOAT> &rVct, const size_t windowLen, const size_t nRepetitions)
  {
    if(rVct.IsInitialized())
    {
      rVct.Destroy();
    }

    size_t vct_len = windowLen * nRepetitions;          
    rVct.Init(vct_len);

    size_t j;
    size_t half  = windowLen / 2;
    if (windowLen > 1)
    {
      for (j = 0; j < half; j++)
        rVct[j] = static_cast<FLOAT>(1.0) / static_cast<FLOAT>(half) * static_cast<FLOAT>(j);
      for (j = half; j < windowLen; j++)
        rVct[j] = (static_cast<FLOAT>(-1.0) / static_cast<FLOAT>(half) * static_cast<FLOAT>(j - half) + static_cast<FLOAT>(1.0));
    }
    else
    {
      rVct[0] = static_cast<FLOAT>(1.0);
    }

    size_t i;
    for(i = 1; i < nRepetitions; i++)
    {
      for(j = 0; j < windowLen; j++)
      {
        rVct[i * windowLen + j] = rVct[j]; 
      }
    }    
  }    
}; // namespace STK
