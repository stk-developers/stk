
#include "Features.h"

namespace STK
{
  //###########################################################################
  //###########################################################################
  // FeatureRepository section
  //###########################################################################
  //###########################################################################

  //***************************************************************************
  //***************************************************************************
  void
  FeatureRepository::
  AddFile(const std::string & rFileName)
  {
    mInputQueue.push_back(rFileName);
  } // AddFile(const std::string & rFileName)

  //***************************************************************************
  //***************************************************************************
  void
  FeatureRepository::
  AddFileList(const std::string & rFileName, const std::string & rFilter)
  {
    IStkStream    l_stream;
    std::string   line;
    
    l_stream.open(rFileName, std::ios::in, rFilter);
    
    if (l_stream.good())
    {
      // read all lines and parse them
      while (!l_stream.eof())
      {
        getline(l_stream, line);
        mInputQueue.push_back(line);
      }
      // close the list file
      l_stream.close();
    }
    else
    {
      //:TODO:
      // Warning or error
    }
  } // AddFileList(const std::string & rFileName)
    
  //***************************************************************************
  //***************************************************************************
  int 
  ReadHTKHeader (FILE * pInFp, HtkHeader * pHeader, bool swap)
  {
    if (!fread(&pHeader->mNSamples,     sizeof(INT_32),  1, pInFp)) return -1;
    if (!fread(&pHeader->mSamplePeriod, sizeof(INT_32),  1, pInFp)) return -1;
    if (!fread(&pHeader->mSampleSize,   sizeof(INT_16),  1, pInFp)) return -1;
    if (!fread(&pHeader->mSampleKind,   sizeof(UINT_16), 1, pInFp)) return -1;
  
    if (swap) 
    {
      swap4(pHeader->mNSamples);
      swap4(pHeader->mSamplePeriod);
      swap2(pHeader->mSampleSize);
      swap2(pHeader->mSampleKind);
    }
  
    if (pHeader->mSamplePeriod < 0 || pHeader->mSamplePeriod > 100000 ||
        pHeader->mNSamples     < 0 || pHeader->mSampleSize   < 0) 
    {
      return -1;
    }
  
    return 0;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  int 
  ReadHTKFeature(
      FILE*     pInFp, 
      FLOAT*    pIn, 
      size_t    feaLen, 
      bool      swap,
      bool      decompress, 
      FLOAT*    pScale, 
      FLOAT*    pBias)
  {
    size_t i;
    
    if (decompress) 
    {
      INT_16 s;
  //    FLOAT pScale = (xmax - xmin) / (2*32767);
  //    FLOAT pBias  = (xmax + xmin) / 2;
  
      for (i = 0; i < feaLen; i++) 
      {
        if (fread(&s, sizeof(INT_16), 1, pInFp) != 1) 
          return -1;
        
        if (swap) swap2(s);
        pIn[i] = (s + pBias[i]) / pScale[i];
      }
      
      return 0;
    }
  
  #if !DOUBLEPRECISION
    if (fread(pIn, sizeof(FLOAT_32), feaLen, pInFp) != feaLen) 
      return -1;
    
    if (swap) 
      for (i = 0; i < feaLen; i++) 
        swap4(pIn[i]);
  #else
    float f;
  
    for (i = 0; i < feaLen; i++) 
    {
      if (fread(&f, sizeof(FLOAT_32), 1, pInFp) != 1)
        return -1;
      
      if (swap) 
        swap4(f);
        
      pIn[i] = f;
    }
  #endif
    return 0;
  }  // int ReadHTKFeature

  
}; // namespace STK
