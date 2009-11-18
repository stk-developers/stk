
#include "Features.h"
#include "Tokenizer.h"

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
  AddFileListToFeatureRepositories(
    const char* pFileName, 
    const char* pFilter, 
    std::queue<FeatureRepository *> &featureRepositoryList)
  {
    IStkStream            l_stream;
    std::string           file_name;
    Tokenizer             file_list(pFileName, ",");
    Tokenizer::iterator   p_file_name;

    //:TODO: error if empty featureRepositoryList
    
    for (p_file_name = file_list.begin(); p_file_name != file_list.end(); ++p_file_name)
    {
      // get rid of initial and trailing blanks
      Trim(*p_file_name);

      // open file name
      l_stream.open(p_file_name->c_str(), std::ios::in, pFilter);
      
      if (!l_stream.good()) {
        //:TODO:
        // Warning or error ... Why warning? -Lukas
        Error("Cannot not open list file %s", p_file_name->c_str());
      }

      // read all lines and parse them
      for(;;)
      {
        l_stream >> file_name;
        //:TODO: if(l_stream.badl()) Error()
        // Reading after last token set the fail bit
        if(l_stream.fail()) 
	  break;
        // we can push_back a std::string as new FileListElem object
        // is created using FileListElem(const std::string&) constructor
        // and logical and physical names are correctly extracted
	featureRepositoryList.front()->mInputQueue.push_back(file_name);
	
	//cycle in the featureRepositoryList
	featureRepositoryList.push(featureRepositoryList.front());
	featureRepositoryList.pop();
      }
      l_stream.close();
    }
  } // AddFileList(const std::string & rFileName)


  //***************************************************************************
  //***************************************************************************
  void
  FeatureRepository::
  Init(
      bool                  swap,
      int                   extLeft,
      int                   extRight,
      int                   targetKind,
      int                   derivOrder,
      int*                  pDerivWinLen,
      const char*           pCmnPath,
      const char*           pCmnMask,
      const char*           pCvnPath,
      const char*           pCvnMask,
      const char*           pCvgFile)
  {
    mSwapFeatures       =   swap;         
    mStartFrameExt      =   extLeft;      
    mEndFrameExt        =   extRight;     
    mTargetKind         =   targetKind;   
    mDerivOrder         =   derivOrder;   
    mDerivWinLengths    =   pDerivWinLen; 
    mpCmnPath           =   pCmnPath;     
    mpCmnMask           =   pCmnMask;     
    mpCvnPath           =   pCvnPath;     
    mpCvnMask           =   pCvnMask;     
    mpCvgFile           =   pCvgFile;     
  } // Init()


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
  AddFileList(const char* pFileName, const char* pFilter)
  {
    IStkStream            l_stream;
    std::string           file_name;
    Tokenizer             file_list(pFileName, ",");
    Tokenizer::iterator   p_file_name;

    for (p_file_name = file_list.begin(); p_file_name != file_list.end(); ++p_file_name)
    {
      // get rid of spaces
      Trim(*p_file_name);

      // open the file
      l_stream.open(p_file_name->c_str(), std::ios::in, pFilter);
      
      if (!l_stream.good())
      {
        //:TODO:
        // Warning or error ... Why warning? -Lukas
        Error("Cannot not open list file %s", p_file_name->c_str());
      }
      // read all lines and parse them
      for(;;)
      {
        l_stream >> file_name;
        //:TODO: if(l_stream.badl()) Error()
        // Reading after last token set the fail bit
        if(l_stream.fail()) 
	  break;
        // we can push_back a std::string as new FileListElem object
        // is created using FileListElem(const std::string&) constructor
        // and logical and physical names are correctly extracted
        mInputQueue.push_back(file_name);
      }
      l_stream.close();
    }
  } // AddFileList(const std::string & rFileName)

  
  //***************************************************************************
  //***************************************************************************
  void
  FeatureRepository::
  MoveNext()
  {
    assert (mInputQueueIterator != mInputQueue.end());
    mInputQueueIterator++;
  } // ReadFullMatrix(Matrix<FLOAT>& rMatrix)


  //***************************************************************************
  //***************************************************************************
  bool
  FeatureRepository::
  ReadFullMatrix(Matrix<FLOAT>& rMatrix)
  {
    // clear the matrix
    rMatrix.Destroy();

    // extract index file name
    if (!mCurrentIndexFileDir.empty())
    {
      char tmp_name[mCurrentIndexFileDir.length() + 
        mCurrentIndexFileExt.length() + 
        mInputQueueIterator->Physical().length()]; 
      
      MakeFileName(tmp_name, mInputQueueIterator->Physical().c_str(), 
          mCurrentIndexFileDir.c_str(), mCurrentIndexFileExt.c_str());
      
      mCurrentIndexFileName = tmp_name;
    }
    else
      mCurrentIndexFileName = "";

    // read the matrix and return the result
    return ReadHTKFeatures(*mInputQueueIterator, rMatrix);
  } // ReadFullMatrix(Matrix<FLOAT>& rMatrix)



  //***************************************************************************
  //***************************************************************************
  // private:
  int 
  FeatureRepository::
  ReadHTKHeader()
  {
    // TODO 
    // Change this... We should read from StkStream
    FILE* fp = mStream.fp();
    
    if (!fread(&mHeader.mNSamples,     sizeof(INT_32),  1, fp)) return -1;
    if (!fread(&mHeader.mSamplePeriod, sizeof(INT_32),  1, fp)) return -1;
    if (!fread(&mHeader.mSampleSize,   sizeof(INT_16),  1, fp)) return -1;
    if (!fread(&mHeader.mSampleKind,   sizeof(UINT_16), 1, fp)) return -1;
  
    if (mSwapFeatures) 
    {
      swap4(mHeader.mNSamples);
      swap4(mHeader.mSamplePeriod);
      swap2(mHeader.mSampleSize);
      swap2(mHeader.mSampleKind);
    }
  
    if (mHeader.mSamplePeriod < 0 
    ||  mHeader.mSamplePeriod > 100000 
    ||  mHeader.mNSamples     < 0 
    ||  mHeader.mSampleSize   < 0) 
    {
      return -1;
    }
  
    return 0;
  }
  
  
  //***************************************************************************
  //***************************************************************************
  // private:
  int 
  FeatureRepository::
  ReadHTKFeature(
      FLOAT*    pIn, 
      size_t    feaLen, 
      bool      decompress, 
      FLOAT*    pScale, 
      FLOAT*    pBias)
  {
    FILE*  fp = mStream.fp();
    
    size_t i;
    
    if (decompress) 
    {
      INT_16 s;
  //    FLOAT pScale = (xmax - xmin) / (2*32767);
  //    FLOAT pBias  = (xmax + xmin) / 2;
  
      for (i = 0; i < feaLen; i++) 
      {
        if (fread(&s, sizeof(INT_16), 1, fp) != 1) 
          return -1;
        
        if (mSwapFeatures) swap2(s);
        pIn[i] = (s + pBias[i]) / pScale[i];
      }
      
      return 0;
    }
  
#if !DOUBLEPRECISION
    if (fread(pIn, sizeof(FLOAT_32), feaLen, fp) != feaLen) 
      return -1;
    
    if (mSwapFeatures) 
      for (i = 0; i < feaLen; i++) 
        swap4(pIn[i]);
#else
    float f;
  
    for (i = 0; i < feaLen; i++) 
    {
      if (fread(&f, sizeof(FLOAT_32), 1, fp) != 1)
        return -1;
      
      if (mSwapFeatures) 
        swap4(f);
        
      pIn[i] = f;
    }
#endif
    return 0;
  }  // int ReadHTKFeature
  
  

  //***************************************************************************
  //***************************************************************************
/*  bool 
  FeatureRepository::
  ReadHTKFeatures(const std::string& rFileName, Matrix<FLOAT>& rFeatureMatrix)
  {
    std::string           file_name(rFileName);
    std::string           cmn_file_name;
    std::string           cvn_file_name;  
    
    int                   ext_left  = mStartFrameExt;
    int                   ext_right = mEndFrameExt;
    int                   from_frame;
    int                   to_frame;
    int                   tot_frames;
    int                   trg_vec_size;
    int                   src_vec_size;
    int                   src_deriv_order;
    int                   lo_src_tgz_deriv_order;
    int                   i;
    int                   j;
    int                   k;
    int                   e;
    int                   coefs;
    int                   trg_E;
    int                   trg_0;
    int                   trg_N;
    int                   src_E;
    int                   src_0;
    int                   src_N;
    int                   comp;
    int                   coef_size;
    char*                 chptr;
  
    
   
    // read frame range definition if any ( physical_file.fea[s,e] )
    if ((chptr = strrchr(file_name.c_str(), '[')) == NULL ||
        ((i=0), sscanf(chptr, "[%d,%d]%n", &from_frame, &to_frame, &i), 
         chptr[i] != '\0')) 
    {
      chptr = NULL;
    }
    
    if (chptr != NULL)                                
      *chptr = '\0';
  
  // Experimental changes...
  // if ((strcmp(file_name.c_str(), "-"))
  // &&  (mpLastFileName != NULL) 
  // &&  (!strcmp(mpLastFileName, file_name.c_str()))) 
  //   {
  //     mHeader = mLastHeader;
  //   } 
  //   else 
  //   {
  //     if (mpLastFileName) 
  //     {
  //       //if (mpFp != stdin) 
  //       //  fclose(mpFp);
  //       mStream.close();
  //       
  //       free(mpLastFileName);
  //       mpLastFileName = NULL;
  //     }

    if ((file_name != "-" )
    &&  (!mLastFileName.empty()) 
    &&  (mLastFileName == file_name)) 
    {
      mHeader = mLastHeader;
    } 
    else 
    {
      if (!mLastFileName.empty()) 
      {
        mStream.close();
        mLastFileName = "";
      }
      
      
      // open the feature file
      mStream.open(file_name.c_str(), ios::binary);
      if (!mStream.good())
      {
        Error("Cannot open feature file: '%s'", file_name.c_str());
      }
      
      
      if (ReadHTKHeader()) 
        Error("Invalid HTK header in feature file: '%s'", file_name.c_str());
      
      if (mHeader.mSampleKind & PARAMKIND_C) 
      {
        // File is in compressed form, scale and pBias vectors
        // are appended after HTK header.
  
        int coefs = mHeader.mSampleSize/sizeof(INT_16);
        mpA = (FLOAT*) realloc(mpA, coefs * sizeof(FLOAT));
        mpB = (FLOAT*) realloc(mpB, coefs * sizeof(FLOAT));
        if (mpA == NULL || mpB == NULL) Error("Insufficient memory");
  
        e  = ReadHTKFeature(mpA, coefs, 0, 0, 0);
        e |= ReadHTKFeature(mpB, coefs, 0, 0, 0);
        
        if (e) 
          Error("Cannot read feature file: '%s'", file_name.c_str());
        
        mHeader.mNSamples -= 2 * sizeof(FLOAT_32) / sizeof(INT_16);
      }
      
      // remember current settings
      mLastFileName = file_name;
      mLastHeader   = mHeader;
    }
    
    if (chptr != NULL) 
      *chptr = '[';
  
    if (chptr == NULL) 
    { // Range [s,e] was not specified
      from_frame = 0;
      to_frame   = mHeader.mNSamples-1;
    }
    
    src_deriv_order = PARAMKIND_T & mHeader.mSampleKind ? 3 :
                      PARAMKIND_A & mHeader.mSampleKind ? 2 :
                      PARAMKIND_D & mHeader.mSampleKind ? 1 : 0;
    src_E =  (PARAMKIND_E & mHeader.mSampleKind) != 0;
    src_0 =  (PARAMKIND_0 & mHeader.mSampleKind) != 0;
    src_N = ((PARAMKIND_N & mHeader.mSampleKind) != 0) * (src_E + src_0);
    comp =    PARAMKIND_C & mHeader.mSampleKind;
    
    mHeader.mSampleKind &= ~PARAMKIND_C;
  
    if (mTargetKind == PARAMKIND_ANON) 
    {
      mTargetKind = mHeader.mSampleKind;
    } 
    else if ((mTargetKind & 077) == PARAMKIND_ANON) 
    {
      mTargetKind &= ~077;
      mTargetKind |= mHeader.mSampleKind & 077;
    }
    
    trg_E = (PARAMKIND_E & mTargetKind) != 0;
    trg_0 = (PARAMKIND_0 & mTargetKind) != 0;
    trg_N =((PARAMKIND_N & mTargetKind) != 0) * (trg_E + trg_0);
  
    coef_size     = comp ? sizeof(INT_16) : sizeof(FLOAT_32);
    coefs         = (mHeader.mSampleSize/coef_size + src_N) / 
                    (src_deriv_order+1) - src_E - src_0;
    src_vec_size  = (coefs + src_E + src_0) * (src_deriv_order+1) - src_N;
  
    //Is coefs dividable by 1 + number of derivatives specified in header
    if (src_vec_size * coef_size != mHeader.mSampleSize) 
    {
      Error("Invalid HTK header in feature file: '%s'. "
            "mSampleSize do not match with parmKind", file_name.c_str());
    }
    
    if (mDerivOrder < 0) 
      mDerivOrder = src_deriv_order;
  
  
    if ((!src_E && trg_E) || (!src_0 && trg_0) || (src_N && !trg_N) ||
        (trg_N && !trg_E && !trg_0) || (trg_N && !mDerivOrder) ||
        (src_N && !src_deriv_order && mDerivOrder) ||
        ((mHeader.mSampleKind & 077) != (mTargetKind & 077) &&
         (mHeader.mSampleKind & 077) != PARAMKIND_ANON)) 
    {
      char srcParmKind[64];
      char trgParmKind[64];
      
      ParmKind2Str(mHeader.mSampleKind, srcParmKind);
      ParmKind2Str(mTargetKind,       trgParmKind);
      Error("Cannot convert %s to %s", srcParmKind, trgParmKind);
    }
  
    lo_src_tgz_deriv_order = LOWER_OF(src_deriv_order, mDerivOrder);
    trg_vec_size  = (coefs + trg_E + trg_0) * (mDerivOrder+1) - trg_N;
    
    i =  LOWER_OF(from_frame, mStartFrameExt);
    from_frame  -= i;
    ext_left     -= i;
  
    i =  LOWER_OF(mHeader.mNSamples-to_frame-1, mEndFrameExt);
    to_frame    += i;
    ext_right    -= i;
  
    if (from_frame > to_frame || from_frame >= mHeader.mNSamples || to_frame< 0)
      Error("Invalid frame range for feature file: '%s'", file_name.c_str());
    
    tot_frames = to_frame - from_frame + 1 + ext_left + ext_right;
    
    // initialize matrix 
    rFeatureMatrix.Init(tot_frames, trg_vec_size);
    
    // fill the matrix with features
    for (i = 0; i <= to_frame - from_frame; i++) 
    {
      FLOAT* A      = mpA;
      FLOAT* B      = mpB;
      FLOAT* mxPtr  = rFeatureMatrix[i+ext_left];
      
      // seek to the desired position
      fseek(mStream.fp(), 
          sizeof(HtkHeader) + (comp ? src_vec_size * 2 * sizeof(FLOAT_32) : 0)
          + (from_frame + i) * src_vec_size * coef_size, 
          SEEK_SET);
  
      e = ReadHTKFeature(mxPtr, coefs, comp, A, B);
      
      mxPtr += coefs; 
      A     += coefs; 
      B     += coefs;
        
      if (src_0 && !src_N) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
      if (trg_0 && !trg_N) mxPtr++;
      if (src_E && !src_N) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
      if (trg_E && !trg_N) mxPtr++;
  
      for (j = 0; j < lo_src_tgz_deriv_order; j++) 
      {
        e |= ReadHTKFeature(mxPtr, coefs, comp, A, B);
        mxPtr += coefs; 
        A     += coefs; 
        B     += coefs;
        
        if (src_0) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
        if (trg_0) mxPtr++;
        if (src_E) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
        if (trg_E) mxPtr++;
      }
      
      if (e) 
        Error("Cannot read feature file: '%s' frame %d/%d", file_name.c_str(),
            i, to_frame - from_frame + 1);
    }
  
    // From now, coefs includes also trg_0 + trg_E !
    coefs += trg_0 + trg_E; 
    
    // If extension of the matrix to the left or to the right is required,
    // perform it here
    for (i = 0; i < ext_left; i++) 
    {
      memcpy(rFeatureMatrix[i],
             rFeatureMatrix[ext_left],
             (coefs * (1+lo_src_tgz_deriv_order) - trg_N) * sizeof(FLOAT));
    }
    
    for (i = tot_frames - ext_right; i < tot_frames; i++) 
    {
      memcpy(rFeatureMatrix[i],
             rFeatureMatrix[tot_frames - ext_right - 1],
             (coefs * (1+lo_src_tgz_deriv_order) - trg_N) * sizeof(FLOAT));
    }

    // Sentence cepstral mean normalization
    if( (mpCmnPath == NULL)
    && !(PARAMKIND_Z & mHeader.mSampleKind) 
    &&  (PARAMKIND_Z & mTargetKind)) 
    {
      // for each coefficient
      for(j=0; j < coefs; j++) 
      {          
        FLOAT norm = 0.0;
        for(i=0; i < tot_frames; i++)      // for each frame
        {
          norm += rFeatureMatrix[i][j - trg_N];
          //norm += fea_mx[i*trg_vec_size - trg_N + j];
        }
        
        norm /= tot_frames;
  
        for(i=0; i < tot_frames; i++)      // for each frame
          rFeatureMatrix[i][j - trg_N] -= norm;
          //fea_mx[i*trg_vec_size - trg_N + j] -= norm;
      }
    }
    
    // Compute missing derivatives
    for (; src_deriv_order < mDerivOrder; src_deriv_order++) 
    { 
      int winLen = mDerivWinLengths[src_deriv_order];
      FLOAT norm = 0.0;
      
      for (k = 1; k <= winLen; k++) 
      {
        norm += 2 * k * k;
      }
      
      // for each frame
      for (i=0; i < tot_frames; i++) 
      {        
        // for each coefficient
        for (j=0; j < coefs; j++) 
        {          
          //FLOAT* src = fea_mx + i*trg_vec_size + src_deriv_order*coefs - trg_N + j;
          FLOAT* src = &rFeatureMatrix[i][src_deriv_order*coefs - trg_N + j];
          
          *(src + coefs) = 0.0;
          
          if (i < winLen || i >= tot_frames-winLen) 
          { // boundaries need special treatment
            for (k = 1; k <= winLen; k++) 
            {  
              *(src+coefs) += k*(src[ LOWER_OF(tot_frames-1-i,k)*rFeatureMatrix.Stride()]
                                -src[-LOWER_OF(i,             k)*rFeatureMatrix.Stride()]);
            }
          } 
          else 
          { // otherwise use more efficient code
            for (k = 1; k <= winLen; k++) 
            {  
              *(src+coefs) += k*(src[ k * rFeatureMatrix.Stride()]
                                -src[-k * rFeatureMatrix.Stride()]);
            }
          }
          *(src + coefs) /= norm;
        }
      }
    }
    
    mHeader.mNSamples    = tot_frames;
    mHeader.mSampleSize  = trg_vec_size * sizeof(FLOAT_32);
    mHeader.mSampleKind  = mTargetKind & ~(PARAMKIND_D | PARAMKIND_A | PARAMKIND_T);
  

    ////////////////////////////////////////////////////////////////////////////
    /////////////// Cepstral mean and variance normalization ///////////////////
    ////////////////////////////////////////////////////////////////////////////
    //.........................................................................
    if (mpCmnPath != NULL
    &&  mpCmnMask != NULL) 
    {
      // retrieve file name
      ProcessMask(file_name, mpCmnMask, cmn_file_name);
      // add the path correctly
      cmn_file_name.insert(0, "/");
      cmn_file_name.insert(0, mpCmnPath);

      // read the file
      ReadCepsNormFile(cmn_file_name.c_str(), &mpLastCmnFile, &mpCmn,
          mHeader.mSampleKind & ~PARAMKIND_Z, CNF_Mean, coefs);
                      
      // recompute feature values
      for (i=0; i < tot_frames; i++) 
      {
        for (j=trg_N; j < coefs; j++) 
        {
          rFeatureMatrix[i][j - trg_N] -= mpCmn[j];
        }
      }
    }
  
    mHeader.mSampleKind |= mDerivOrder==3 ? PARAMKIND_D | PARAMKIND_A | PARAMKIND_T :
                           mDerivOrder==2 ? PARAMKIND_D | PARAMKIND_A :
                           mDerivOrder==1 ? PARAMKIND_D : 0;
  
    //.........................................................................
    if (mpCvnPath != NULL
    &&  mpCvnMask != NULL) 
    {
      // retrieve file name
      ProcessMask(file_name, mpCvnMask, cvn_file_name);
      // add the path correctly
      cvn_file_name.insert(0, "/");
      cvn_file_name.insert(0, mpCvnPath);

      // read the file
      ReadCepsNormFile(cvn_file_name.c_str(), &mpLastCvnFile, &mpCvn,
          mHeader.mSampleKind, CNF_Variance, trg_vec_size);
                      
      // recompute feature values
      for (i=0; i < tot_frames; i++) 
      {
        for (j=trg_N; j < trg_vec_size; j++) 
        {
          rFeatureMatrix[i][j - trg_N] *= mpCvn[j];
        }
      }
    }
    
    //.........................................................................
    // process the global covariance file
    if (mpCvgFile != NULL) 
    {
      ReadCepsNormFile(mpCvgFile, &mpLastCvgFile, &mpCvg,
                      -1, CNF_VarScale, trg_vec_size);
                      
      // recompute feature values
      for (i=0; i < tot_frames; i++) 
      {
        for (j=trg_N; j < trg_vec_size; j++) 
        {
          rFeatureMatrix[i][j - trg_N] *= mpCvg[j];
        }
      }
    }
    
    return true;
  }
*/

  
  //***************************************************************************
  //***************************************************************************
  bool 
  FeatureRepository::
  ReadHTKFeatures(const FileListElem&    rFileNameRecord, 
                        Matrix<FLOAT>&        rFeatureMatrix)
  {
    std::string           file_name(rFileNameRecord.Physical());
    std::string           cmn_file_name;
    std::string           cvn_file_name;  
    
    int                   ext_left  = mStartFrameExt;
    int                   ext_right = mEndFrameExt;
    int                   from_frame;
    int                   to_frame;
    int                   tot_frames;
    int                   trg_vec_size;
    int                   src_vec_size;
    int                   src_deriv_order;
    int                   lo_src_tgz_deriv_order;
    int                   i;
    int                   j;
    int                   k;
    int                   e;
    int                   coefs;
    int                   trg_E;
    int                   trg_0;
    int                   trg_N;
    int                   src_E;
    int                   src_0;
    int                   src_N;
    int                   comp;
    int                   coef_size;
    char*                 chptr;
  
   
    // read frame range definition if any ( physical_file.fea[s,e] )
    if ((chptr =  const_cast<char *> (strrchr(file_name.c_str(), '['))) == NULL ||
        ((i=0), sscanf(chptr, "[%d,%d]%n", &from_frame, &to_frame, &i), 
         chptr[i] != '\0')) 
    {
      chptr = NULL;
    }
    
    if (chptr != NULL)                                
      *chptr = '\0';
  

    if ((file_name != "-" )
    &&  (!mLastFileName.empty()) 
    &&  (mLastFileName == file_name)) 
    {
      mHeader = mLastHeader;
    } 
    else 
    {
      if (!mLastFileName.empty()) 
      {
        mStream.close();
        mLastFileName = "";
      }
      
      
      // open the feature file
      mStream.open(file_name.c_str(), std::ios::binary);
      if (!mStream.good())
      {
        Error("Cannot open feature file: '%s'", file_name.c_str());
      }
      
      
      if (ReadHTKHeader())  {
        Error("Invalid HTK header in feature file: '%s'", file_name.c_str());
      }
      
      if (mHeader.mSampleKind & PARAMKIND_C) 
      {
        // File is in compressed form, scale and pBias vectors
        // are appended after HTK header.
        int coefs = mHeader.mSampleSize/sizeof(INT_16);

        mpA = (FLOAT*) realloc(mpA, coefs * sizeof(FLOAT));
        mpB = (FLOAT*) realloc(mpB, coefs * sizeof(FLOAT));

        if (mpA == NULL || mpB == NULL) {
          Error("Insufficient memory");
        }
  
        e  = ReadHTKFeature(mpA, coefs, 0, 0, 0);
        e |= ReadHTKFeature(mpB, coefs, 0, 0, 0);
        
        if (e) {
          Error("Cannot read feature file: '%s'", file_name.c_str());
        }
        
        mHeader.mNSamples -= 2 * sizeof(FLOAT_32) / sizeof(INT_16);
      }
      
      // remember current settings
      mLastFileName = file_name;
      mLastHeader   = mHeader;
    }
    
    if (chptr != NULL) {
      *chptr = '[';
    }
  
    if (chptr == NULL) { 
      // Range [s,e] was not specified
      from_frame = 0;
      to_frame   = mHeader.mNSamples-1;
    }
    
    src_deriv_order = PARAMKIND_T & mHeader.mSampleKind ? 3 :
                      PARAMKIND_A & mHeader.mSampleKind ? 2 :
                      PARAMKIND_D & mHeader.mSampleKind ? 1 : 0;
    src_E =  (PARAMKIND_E & mHeader.mSampleKind) != 0;
    src_0 =  (PARAMKIND_0 & mHeader.mSampleKind) != 0;
    src_N = ((PARAMKIND_N & mHeader.mSampleKind) != 0) * (src_E + src_0);
    comp =    PARAMKIND_C & mHeader.mSampleKind;
    
    mHeader.mSampleKind &= ~PARAMKIND_C;
  
    if (mTargetKind == PARAMKIND_ANON) 
    {
      mTargetKind = mHeader.mSampleKind;
    } 
    else if ((mTargetKind & 077) == PARAMKIND_ANON) 
    {
      mTargetKind &= ~077;
      mTargetKind |= mHeader.mSampleKind & 077;
    }
    
    trg_E = (PARAMKIND_E & mTargetKind) != 0;
    trg_0 = (PARAMKIND_0 & mTargetKind) != 0;
    trg_N =((PARAMKIND_N & mTargetKind) != 0) * (trg_E + trg_0);
  
    coef_size     = comp ? sizeof(INT_16) : sizeof(FLOAT_32);
    coefs         = (mHeader.mSampleSize/coef_size + src_N) / 
                    (src_deriv_order+1) - src_E - src_0;
    src_vec_size  = (coefs + src_E + src_0) * (src_deriv_order+1) - src_N;
  
    //Is coefs dividable by 1 + number of derivatives specified in header
    if (src_vec_size * coef_size != mHeader.mSampleSize) 
    {
      Error("Invalid HTK header in feature file: '%s'. "
            "mSampleSize do not match with parmKind", file_name.c_str());
    }
    
    if (mDerivOrder < 0) 
      mDerivOrder = src_deriv_order;
  
  
    if ((!src_E && trg_E) || (!src_0 && trg_0) || (src_N && !trg_N) ||
        (trg_N && !trg_E && !trg_0) || (trg_N && !mDerivOrder) ||
        (src_N && !src_deriv_order && mDerivOrder) ||
        ((mHeader.mSampleKind & 077) != (mTargetKind & 077) &&
         (mHeader.mSampleKind & 077) != PARAMKIND_ANON)) 
    {
      char srcParmKind[64];
      char trgParmKind[64];
      
      ParmKind2Str(mHeader.mSampleKind, srcParmKind);
      ParmKind2Str(mTargetKind,       trgParmKind);
      Error("Cannot convert %s to %s", srcParmKind, trgParmKind);
    }
  
    lo_src_tgz_deriv_order = LOWER_OF(src_deriv_order, mDerivOrder);
    trg_vec_size  = (coefs + trg_E + trg_0) * (mDerivOrder+1) - trg_N;
    
    i =  LOWER_OF(from_frame, mStartFrameExt);
    from_frame  -= i;
    ext_left     -= i;
  
    i =  LOWER_OF(mHeader.mNSamples-to_frame-1, mEndFrameExt);
    to_frame    += i;
    ext_right    -= i;
  
    if (from_frame > to_frame || from_frame >= mHeader.mNSamples || to_frame< 0)
      Error("Invalid frame range for feature file: '%s'", file_name.c_str());
    
    tot_frames = to_frame - from_frame + 1 + ext_left + ext_right;
    
    // initialize matrix 
    rFeatureMatrix.Init(tot_frames, trg_vec_size);
    
    // fill the matrix with features
    for (i = 0; i <= to_frame - from_frame; i++) 
    {
      FLOAT* A      = mpA;
      FLOAT* B      = mpB;
      FLOAT* mxPtr  = rFeatureMatrix[i+ext_left];
      
      // seek to the desired position
      fseek(mStream.fp(), 
          sizeof(HtkHeader) + (comp ? src_vec_size * 2 * sizeof(FLOAT_32) : 0)
          + (from_frame + i) * src_vec_size * coef_size, 
          SEEK_SET);
  
      e = ReadHTKFeature(mxPtr, coefs, comp, A, B);
      
      mxPtr += coefs; 
      A     += coefs; 
      B     += coefs;
        
      if (src_0 && !src_N) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
      if (trg_0 && !trg_N) mxPtr++;
      if (src_E && !src_N) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
      if (trg_E && !trg_N) mxPtr++;
  
      for (j = 0; j < lo_src_tgz_deriv_order; j++) 
      {
        e |= ReadHTKFeature(mxPtr, coefs, comp, A, B);
        mxPtr += coefs; 
        A     += coefs; 
        B     += coefs;
        
        if (src_0) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
        if (trg_0) mxPtr++;
        if (src_E) e |= ReadHTKFeature(mxPtr, 1, comp, A++, B++);
        if (trg_E) mxPtr++;
      }
      
      if (e) 
        Error("Cannot read feature file: '%s' frame %d/%d", file_name.c_str(),
            i, to_frame - from_frame + 1);
    }
  
    // From now, coefs includes also trg_0 + trg_E !
    coefs += trg_0 + trg_E; 
    
    // If extension of the matrix to the left or to the right is required,
    // perform it here
    for (i = 0; i < ext_left; i++) 
    {
      memcpy(rFeatureMatrix[i],
             rFeatureMatrix[ext_left],
             (coefs * (1+lo_src_tgz_deriv_order) - trg_N) * sizeof(FLOAT));
    }
    
    for (i = tot_frames - ext_right; i < tot_frames; i++) 
    {
      memcpy(rFeatureMatrix[i],
             rFeatureMatrix[tot_frames - ext_right - 1],
             (coefs * (1+lo_src_tgz_deriv_order) - trg_N) * sizeof(FLOAT));
    }

    // Sentence cepstral mean normalization
    if( (mpCmnPath == NULL)
    && !(PARAMKIND_Z & mHeader.mSampleKind) 
    &&  (PARAMKIND_Z & mTargetKind)) 
    {
      // for each coefficient
      for(j=0; j < coefs; j++) 
      {          
        FLOAT norm = 0.0;
        for(i=0; i < tot_frames; i++)      // for each frame
        {
          norm += rFeatureMatrix[i][j - trg_N];
          //norm += fea_mx[i*trg_vec_size - trg_N + j];
        }
        
        norm /= tot_frames;
  
        for(i=0; i < tot_frames; i++)      // for each frame
          rFeatureMatrix[i][j - trg_N] -= norm;
          //fea_mx[i*trg_vec_size - trg_N + j] -= norm;
      }
    }
    
    // Compute missing derivatives
    for (; src_deriv_order < mDerivOrder; src_deriv_order++) 
    { 
      int winLen = mDerivWinLengths[src_deriv_order];
      FLOAT norm = 0.0;
      
      for (k = 1; k <= winLen; k++) 
      {
        norm += 2 * k * k;
      }
      
      // for each frame
      for (i=0; i < tot_frames; i++) 
      {        
        // for each coefficient
        for (j=0; j < coefs; j++) 
        {          
          //FLOAT* src = fea_mx + i*trg_vec_size + src_deriv_order*coefs - trg_N + j;
          FLOAT* src = &rFeatureMatrix[i][src_deriv_order*coefs - trg_N + j];
          
          *(src + coefs) = 0.0;
          
          if (i < winLen || i >= tot_frames-winLen) 
          { // boundaries need special treatment
            for (k = 1; k <= winLen; k++) 
            {  
              *(src+coefs) += k*(src[ LOWER_OF(tot_frames-1-i,k)*rFeatureMatrix.Stride()]
                                -src[-LOWER_OF(i,             k)*rFeatureMatrix.Stride()]);
            }
          } 
          else 
          { // otherwise use more efficient code
            for (k = 1; k <= winLen; k++) 
            {  
              *(src+coefs) += k*(src[ k * rFeatureMatrix.Stride()]
                                -src[-k * rFeatureMatrix.Stride()]);
            }
          }
          *(src + coefs) /= norm;
        }
      }
    }
    
    mHeader.mNSamples    = tot_frames;
    mHeader.mSampleSize  = trg_vec_size * sizeof(FLOAT_32);
    mHeader.mSampleKind  = mTargetKind & ~(PARAMKIND_D | PARAMKIND_A | PARAMKIND_T);
  

    ////////////////////////////////////////////////////////////////////////////
    /////////////// Cepstral mean and variance normalization ///////////////////
    ////////////////////////////////////////////////////////////////////////////
    //.........................................................................
    if (mpCmnPath != NULL
    &&  mpCmnMask != NULL) 
    {
      // retrieve file name
      ProcessMask(rFileNameRecord.Logical(), mpCmnMask, cmn_file_name);
      // add the path correctly
      cmn_file_name.insert(0, "/");
      cmn_file_name.insert(0, mpCmnPath);

      // read the file
      ReadCepsNormFile(cmn_file_name.c_str(), &mpLastCmnFile, &mpCmn,
          mHeader.mSampleKind & ~PARAMKIND_Z, CNF_Mean, coefs);
                      
      // recompute feature values
      for (i=0; i < tot_frames; i++) 
      {
        for (j=trg_N; j < coefs; j++) 
        {
          rFeatureMatrix[i][j - trg_N] -= mpCmn[j];
        }
      }
    }
  
    mHeader.mSampleKind |= mDerivOrder==3 ? PARAMKIND_D | PARAMKIND_A | PARAMKIND_T :
                           mDerivOrder==2 ? PARAMKIND_D | PARAMKIND_A :
                           mDerivOrder==1 ? PARAMKIND_D : 0;
  
    //.........................................................................
    if (mpCvnPath != NULL
    &&  mpCvnMask != NULL) 
    {
      // retrieve file name
      ProcessMask(rFileNameRecord.Logical(), mpCvnMask, cvn_file_name);
      // add the path correctly
      cvn_file_name.insert(0, "/");
      cvn_file_name.insert(0, mpCvnPath);

      // read the file
      ReadCepsNormFile(cvn_file_name.c_str(), &mpLastCvnFile, &mpCvn,
          mHeader.mSampleKind, CNF_Variance, trg_vec_size);
                      
      // recompute feature values
      for (i=0; i < tot_frames; i++) 
      {
        for (j=trg_N; j < trg_vec_size; j++) 
        {
          rFeatureMatrix[i][j - trg_N] *= mpCvn[j];
        }
      }
    }
    
    //.........................................................................
    // process the global covariance file
    if (mpCvgFile != NULL) 
    {
      ReadCepsNormFile(mpCvgFile, &mpLastCvgFile, &mpCvg,
                      -1, CNF_VarScale, trg_vec_size);
                      
      // recompute feature values
      for (i=0; i < tot_frames; i++) 
      {
        for (j=trg_N; j < trg_vec_size; j++) 
        {
          rFeatureMatrix[i][j - trg_N] *= mpCvg[j];
        }
      }
    }
    
    return true;
  }
}; // namespace STK
