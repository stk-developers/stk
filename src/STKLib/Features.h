//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STK_Features_h
#define STK_Features_h


#define FEATURE_REPOSITORY_PTHREAD_ENABLE 0
#define FEATURE_REPOSITORY_DEFAULT_READ_OVERHEAD

//*****************************************************************************
//*****************************************************************************
// Specific includes
//
#include "common.h"
#include "fileio.h"
#include "Matrix.h"
#include "stkstream.h"


//*****************************************************************************
//*****************************************************************************
// Standard includes
//
#include <list>
#include <queue>
#include <string>


//*****************************************************************************
//*****************************************************************************
// Code ...
//

namespace STK
{
  /** *************************************************************************
   * @brief 
   */
  class FeatureRepository
  {

  //////////////////////////////////////////////////////////////////////////////
  //  PUBLIC SECTION
  //////////////////////////////////////////////////////////////////////////////
  public:
    /// Iterates through the list of feature file records
    typedef   std::list<FileListElem>::iterator  ListIterator;

    // some params for loading features
    bool                        mSwapFeatures;
    int                         mStartFrameExt;
    int                         mEndFrameExt;
    int                         mTargetKind;
    int                         mDerivOrder;
    int*                        mDerivWinLengths;
    const char*                 mpCvgFile;
    //:TODO: get rid of these
    const char*                 mpCmnPath;
    const char*                 mpCmnMask;
    const char*                 mpCvnPath;
    const char*                 mpCvnMask;
    
    
    // Constructors and destructors
    /**
     * @brief Default constructor that creates an empty repository
     */
    FeatureRepository() : mDerivWinLengths(NULL), mpCvgFile(NULL), 
       mpCmnPath(NULL), mpCmnMask(NULL), mpCvnPath(NULL), mpCvnMask(NULL),
       mpLastFileName(NULL), mLastFileName(""), mpLastCmnFile (NULL), 
       mpLastCvnFile (NULL), mpLastCvgFile (NULL), mpCmn(NULL), 
       mpCvn(NULL), mpCvg(NULL), mpA(NULL), mpB(NULL)
    { 
      mInputQueueIterator        = mInputQueue.end(); 

#if FEATURE_REPOSITORY_PTHREAD_ENABLE
      mpRoMutex        = new pthread_mutex_t;
      mpRoCondNotFull  = new pthread_cond_t;
      mpRoCondNotEmpty = new pthread_cond_t;

      pthread_mutex_init(mpRoMutex, NULL);
      pthread_cond_init(mpRoCondNotFull, NULL);
      pthread_cond_init(mpRoCondNotEmpty, NULL);
#endif
    }

    /**
     * @brief Destroys the repository
     */
    ~FeatureRepository()
    {
      if (NULL != mpA) {
        free(mpA);
      }

      if (NULL != mpB) {
        free(mpB);
      }

#if FEATURE_REPOSITORY_PTHREAD_ENABLE
      delete mpRoMutex;
      delete mpRoCondNotFull;
      delete mpRoCondNotEmpty;
#endif
    }


    /**
     * @brief Initializes the object using the given parameters
     *
     * @param swap          Boolean value specifies whether to swap bytes 
     *                      when reading file or not. 
     * @param extLeft       Features read from file are extended with extLeft 
     *                      initial frames. Normally, these frames are 
     *                      repetitions of the first feature frame in file 
     *                      (with its derivative, if derivatives are preset in
     *                      the file). However, if segment of feature frames 
     *                      is extracted according to range specification, the 
     *                      true feature frames from beyond the segment boundary
     *                      are used, wherever it is possible. Note that value 
     *                      of extLeft can be also negative. In such case
     *                      corresponding number of initial frames is discarded. 
     * @param extRight      The paramerer is complementary to parameter extLeft 
     *                      and has obvious meaning. (Controls extensions over
     *                      the last frame, last frame from file is repeated 
     *                      only if necessary).
     * @param targetKind    The parameters is used to check whether 
     *                      pHeader->mSampleKind match to requited targetKind 
     *                      and to control suppression of 0'th cepstral or 
     *                      energy coefficients accorging to modifiers _E, _0, 
     *                      and _N. Modifiers _D, _A and _T are ignored; 
     *                      Computation of derivatives is controled by parameters
     *                      derivOrder and derivWinLen. Value PARAMKIND_ANON 
     *                      ensures that function do not result in targetKind 
     *                      mismatch error and cause no _E or _0 suppression.
     * @param derivOrder    Final features will be augmented with their 
     *                      derivatives up to 'derivOrder' order. If 'derivOrder'
     *                      is negative value, no new derivatives are appended 
     *                      and derivatives that already present in feature file
     *                      are preserved.  Straight features are considered 
     *                      to be of zero order. If some derivatives are already 
     *                      present in feature file, these are not computed 
     *                      again, only higher order derivatives are appended 
     *                      if required. Note, that HTK feature file cannot 
     *                      contain higher order derivatives (e.g. double delta)
     *                      without containing lower ones (e.g. delta). 
     *                      Derivative present in feature file that are of 
     *                      higher order than is required are discarded.  
     *                      Derivatives are computed in the final stage from 
     *                      (extracted segment of) feature frames possibly 
     *                      extended by repeated frames. Derivatives are 
     *                      computed using the same formula that is employed 
     *                      also by HTK tools. Lengths of windows used for 
     *                      computation of derivatives are passed in parameter 
     *                      derivWinLen. To compute derivatives for frames close 
     *                      to boundaries, frames before the first and after the 
     *                      last frame (of the extracted segment) are considered 
     *                      to be (yet another) repetitions of the first and the 
     *                      last frame, respectively. If the segment of frames 
     *                      is extracted according to range specification and 
     *                      parameters extLeft and extLeft are set to zero, the 
     *                      first and the last frames of the segment are 
     *                      considered to be repeated, eventough the true feature
     *                      frames from beyond the segment boundary can be
     *                      available in the file. Therefore, segment extracted 
     *                      from features that were before augmented with 
     *                      derivatives will differ 
     *                      from the same segment augmented with derivatives by 
     *                      this function. Difference will be of course only on 
     *                      boundaries and only in derivatives. This "incorrect" 
     *                      behavior was chosen to fully simulate behavior of 
     *                      HTK tools. To obtain more correct computation of 
     *                      derivatives, use parameters extLeft and extRight, 
     *                      which correctly extend segment with the true frames 
     *                      (if possible) and in resulting feature matrix ignore 
     *                      first extLeft and last extRight frames. For this 
     *                      purpose, both extLeft and extRight should be set to 
     *                      sum of all values in the array derivWinLen.
     * @param pDerivWinLen  Array of size derivOrder specifying lengths of 
     *                      windows used for computation of derivatives. 
     *                      Individual values represents one side context 
     *                      used in the computation. The each window length is 
     *                      therefore twice the value from array plus one. 
     *                      Value at index zero specify window length for first 
     *                      order derivatives (delta), higher indices 
     *                      corresponds to higher order derivatives.
     * @param pCmnPath      Cepstral mean normalization path
     * @param pCmnMask      Cepstral mean normalization mask
     * @param pCvnPath      Cepstral variance normalization path
     * @param pCvnMask      Cepstral variance normalization mask
     * @param pCvgFile      Global variance file to be parsed
     *
     * The given parameters are necessary for propper feature extraction 
     */
    void
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
        const char*           pCvgFile);
    
        
    /** 
     * @brief Returns a refference to the current file header
     */
    const HtkHeader&
    CurrentHeader() const 
    { return mHeader; }

    /** 
     * @brief Returns a refference to the current file header
     */
    const HtkHeaderExt&
    CurrentHeaderExt() const 
    { return mHeaderExt; }

    /**
     * @brief Returns the current file details
     *
     * @return Refference to a class @c FileListElem
     *
     * Logical and physical file names are stored in @c FileListElem class
     */
    const std::list<FileListElem>::iterator&
    pCurrentRecord() const
    { return mInputQueueIterator; }


    /**
     * @brief Returns the following file details
     *
     * @return Refference to a class @c FileListElem
     *
     * Logical and physical file names are stored in @c FileListElem class
     */
    const std::list<FileListElem>::iterator&
    pFollowingRecord() const
    { return mInputQueueIterator; }


    void
    Rewind()
    { mInputQueueIterator = mInputQueue.begin(); }
    
    
    /**
     * @brief Adds a single feature file to the repository
     * @param rFileName file to read features from
     */
    void
    AddFile(const std::string & rFileName);
    

    /**
     * @brief Adds a list of feature files to the repository
     * @param rFileName feature list file to read from
     */
    void
    AddFileList(const char* pFileName, const char* pFilter = "");
  
    
    const FileListElem&
    Current() const
    { return *mInputQueueIterator; }

    
    /** 
     * @brief Moves to the next record
     */
    void
    MoveNext();
    
    /**
     * @brief Reads full feature matrix from a feature file
     * @param rMatrix matrix to be created and filled with read data
     * @return number of successfully read feature vectors
     */
    bool
    ReadFullMatrix(Matrix<FLOAT>& rMatrix);    
    
    const size_t
    QueueSize() const {return mInputQueue.size(); }

    /**
     * @brief Reads feature vectors from a feature file
     * @param rMatrix matrix to be (only!) filled with read data. 
     * @return number of successfully read feature vectors
     * 
     * The function tries to fill @c pMatrix with feature vectors comming from
     * the current stream. If there are less vectors left in the stream, 
     * they are used and true number of successfuly read vectors is returned.
     */
    int
    ReadPartialMatrix(Matrix<FLOAT>& rMatrix);    
    

    /**
     * @brief Returns true if there are no feature files left on input
     */
    const bool
    EndOfList() const 
    { return mInputQueueIterator == mInputQueue.end(); }

    const std::string&
    CurrentIndexFileName() const
    { return mCurrentIndexFileName; }
    
    friend
    void
    AddFileListToFeatureRepositories(
    const char* pFileName,
    const char* pFilter,
    std::queue<FeatureRepository *> &featureRepositoryList);

////////////////////////////////////////////////////////////////////////////////
//  PRIVATE SECTION
////////////////////////////////////////////////////////////////////////////////
  private:
    /// List (queue) of input feature files
    std::list<FileListElem>             mInputQueue;
    std::list<FileListElem>::iterator   mInputQueueIterator;
    
    std::string                         mCurrentIndexFileName;
    std::string                         mCurrentIndexFileDir;
    std::string                         mCurrentIndexFileExt;

    /// current stream
    IStkStream                  mStream;
      
    // stores feature file's HTK header
    HtkHeader                   mHeader;
    HtkHeaderExt                mHeaderExt;


    // this group of variables serve for working withthe same physical
    // file name more than once
    char*                       mpLastFileName;
    std::string                 mLastFileName;
    char*                       mpLastCmnFile;
    char*                       mpLastCvnFile;
    char*                       mpLastCvgFile;
    FLOAT*                      mpCmn;
    FLOAT*                      mpCvn;
    FLOAT*                      mpCvg;
    HtkHeader                   mLastHeader;
    FLOAT*                      mpA;
    FLOAT*                      mpB;

#if FEATURE_REPOSITORY_PTHREAD_ENABLE
    // read overhead paralelism
    pthread_mutex_t*            mpRoMutex;
    pthread_cond_t*             mpRoCondNotFull;
    pthread_cond_t*             mpRoCondNotEmpty;
#endif 

    // Reads HTK feature file header
    int 
    ReadHTKHeader();

    int 
    ReadHTKFeature(FLOAT*    pIn, 
      size_t    feaLen, 
      bool      decompress, 
      FLOAT*    pScale, 
      FLOAT*    pBias);

    
    bool 
    ReadHTKFeatures(const std::string& rFileName, Matrix<FLOAT>& rFeatureMatrix);
    
    bool 
    ReadHTKFeatures(const FileListElem& rFileNameRecord, Matrix<FLOAT>& rFeatureMatrix);
  }; // class FeatureStream
}; //namespace STK

#endif // STK_Features_h
