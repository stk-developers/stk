#ifndef NNET_H
#define NNET_H

#include"nlayer.h"
#include "timers.h"

namespace SNet{
  // Main neural net class
  class NNet{
    private:
      int mNLayers;          ///< Number of neural net layers
      int mActualCache;      ///< Actual number of frames in cache
      int mCacheSize;        ///< Maximum cache size
      int mBunchSize;        ///< Size of one bunch - after bunch weights are updated
      int mActualNOfBunch;   ///< Actual number of bunches based on actual cache
      int mVectors;          ///< Number of training frames used for training so far
      int mGood;             ///< Number of training frames which are classified as in example
      int mDiscarded;        ///< Number of training frames discarded and not used
      float mLearnRate;      ///< Neural net training rate
      bool mCrossValidation; ///< Only compute cross-validation
      int mNCache;           ///< Actual cache number
      
      Timers *mpTimers;                 ///< Help timers/counters
      WindowMatrix<FLOAT>* mpInCache;   ///< Matrix of input vectors
      WindowMatrix<FLOAT>* mpOutCache;  ///< Matrix of example output vectors
      Matrix<FLOAT>* mpCompCachePart;   ///< Matrix of computed output vectors
      Matrix<FLOAT>* mpError;           ///< Neural net global output error
      NLayer** mpLayers;                ///< Pointer to neural net layers
    public:
      NNet(CompositeXform* nn, int cacheSize, int bunchSize, bool crossValidation, ///< Constructor
           float learningRate); 
      ~NNet();                                                                     ///< Destructor
      
      void AddToCache(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize); ///< Copy input vector to cache
      bool CacheFull() const {return (mActualCache == mCacheSize);};               ///< Test if cache full
      void RandomizeCache();                                                       ///< Randomly shuffle cache     
      void ComputeCache();                                                         ///< Compute whole cache
      void ComputeBunch();                                                         ///< Compute one bunch
      void GetAccuracy();                                                          ///< Find how many vectors in bunch are accurate
      void ComputeGlobalError();                                                   ///< Compute error on last layer
      void ComputeUpdates();                                                       ///< Compute weights changes using back-propagation
      void ChangeWeights();                                                        ///< Change weights
      void PrintInfo();                                                            ///< Prit informations about training/cross-validation
      
      // Accessors
      bool CrossValidation() const {return mCrossValidation;};
  };
} // namespace
#endif
