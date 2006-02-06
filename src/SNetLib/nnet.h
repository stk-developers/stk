#ifndef NNET_H
#define NNET_H

#include"nlayer.h"

namespace SNet{
  /// Main neural net class
  class NNet{
    private:
      int mNLayers;             ///< Number of neural net layers
      int mActualCache;          ///< Actual number of frames in cache
      int mCacheSize;            ///< Maximum cache size
      int mBunchSize;            ///< Size of one bunch - after bunch weights are updated
      int mActualNOfBunch;       ///< Actual number of bunches based on actual cache
      int mFrames;               ///< Number of training frames used for training so far
      int mGood;                 ///< Number of training frames which are classified as in example
      int mDicarded;             ///< Number of training frames discarded and not used
      int mLearnRate;            ///< Neural net training rate
      bool mCrossValidation;     ///< Only compute cross-validation

      WindowMatrix<FLOAT>* mpInCache;   ///< Matrix of input vectors
      WindowMatrix<FLOAT>* mpOutCache;  ///< Matrix of example output vectors
      Matrix<FLOAT>* mpCompCachePart;   ///< Matrix of computed output vectors
      Matrix<FLOAT>* mpError;           ///< Neural net global output error
      NLayer** mpLayers;                ///< Pointer to neural net layers
    public:
      NNet(CompositeXform* nn, int cache_size, int bunch_size, bool cross_validation);
      ~NNet();
      void AddToCache(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize);
      bool CacheFull() {return (mActualCache == mCacheSize);};
      void RandomizeCache();
      void ComputeCache();
      void ComputeBunch();
  };
} // namespace
#endif
