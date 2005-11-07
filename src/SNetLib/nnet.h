#ifndef NNET_H
#define NNET_H

#include"nlayer.h"

namespace SNet{
  /// Main neural net class
  class NNet{
    private:
      int mNoLayers;             ///< Number of neural net layers
      int mActualCache;          ///< Actual number of frames in cache
      int mCacheSize;            ///< Maximum cache size
      int mBunchSize;            ///< Size of one bunch - after bunch weights are updated
      int mActualNOfBunch;       ///< Actual number of bunches based on actual cache
      int mFrames;               ///< Number of training frames used for training so far
      int mGood;                 ///< Number of training frames which are classified as in example
      int mDicarded;             ///< Number of training frames discarded and not used
      int mLearnRate;            ///< Neural net training rate

      WindowMatrix* mpInCache;   ///< Matrix of input vectors
      WindowMatrix* mpOutCache;  ///< Matrix of example output vectors
      WindowMatrix* mpCompCache; ///< Matrix of computed output vectors
      Matrix* mpError;           ///< Neural net global output error
      NLayer* mpLayers;          ///< Pointer to neural net layers
    public:

  };
} // namespace
#endif
