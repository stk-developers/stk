#ifndef NNET_H
#define NNET_H

#include "nlayer.h"
#include "timers.h"
#include "socketobj.h"
#include<cstdlib>
#include<queue>

#define DEBUG_PROG 0

namespace SNet{
  //! Main neural net class
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
      float *mpLearningRateList;      ///< Neural net training rates
      bool mCrossValidation; ///< Only compute cross-validation
      int mNCache;           ///< Actual cache number
      
      Timers *mpTimers;                         ///< Help timers/counters
      STK::WindowMatrix<FLOAT>* mpInCache;           ///< Matrix of input vectors
      STK::WindowMatrix<FLOAT>* mpOutCache;          ///< Matrix of example output vectors
      STK::Matrix<FLOAT>* mpCompCachePart;           ///< Matrix of computed output vectors
      STK::Matrix<FLOAT>* mpError;                   ///< Neural net global output error
      NLayer** mpLayers;                        ///< Pointer to neural net layers
      Element *mpUpdateElement;                 ///< Pointer to element used for update matrixes which are sent to server
      Socket::Client *mpClient;                 ///< Client - used for sending updates
      std::queue<Element*> *mpReceivedElements; ///< Pointer to queue of received elements of new weights (from progobj)
      std::queue<Element*> *mpFreeElements;     ///< Pointer to queue of empty elements of new weights (from progobj)
      pthread_mutex_t *mpFreeMutex;             ///< Pointer to mutex for queue of free elements (from progobj)
      pthread_mutex_t *mpReceivedMutex;         ///< Pointer to mutex for queue of received elements (from progobj)
      barrier_t *mpBarrier;                     ///< Pointer to barrier (from progobj)
      bool *mpSync;                             ///< True if synchronization enabled (from progobj)
      
      void RandomizeIndices(int *randind, int n); ///< Makes random list 0,1,2,3...
    public:
      NNet(STK::CompositeXform* nn, int cacheSize, int bunchSize, bool crossValidation, ///< Constructor
           float *learningRateList); 
      ~NNet();                                                                     ///< Destructor
      
      void AddToCache(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize); ///< Copy input vector to cache
      bool CacheFull() const {return (mActualCache == mCacheSize);};               ///< Test if cache full
      void RandomizeCache();                                                       ///< Randomly shuffle cache     
      void ComputeCache(bool last);                                                ///< Compute whole cache
      void ComputeBunch();                                                         ///< Compute one bunch
      void GetAccuracy();                                                          ///< Find how many vectors in bunch are accurate
      void ComputeGlobalError();                                                   ///< Compute error on last layer
      void ComputeUpdates();                                                       ///< Compute weights changes using back-propagation
      void ChangeWeights();                                                        ///< Change weights
      void PrintInfo();                                                            ///< Print informations about training/cross-validation
      void PrepareUpdateElement();                                                 ///< Creates new element and connects it to update matrixes
      void ChangeToElement(Element *element);                                      ///< Changes NN using element
      void WaitForStartingWeights();                                               ///< Wait, receive and change to starting weights from server
      
      // Accessors
      bool CrossValidation() const {return mCrossValidation;};
      int NLayers()          const {return mNLayers;};
      NLayer* Layers(int i)  const {return mpLayers[i];};
      float LearningRate(int i)      const {return mpLearningRateList[i];};
      float* LearningRate()      const {return mpLearningRateList;};
      int Good()             const {return mGood;};
      int Vectors()          const {return mVectors;};
      int Discarded()        const {return mDiscarded;};
      Timers* TimersGet()       const {return mpTimers;};
      void Client(Socket::Client *client){mpClient = client;};
      void ReceivedElements(std::queue<Element*> *receivedElements){mpReceivedElements = receivedElements;};
      void FreeElements(std::queue<Element*> *freeElements){mpFreeElements = freeElements;};
      void Mutexes(pthread_mutex_t *free, pthread_mutex_t *received, barrier_t *barrier, bool *sync){ ///< Sets sync objects  
                   mpFreeMutex = free; mpReceivedMutex = received; mpBarrier = barrier; mpSync = sync;}; 
      void Good(int good){mGood = good;};
      void Vectors(int vectors){mVectors = vectors;}
      void Discarded(int discarded){mDiscarded = discarded;};
  };
} // namespace
#endif
