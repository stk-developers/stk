#include "nnet.h"
#include "barrier.h"

using namespace STK;
SNet::NNet::NNet(CompositeXform* nn, int cacheSize, int bunchSize, bool crossValidation, float *learningRateList){
  if(nn->mNLayers % 3 != 0) // linear / biases / non-linear
    Error("NN has to have 3 Xform layers for one NN layer");
  mNLayers = nn->mNLayers / 3; // true NN layers consisting of 3 STK layers
  mpLayers = new NLayer*[mNLayers]; // create new layers 
  
  // Copy STK matrixes to SNet
  for(int i=0; i<mNLayers; i++){ 
    if(nn->mpLayer[3*i].mpBlock[0]->mXformType != XT_LINEAR)
      Error("First Xform in each NN layer has to be LinearXform");
    if(nn->mpLayer[3*i+1].mpBlock[0]->mXformType != XT_BIAS)
      Error("First Xform in each NN layer has to be BiasXform");
    if(nn->mpLayer[3*i+2].mpBlock[0]->mXformType != XT_FUNC)
      Error("First Xform in each NN layer has to be FuncXform");
    
    static_cast<LinearXform*>(nn->mpLayer[3*i].mpBlock[0])->mPredefinedID = PLTID_NONE;  // do not save this transform as predefined macro
    
    mpLayers[i] = new NLayer(&(static_cast<LinearXform*>(nn->mpLayer[3*i].mpBlock[0]))->mMatrix,
                             &(static_cast<BiasXform*>(nn->mpLayer[3*i+1].mpBlock[0]))->mVector,
                              (static_cast<FuncXform*>(nn->mpLayer[3*i+2].mpBlock[0]))->mFuncId);     
  }
  
  mBunchSize = bunchSize;
  mCacheSize = cacheSize;
  mCrossValidation = crossValidation;
  mpLearningRateList = learningRateList;
  mActualCache = 0;
  mVectors = 0;
  mDiscarded = 0;
  mGood = 0;
  mNCache = 0;  
  mpUpdateElement = NULL;
  
  // Create cache
  int in_cols =  (static_cast<LinearXform*>(nn->mpLayer[0].mpBlock[0]))->mMatrix.Cols();              // get number of columns
  int out_cols = (static_cast<LinearXform*>(nn->mpLayer[3*(mNLayers-1)].mpBlock[0]))->mMatrix.Rows(); // for cache size
  mpInCache =  new WindowMatrix<FLOAT>(cacheSize, in_cols);
  mpOutCache = new WindowMatrix<FLOAT>(cacheSize, out_cols);
   
  // Create partial cache for forwarded vectors
  mpCompCachePart = new Matrix<FLOAT>(bunchSize, mpLayers[mNLayers-1]->Weights()->Rows());
  
  // Create matrix for global NN errors
  mpError = new Matrix<FLOAT>(bunchSize, mpLayers[mNLayers-1]->Weights()->Rows());
  
  // Connect layers using their inputs, outputs and errors
  // If you would like to make not full weight connections, you should have separate weights
  for(int i=0; i<mNLayers; i++){
    mpLayers[i]->In((i == 0)           ? mpInCache  : mpLayers[i-1]->Out());
    mpLayers[i]->Out((i == mNLayers-1) ? mpCompCachePart : new Matrix<FLOAT>(bunchSize, mpLayers[i]->Weights()->Rows()));
    mpLayers[i]->Err((i == mNLayers-1) ? mpError : new Matrix<FLOAT>(bunchSize, mpLayers[i]->Weights()->Rows()));
  }
  
  // Get next layer's weights and errors - for error propagation
  for(int i=0; i<mNLayers; i++){
    mpLayers[i]->NextErr((i == mNLayers-1)     ? NULL : mpLayers[i+1]->Err());
    mpLayers[i]->NextWeights((i == mNLayers-1) ? NULL : mpLayers[i+1]->Weights());
  } 
  
  // Timers initialize
  mpTimers = new Timers(1, 3);
}

SNet::NNet::~NNet(){
  // Delete alocated matrixes - not last ones 
  // Also delete layers
  for(int i=0; i<mNLayers; i++){
    if(i != mNLayers-1){
      delete mpLayers[i]->Out();
      delete mpLayers[i]->Err();
    }
    delete mpLayers[i];
  }
  delete[] mpLayers;
  
  // Delete global matrixes
  delete mpInCache;
  delete mpOutCache;
  delete mpCompCachePart;
  delete mpError;
  delete mpTimers;
}

void SNet::NNet::AddToCache(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize){
  memcpy((*mpInCache)[mActualCache], inVector, inSize * sizeof(FLOAT));
  memcpy((*mpOutCache)[mActualCache], outVector, outSize * sizeof(FLOAT));
  mActualCache++;
}

int CmpRandRand (const void *dummy1, const void *dummy2) {
  double a = drand48();
  double b = drand48();
  if (a>b) return 1;
  else return -1;
}

void SNet::NNet::RandomizeIndices(int *randind, int n) {
  int i;
  for (i=0; i < n; i++)
    randind[i] = i;
  qsort((void*)randind, n, sizeof(int), &CmpRandRand);
}

void SNet::NNet::RandomizeCache(){  
  int randind[mActualCache];
  RandomizeIndices(randind, mActualCache); // make randomized list
  
  // Copy matrixes
  Matrix<FLOAT> old_in = *mpInCache;
  Matrix<FLOAT> old_out = *mpOutCache;

  // Move elements
  for(int i=0; i < mActualCache; i++){
    memcpy((*mpInCache)[i], old_in.Row(randind[i]), mpInCache->Cols() * sizeof(FLOAT)); // element i have to be at possition randind[i]
    memcpy((*mpOutCache)[i], old_out.Row(randind[i]), mpOutCache->Cols() * sizeof(FLOAT));
  }
}

void SNet::NNet::ComputeCache(bool last){
  std::cout << "Computing cache #" << mNCache << " with " << mActualCache << " vectors ... " << std::flush;
  mActualNOfBunch = mActualCache / mBunchSize; // number of full bunches 
  mDiscarded += mActualCache % mBunchSize;     // remains will be discarded
  
  if(mpUpdateElement != NULL){
    if(mActualNOfBunch == 0){ 
      if(DEBUG_PROG) std::cerr << "Problem SOLVED - cache with 0 bunches - should write only in Parallel version!";
      mpUpdateElement->mLast = 1;
      mpClient->SendElement(mpUpdateElement);
      if(DEBUG_PROG) if(*mpSync) std::cerr << "Waiting on barrier\n";
      if(*mpSync) barrier_wait(mpBarrier);
    }
  }
  
  // Compute all bunches
  for(int i=0; i < mActualNOfBunch; i++){
    mpInCache->SetSize(i*mBunchSize, mBunchSize);  // set window matrix
    mpOutCache->SetSize(i*mBunchSize, mBunchSize); //  --||--
    ComputeBunch(); // compute this bunch using window matrixes
    GetAccuracy(); // compute how good training or cross-validation is 
    if(!mCrossValidation){    
      ComputeGlobalError(); // compute global error - last layer
      ComputeUpdates(); // back-propagation
      if(mpUpdateElement != NULL){ // server-client, not 1 CPU
        if(i == mActualNOfBunch-1 && last){ // last cache & last bunch
          mpUpdateElement->mLast = 1;
        }
        mpClient->SendElement(mpUpdateElement);
        if(mpUpdateElement->mLast == 1){
          *mpSync = false;
          if(DEBUG_PROG) if(mpBarrier->counter == 1) std::cerr << "Waiting on barrier\n";
          if(mpBarrier->counter == 1) barrier_wait(mpBarrier);
        }
        mpTimers->Count(1);
        if(DEBUG_PROG) {
          if(mpUpdateElement->mLast == 1) std::cerr << "Element sent "<< mpTimers->Counter(1) <<" LAST\n";
          else std::cerr << "Element sent "<< mpTimers->Counter(1) << "\n";
        }
        if(DEBUG_PROG) if(*mpSync) std::cerr << "Waiting on barrier\n";
        if(*mpSync) barrier_wait(mpBarrier);
      }
      ChangeWeights(); // update 
    }
    mVectors += mBunchSize; // number of used vectors increased
  }
  
  // Reset - important!
  mpInCache->Reset();
  mpOutCache->Reset();
  mActualCache = 0;
  
  mNCache++;
  std::cout << "DONE! \n" << std::flush;
  if(DEBUG_PROG) std::cout << "GOOD " << mGood << "\n" << std::flush;
}

void SNet::NNet::ComputeBunch(){
  for(int i=0; i < mNLayers; i++){
    mpLayers[i]->BunchBias();
    mpLayers[i]->BunchLinear();
    mpLayers[i]->BunchNonLinear();
    // If you would like to make not full weight connections, you should copy results of one layer to another
  }
}

void SNet::NNet::GetAccuracy(){
  assert(mpCompCachePart->Rows() == mpOutCache->Rows());
  assert(mpCompCachePart->Cols() == mpOutCache->Cols());
  
  // Looking for maxs in rows of computed and example vectors
  for(unsigned r=0; r < mpCompCachePart->Rows(); r++){
    FLOAT* a1 = (*mpCompCachePart)[r];
    FLOAT* a2 = (*mpOutCache)[r];
    int maxPos1 = 0;
    int maxPos2 = 0;
    
    // Look for maximum in row
    for(unsigned i = 1; i < mpCompCachePart->Cols(); i++){
      if(a1[i] > a1[maxPos1]){
        maxPos1 = i;
      }
      if(a2[i] > a2[maxPos2]){
        maxPos2 = i;
      }
    }
    
    if(maxPos1 == maxPos2) {
      mGood++;
    }
  }
}

void SNet::NNet::ComputeGlobalError(){
  // Matrix computation -- E = C - O
  mpError->RepMMSub(*mpCompCachePart, *mpOutCache);
}

void SNet::NNet::ComputeUpdates(){
  for(int i=mNLayers-1; i >= 0; i--){ // go from last layer - back-propagation
    mpLayers[i]->ErrorPropagation();
    mpLayers[i]->ComputeLayerUpdates(); // needs previous layer computed - not a problem
  }
}
int pocitadlo = 0;
void SNet::NNet::ChangeWeights(){
  if(mpUpdateElement == NULL){ // 1 CPU, not server-client
    for(int i=0; i < mNLayers; i++){
      mpLayers[i]->ChangeLayerWeights(mpLearningRateList[i]); // learning rate needed
    }  
  }
  else{ // server-client, not 1 CPU
    Element *element = NULL;
    int size = 0;
    pthread_mutex_lock(mpReceivedMutex);
     size = mpReceivedElements->size();
    pthread_mutex_unlock(mpReceivedMutex);
    
    // If there are new weights
    if(size > 0){
      if(DEBUG_PROG) std::cerr << "Have new weights\n";
      pthread_mutex_lock(mpReceivedMutex);
      while(mpReceivedElements->size() > 0){ // get last weights in queue
        element = mpReceivedElements->front();
        mpReceivedElements->pop();
        if(mpReceivedElements->size() > 0){
          pthread_mutex_lock(mpFreeMutex);
           mpFreeElements->push(element);
          pthread_mutex_unlock(mpFreeMutex);
        }
      }
      pthread_mutex_unlock(mpReceivedMutex);
      this->ChangeToElement(element); // change NN using last received weights
      if(DEBUG_PROG) std::cerr << "Weights changed\n";
      mpTimers->Count(2);
      pthread_mutex_lock(mpFreeMutex);
      mpFreeElements->push(element);
      pthread_mutex_unlock(mpFreeMutex);
    }
    
  }
}

void SNet::NNet::PrintInfo(){
  std::cout << "===== SNET FINISHED (" << mpTimers->Timer(0) << "s) ===== \n";
  if(mCrossValidation)
    std::cout << "-- CV correct: >> ";
  else 
    std::cout << "-- TR correct: >> ";
  std::cout << 100.0*mGood / mVectors << "% << (Vectors " << mVectors << ", Good " << mGood << ", Discarded " << mDiscarded;
  if(mpUpdateElement != NULL){
    std::cout << " SRR=" << (float) mpTimers->Counter(2) / mpTimers->Counter(1);
  }
  
  std::cout << ") \n";
  
  std::cout << "\n";
}

void SNet::NNet::PrepareUpdateElement(){
  mpUpdateElement = new Element(this, false);
  mpUpdateElement->ReferenceUpdate(this); // set element as reference to update matrixes
}

void SNet::NNet::ChangeToElement(Element *element){
  Matrix<FLOAT> *pom;
  for(int i=mNLayers-1; i >= 0; i--){
    pom = this->Layers(i)->Weights();
    this->Layers(i)->Weights(element->mpWeights[i]);
    element->mpWeights[i] = pom;
    pom = this->Layers(i)->Biases();
    this->Layers(i)->Biases(element->mpBiases[i]);
    element->mpBiases[i] = pom;
    
    // For every not-last layer change NextErr and NextWeights - REALLY IMPORTANT!!!
    if(i < mNLayers-1){
      this->Layers(i)->NextErr(this->Layers(i+1)->Err());
      this->Layers(i)->NextWeights(this->Layers(i+1)->Weights());
    }
    
  }
}

void SNet::NNet::WaitForStartingWeights(){
  Element *element;
  int size = 0;
  if(DEBUG_PROG) std::cerr << "Waiting for first weights\n";
  do{
    pthread_mutex_lock(mpReceivedMutex);
     size = mpReceivedElements->size();
    pthread_mutex_unlock(mpReceivedMutex);
    if(size != 0){
      if(DEBUG_PROG) std::cerr << "Have first weights\n";
      pthread_mutex_lock(mpReceivedMutex);
       element = mpReceivedElements->front();
       mpReceivedElements->pop();
      pthread_mutex_unlock(mpReceivedMutex);
      this->ChangeToElement(element); // change NN using last received weights
      if(DEBUG_PROG) std::cerr << "Weights changed to first\n";
      pthread_mutex_lock(mpFreeMutex);
       mpFreeElements->push(element);
      pthread_mutex_unlock(mpFreeMutex);
    }
  } while (size == 0);
  if(DEBUG_PROG) std::cerr << "End of waiting for first weights\n";
}
