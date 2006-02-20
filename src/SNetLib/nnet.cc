#include"nnet.h"

SNet::NNet::NNet(CompositeXform* nn, int cacheSize, int bunchSize, bool crossValidation, float learningRate){
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
    mpLayers[i] = new NLayer(&(static_cast<LinearXform*>(nn->mpLayer[3*i].mpBlock[0]))->mMatrix,
                             &(static_cast<BiasXform*>(nn->mpLayer[3*i+1].mpBlock[0]))->mVector,
                              (static_cast<FuncXform*>(nn->mpLayer[3*i+2].mpBlock[0]))->mFuncId);     
  }
  
  mBunchSize = bunchSize;
  mCacheSize = cacheSize;
  mCrossValidation = crossValidation;
  mLearnRate = learningRate;
  mActualCache = 0;
  mVectors = 0;
  mDiscarded = 0;
  mGood = 0;
  mNCache = 0;  
  
  // Create cache
  int in_cols =  (static_cast<LinearXform*>(nn->mpLayer[0].mpBlock[0]))->mMatrix.Rows();              // get number of columns
  int out_cols = (static_cast<LinearXform*>(nn->mpLayer[3*(mNLayers-1)].mpBlock[0]))->mMatrix.Cols(); // for cache size
  mpInCache =  new WindowMatrix<FLOAT>(cacheSize, in_cols);
  mpOutCache = new WindowMatrix<FLOAT>(cacheSize, out_cols);
   
  // Create partial cache for forwarded vectors
  mpCompCachePart = new Matrix<FLOAT>(bunchSize, mpLayers[mNLayers-1]->Weights()->Cols());
  
  // Create matrix for global NN errors
  mpError = new Matrix<FLOAT>(bunchSize, mpLayers[mNLayers-1]->Weights()->Cols());
  
  // Connect layers using their inputs, outputs and errors
  for(int i=0; i<mNLayers; i++){
    mpLayers[i]->In((i == 0)           ? mpInCache  : mpLayers[i-1]->Out());
    mpLayers[i]->Out((i == mNLayers-1) ? mpCompCachePart : new Matrix<FLOAT>(bunchSize, mpLayers[i]->Weights()->Cols()));
    mpLayers[i]->Err((i == mNLayers-1) ? mpError : new Matrix<FLOAT>(bunchSize, mpLayers[i]->Weights()->Cols()));
  }
  
  // Get next layer's weights and errors - for error propagation
  for(int i=0; i<mNLayers; i++){
    mpLayers[i]->NextErr((i == mNLayers-1)     ? NULL : mpLayers[i+1]->Err());
    mpLayers[i]->NextWeights((i == mNLayers-1) ? NULL : mpLayers[i+1]->Weights());
  } 
  
  // Timers initialize
  mpTimers = new Timers(1, 0);
  mpTimers->Start(0);
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
}

void SNet::NNet::AddToCache(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize){
  memcpy(mpInCache->Row(mActualCache), inVector, inSize * sizeof(FLOAT));
  memcpy(mpOutCache->Row(mActualCache), outVector, outSize * sizeof(FLOAT));
  mActualCache++;
}

void SNet::NNet::RandomizeCache(){
  Warning("Randomization not implemented!");
}

void SNet::NNet::ComputeCache(){
  std::cout << "Computing cache #" << mNCache << " with " << mActualCache << " vectors ... " << std::flush;
  mActualNOfBunch = mActualCache / mBunchSize; // number of full bunches 
  mDiscarded += mActualCache % mBunchSize;     // remains will be discarded
  
  // Compute all bunches
  for(int i=0; i < mActualNOfBunch; i++){
    mpInCache->SetSize(i*mBunchSize, mBunchSize);  // set window matrix
    mpOutCache->SetSize(i*mBunchSize, mBunchSize); //  --||--
    ComputeBunch(); // compute this bunch using window matrixes
    GetAccuracy(); // compute how good training or cross-validation is 
    if(!mCrossValidation){
      ComputeGlobalError(); // compute global error - last layer
      ComputeUpdates(); // back-propagation
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
}

void SNet::NNet::ComputeBunch(){
  for(int i=0; i < mNLayers; i++){
    mpLayers[i]->BunchBias();
    mpLayers[i]->BunchLinear();
    mpLayers[i]->BunchNonLinear();
  }
}

void SNet::NNet::GetAccuracy(){
  assert(mpCompCachePart->Rows() == mpOutCache->Rows());
  assert(mpCompCachePart->Cols() == mpOutCache->Cols());
  
  // Looking for maxs in rows of computed and example vectors
  for(unsigned r=0; r < mpCompCachePart->Rows(); r++){
    FLOAT* a1 = mpCompCachePart->Row(r);
    FLOAT* a2 = mpOutCache->Row(r);
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
  for(int i=mNLayers-1; i >= 0; i--){
    mpLayers[i]->ErrorPropagation();
  }
  for(int i=0; i < mNLayers; i++){
    mpLayers[i]->ComputeLayerUpdates();
  }
}

void SNet::NNet::ChangeWeights(){
  for(int i=0; i < mNLayers; i++){
    mpLayers[i]->ChangeLayerWeights(mLearnRate);
  }
}

void SNet::NNet::PrintInfo(){
  mpTimers->End(0);
  std::cout << "===== SNET FINISHED (" << mpTimers->Timer(0) << "s) ===== \n";
  if(mCrossValidation)
    std::cout << "-- CV correct: >> ";
  else 
    std::cout << "-- TR correct: >> ";
  std::cout << 100.0*mGood / mVectors << "% << (Vectors " << mVectors << ", Good " << mGood << ", Discarded " << mDiscarded << ") \n";
  std::cout << "\n";
}


