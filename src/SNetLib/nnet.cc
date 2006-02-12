#include"nnet.h"
/*
int SNet::NNet::FindMaxInVector(FLOAT *vector, int size){
  FLOAT *p = vector+1;
  FLOAT max = vector[0];
  int max_pos = 0;
  for(int i=1; i<size; i++){
    if(*p > max){
      max = *p;
      max_pos = i;
    }
  }
  return max_pos;
}
*/
SNet::NNet::NNet(CompositeXform* nn, int cache_size, int bunch_size, bool cross_validation, float learning_rate){
  if(nn->mNLayers % 3 != 0)
    Error("NN has to have 3 Xform layers for one NN layer");

  mNLayers = nn->mNLayers / 3;
  mpLayers = new NLayer*[mNLayers];  
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
  mBunchSize = bunch_size;
  mCacheSize = cache_size;
  mCrossValidation = cross_validation;
  mLearnRate = learning_rate;
  int in_cols =  (static_cast<LinearXform*>(nn->mpLayer[0].mpBlock[0]))->mMatrix.Rows();
  int out_cols = (static_cast<LinearXform*>(nn->mpLayer[3*(mNLayers-1)].mpBlock[0]))->mMatrix.Cols();
  mpInCache =  new WindowMatrix<FLOAT>(cache_size, in_cols);
  mpOutCache = new WindowMatrix<FLOAT>(cache_size, out_cols);
  mActualCache = 0;
  mVectors = 0;
  mDiscarded = 0;
  mGood = 0;
  mNCache = 0;
  
  mpCompCachePart = new Matrix<FLOAT>(bunch_size, mpLayers[mNLayers-1]->Weights()->Cols());
  mpError = new Matrix<FLOAT>(bunch_size, mpLayers[mNLayers-1]->Weights()->Cols());
  
  /// Connect layers
  for(int i=0; i<mNLayers; i++){
    mpLayers[i]->In((i == 0)             ? mpInCache  : mpLayers[i-1]->Out());
    mpLayers[i]->Out((i == mNLayers-1)   ? mpCompCachePart : new Matrix<FLOAT>(bunch_size, mpLayers[i]->Weights()->Cols()));
    mpLayers[i]->Err((i == mNLayers-1)   ? mpError : new Matrix<FLOAT>(bunch_size, mpLayers[i]->Weights()->Cols()));
    mpLayers[i]->PreviousErr((i == 0)    ? NULL : mpLayers[i-1]->Err());
  } 
  
  mpTimers = new Timers(1, 0);
  mpTimers->Start(0);
}


SNet::NNet::~NNet(){
  for(int i=0; i<mNLayers; i++){
    if(i != mNLayers-1){
      delete mpLayers[i]->Out();
      delete mpLayers[i]->Err();
    }
    delete mpLayers[i];
  }
  delete[] mpLayers;
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
  Error("Randomization not implemented!");
}

void SNet::NNet::ComputeCache(){
  std::cout << "Computing cache #" << mNCache << " with " << mActualCache << " vectors ... " << std::flush;
  mActualNOfBunch = mActualCache / mBunchSize;
  mDiscarded += mActualCache % mBunchSize;
  for(int i=0; i < mActualNOfBunch; i++){
    mpInCache->SetSize(i*mBunchSize, mBunchSize);
    mpOutCache->SetSize(i*mBunchSize, mBunchSize);
    ComputeBunch();
    GetAccuracy();
    if(!mCrossValidation){
      ComputeGlobalError();
      ComputeUpdates();
      ChangeWeights();
    }
      
    mVectors += mBunchSize;
  }
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
  for(unsigned r=0; r < mpCompCachePart->Rows(); r++){
    FLOAT* a1 = mpCompCachePart->Row(r);
    FLOAT* a2 = mpOutCache->Row(r);
    int maxPos1 = 0;
    int maxPos2 = 0;
    for(unsigned i = 1; i < mpCompCachePart->Cols(); i++){
      if(a1[i] > a1[maxPos1]){
        maxPos1 = i;
      }
      if(a2[i] > a2[maxPos2]){
        maxPos2 = i;
      }
    }
    if(maxPos1 == maxPos2) mGood++;
  }
}


void SNet::NNet::ComputeGlobalError(){
  mpError->RepMMSub(*mpCompCachePart, *mpOutCache);
}

void SNet::NNet::ComputeUpdates(){
  for(int i=mNLayers-1; i >= 0; i--){
    mpLayers[i]->ErrorPropagation();
  }
  for(int i=mNLayers-1; i >= 0; i--){
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


