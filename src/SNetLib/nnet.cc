#include"nnet.h"

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

SNet::NNet::NNet(CompositeXform* nn, int cache_size, int bunch_size, bool cross_validation){
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
  
  /// Connect layers
  for(int i=0; i<mNLayers; i++){
    mpLayers[i]->In((i == 0)           ? mpInCache  : mpLayers[i-1]->Out());
    mpLayers[i]->Out((i == mNLayers-1) ? mpCompCachePart : new Matrix<FLOAT>(bunch_size, mpLayers[i]->Weights()->Cols()));
  } 
  
  mpTimers = new Timers(1, 0);
  mpTimers->Start(0);
}


SNet::NNet::~NNet(){
  for(int i=0; i<mNLayers; i++){
    if(i != mNLayers-1)
      delete mpLayers[i]->Out();
    delete mpLayers[i];
  }
  delete[] mpLayers;
  delete mpInCache;
  delete mpOutCache;
  delete mpCompCachePart;
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
  mDiscarded += mActualCache % mCacheSize;
  for(int i=0; i < mActualNOfBunch; i++){
    mpInCache->SetSize(i*mBunchSize, mBunchSize);
    mpOutCache->SetSize(i*mBunchSize, mBunchSize);
    ComputeBunch();
    GetAccuracy();
    mVectors += mBunchSize;
  }
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
  for(unsigned r=0; r < mpCompCachePart->Rows(); r++){
    if(FindMaxInVector(mpCompCachePart->Row(r), mpCompCachePart->Rows()) == FindMaxInVector(mpInCache->Row(r), mpInCache->Rows())){
      mGood++;
    }
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


