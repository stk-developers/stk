#include"nnet.h"

SNet::NNet::NNet(CompositeXform* nn, int cache_size, int bunch_size, bool cross_validation){
  if(nn->mNLayers % 3 != 0)
    Error("NN has to have 3 Xform layers for one NN layer");

  mNLayers = nn->mNLayers;
  mpLayers = new NLayer*[mNLayers];  
  for(int i=0; i<mNLayers; i++){
    if(nn->mpLayer[i].mpBlock[0]->mXformType != XT_LINEAR)
      Error("First Xform in each NN layer has to be LinearXform");
    if(nn->mpLayer[i+1].mpBlock[0]->mXformType != XT_BIAS)
      Error("First Xform in each NN layer has to be BiasXform");
    if(nn->mpLayer[i+2].mpBlock[0]->mXformType != XT_FUNC)
      Error("First Xform in each NN layer has to be FuncXform");
    mpLayers[i] = new NLayer(&(static_cast<LinearXform*>(nn->mpLayer[i].mpBlock[0]))->mMatrix,
                             &(static_cast<BiasXform*>(nn->mpLayer[i].mpBlock[0]))->mVector,
			     (static_cast<FuncXform*>(nn->mpLayer[i].mpBlock[0]))->mFuncId);     
  }
  mBunchSize = bunch_size;
  mCacheSize = cache_size;
  mCrossValidation = cross_validation;
  int in_cols =  (static_cast<LinearXform*>(nn->mpLayer[0].mpBlock[0]))->mMatrix.Cols();
  int out_cols = (static_cast<LinearXform*>(nn->mpLayer[mNLayers].mpBlock[0]))->mMatrix.Cols();
  mpInCache =  new WindowMatrix<FLOAT>(cache_size, in_cols);
  mpOutCache = new WindowMatrix<FLOAT>(cache_size, out_cols);
  mActualCache = 0;
}


SNet::NNet::~NNet(){
  for(int i=0; i<mNLayers; i++){
    delete mpLayers[i];
  }
  delete[] mpLayers;
  delete mpInCache;
  delete mpOutCache;
}

void SNet::NNet::AddToCache(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize){
  memcpy(mpInCache->Row(mActualCache), inVector, inSize * sizeof(FLOAT));
  memcpy(mpOutCache->Row(mActualCache), outVector, outSize * sizeof(FLOAT));
  mActualCache++;
}

void SNet::NNet::RandomizeCache(){

}

void SNet::NNet::ComputeCache(){
  mActualNOfBunch = mActualCache / mCacheSize;
  mDicarded += mActualCache % mCacheSize;
  for(int i=0; i < mActualNOfBunch; i++){
    mpInCache->SetSize(i*mBunchSize, mBunchSize);
    mpOutCache->SetSize(i*mBunchSize, mBunchSize);
  }
  ComputeBunch();
  mActualCache = 0;
}

void SNet::NNet::ComputeBunch(){

}

