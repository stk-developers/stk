#include"nlayer.h"

SNet::NLayer::NLayer(Matrix<FLOAT>* weights, Matrix<FLOAT>* biases, int outFunc){
  mpWeights = weights;
  mpBiases = biases;
  mOutFunc = gFuncTable[outFunc].KID;
}

void SNet::NLayer::BunchBias(){
  for(unsigned i=0; i < mpOut->Rows(); i++){
    assert(mpBiases->Cols() == mpOut->Cols());
    memcpy(mpOut->Row(i), mpBiases->Row(0), mpBiases->Cols() * sizeof(FLOAT));
  }
}
      
void SNet::NLayer::BunchLinear(){
  mpOut->AddMMMul(*mpIn, *mpWeights);
}
      
void SNet::NLayer::BunchNonLinear(){
  if(mOutFunc == KID_Sigmoid){
    mpOut->FastRowSigmoid();
  }
  else if(mOutFunc == KID_SoftMax){
    mpOut->FastRowSoftmax();
  }
  else {
    Error("This out function is not implemented. Only Sigmoid and SoftMax supported so far.");
  }
}

void SNet::NLayer::ChangeLayerWeights(FLOAT learnRate){
  mpWeights->AddMCMul(*mpChangesWeights, (FLOAT)-1.0*learnRate);
  mpBiases->AddMCMul(*mpChangesBiases, (FLOAT)-1.0*learnRate);
}

void SNet::NLayer::ErrorPropagation(){
  DerivateError();
  if(mpPreviousErr != NULL){
    mpPreviousErr->RepMMTMul(*mpOut, *mpWeights);
  }
}

void SNet::NLayer::DerivateError(){
  if(mOutFunc == KID_SoftMax){
    // nothing
  }
  else if(mOutFunc == KID_Sigmoid){
    for(unsigned row = 0; row < mpErr->Rows(); row++){
      for(unsigned col = 0; col < mpErr->Cols(); col++){
        (*mpErr)(row, col) = (1 - (*mpOut)(row, col)) * (*mpOut)(row, col) * (*mpErr)(row, col);
      }
    }
  }
}

void SNet::NLayer::ComputeLayerUpdates(){
  mpChangesWeights->RepMTMMul();
  // !!! rowsum a DNESKA DOM
  
}




