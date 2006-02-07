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
  mpOut->AddMatMult(*mpIn, *mpWeights);
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
