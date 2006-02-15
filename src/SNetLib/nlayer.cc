#include"nlayer.h"

SNet::NLayer::NLayer(Matrix<FLOAT>* weights, Matrix<FLOAT>* biases, int outFunc){
  mpWeights = weights;
  mpBiases = biases;
  mOutFunc = gFuncTable[outFunc].KID;
  mpChangesWeights = new Matrix<FLOAT>(mpWeights->Rows(), mpWeights->Cols(), mpWeights->Storage());
  mpChangesBiases = new Matrix<FLOAT>(mpBiases->Rows(), mpBiases->Cols(), mpBiases->Storage());
}

SNet::NLayer::~NLayer(){
  delete mpChangesWeights;
  delete mpChangesBiases;
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
    mpPreviousErr->RepMMTMul(*mpErr, *mpWeights);
    std::cout << "\n" << *mpPreviousErr;
  }
  
}

void SNet::NLayer::DerivateError(){
  if(mOutFunc == KID_SoftMax){
    // nothing
  }
  else if(mOutFunc == KID_Sigmoid){
    FLOAT *pErr = NULL;
    FLOAT *pOut = NULL;
    for(unsigned row = 0; row < mpErr->Rows(); row++){
      pErr = mpErr->Row(row);
      pOut = mpOut->Row(row);
      for(unsigned col = 0; col < mpErr->Cols(); col++){
        //(*mpErr)(row, col) = (1 - (*mpOut)(row, col)) * (*mpOut)(row, col) * (*mpErr)(row, col);
	(*pErr) = (1 - (*pOut)) * (*pOut) * (*pErr);
	pOut++;
	pErr++;
      }
    }
  }
}

void SNet::NLayer::ComputeLayerUpdates(){
  mpChangesWeights->RepMTMMul(*mpErr, *mpIn);
  mpChangesBiases->Clear();
  FLOAT *pChangesBiases = NULL;
  FLOAT *pErr = NULL;
  for(unsigned r=0; r < mpErr->Rows(); r++){
    pChangesBiases = mpChangesBiases->Row(0);
    pErr = mpErr->Row(r);
    for(unsigned c=0; c < mpErr->Cols(); c++){
      //(*mpChangesBiases)(0, c) += (*mpErr)(r, c);
      (*pChangesBiases) += (*pErr);
      pChangesBiases++;
      pErr++;
    }
  } 
}




