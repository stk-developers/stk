#include"nlayer.h"

SNet::NLayer::NLayer(Matrix<FLOAT>* weights, Matrix<FLOAT>* biases, int outFunc){
  mpWeights = weights;
  mpBiases = biases;
  mOutFunc = gFuncTable[outFunc].KID; // STK struct - only function number is needed
  
  // Make weights and biases copy => update matrixes
  mpChangesWeights = new Matrix<FLOAT>(mpWeights->Rows(), mpWeights->Cols(), mpWeights->Storage());
  mpChangesBiases = new Matrix<FLOAT>(mpBiases->Rows(), mpBiases->Cols(), mpBiases->Storage());
}

SNet::NLayer::~NLayer(){
  delete mpChangesWeights;
  delete mpChangesBiases;
}

void SNet::NLayer::BunchBias(){
  // For every row precopy biases FIRST! (see BunchLinear for details)
  for(unsigned i=0; i < mpOut->Rows(); i++){  
    assert(mpBiases->Cols() == mpOut->Cols());
    memcpy(mpOut->Row(i), mpBiases->Row(0), mpBiases->Cols() * sizeof(FLOAT));
  }
}
      
void SNet::NLayer::BunchLinear(){
  // Matrix computation -- O += I * W  
  mpOut->AddMMTMul(*mpIn, *mpWeights);
  // Biases ready, added by BLAS operation - it is optimization
}
      
void SNet::NLayer::BunchNonLinear(){
  // Using fast exponentials
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
  // Matrix computation -- W += CW * l 
  mpWeights->AddMCMul(*mpChangesWeights, (FLOAT)-1.0*learnRate);
  // Matrix computation -- B += CB * l
  mpBiases->AddMCMul(*mpChangesBiases, (FLOAT)-1.0*learnRate);
}

void SNet::NLayer::ErrorPropagation(){
  if(mpNextErr != NULL){ // last layer error = global error => do not compute it
    // Matrix computation -- E = NE * NW 
    mpErr->RepMMMul(*mpNextErr, *mpNextWeights);
  }
  DerivateError(); // derivate error of this layer 
}

void SNet::NLayer::DerivateError(){
  if(mOutFunc == KID_SoftMax){
    // Nothing
  }
  else if(mOutFunc == KID_Sigmoid){
    FLOAT *pErr = NULL;
    FLOAT *pOut = NULL;
    
    // For every element of matrix compute  e = (1 - o) * o * e
    for(unsigned row = 0; row < mpErr->Rows(); row++){
      pErr = mpErr->Row(row);
      pOut = mpOut->Row(row);
      for(unsigned col = 0; col < mpErr->Cols(); col++){
	(*pErr) = (1 - (*pOut)) * (*pOut) * (*pErr);
	pOut++;
	pErr++;
      }
    }
    
  }
}

void SNet::NLayer::ComputeLayerUpdates(){
  // Matrix computation -- CW = E^T * I  
  mpChangesWeights->RepMTMMul(*mpErr, *mpIn);
  
  mpChangesBiases->Clear();
  FLOAT *pChangesBiases = NULL;
  FLOAT *pErr = NULL;
  
  // Biases changes are computed as sum of error rows
  for(unsigned r=0; r < mpErr->Rows(); r++){
    pChangesBiases = mpChangesBiases->Row(0);
    pErr = mpErr->Row(r);
    for(unsigned c=0; c < mpErr->Cols(); c++){
      (*pChangesBiases) += (*pErr);
      pChangesBiases++;
      pErr++;
    }
  } 
  
}
