#include "element.h"

SNet::Element::Element(NNet *nnet){
  mNLayers = nnet->NLayers();
  mpWeights = new Matrix<FLOAT>* [mNLayers];
  mpBiases = new Matrix<FLOAT>* [mNLayers];
  for(int i=0; i < mNLayers; i++){
    Matrix<FLOAT> *pomMatrix = nnet->Layers(i)->Weights();
    mpWeights[i] = new Matrix<FLOAT>(pomMatrix->Rows(), pomMatrix->Cols(), pomMatrix->Storage());
    pomMatrix = nnet->Layers(i)->Biases();
    mpBiases[i] = new Matrix<FLOAT>(pomMatrix->Rows(), pomMatrix->Cols(), pomMatrix->Storage());
  }
}

SNet::Element::~Element(){
  for(int i=0; i < mNLayers; i++){
    delete mpWeights[i];
    delete mpBiases[i];
  }
  delete [] mpWeights;
  delete [] mpBiases;
}
