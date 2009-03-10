#include "element.h"
#include "nnet.h"

using namespace STK;

SNet::Element::Element(NNet *nnet, bool allocate){
  mNLayers = nnet->NLayers();
  mLast = 0;
  mpWeights = new Matrix<FLOAT>* [mNLayers];
  mpBiases = new Matrix<FLOAT>* [mNLayers];
  if(allocate){
    for(int i=0; i < mNLayers; i++){
      Matrix<FLOAT> *pomMatrix = nnet->Layers(i)->Weights();
      mpWeights[i] = new Matrix<FLOAT>(pomMatrix->Rows(), pomMatrix->Cols());
      pomMatrix = nnet->Layers(i)->Biases();
      mpBiases[i] = new Matrix<FLOAT>(pomMatrix->Rows(), pomMatrix->Cols());
    }
  } 
  else{
    for(int i=0; i < mNLayers; i++){
      mpWeights[i] = NULL;
      mpBiases[i] = NULL;
    } 
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

void SNet::Element::ReferenceUpdate(NNet *nnet){
  for(int i=0; i < mNLayers; i++){
    assert(mpWeights[i] == NULL);
    assert(mpBiases[i] == NULL);
    mpWeights[i] = nnet->Layers(i)->ChangesWeights();
    mpBiases[i] = nnet->Layers(i)->ChangesBiases();
  }
}

void SNet::Element::Reference(NNet *nnet){
  for(int i=0; i < mNLayers; i++){
    assert(mpWeights[i] == NULL);
    assert(mpBiases[i] == NULL);
    mpWeights[i] = nnet->Layers(i)->Weights();
    mpBiases[i] = nnet->Layers(i)->Biases();
  }
}

void SNet::Element::Clear(){
  for(int i=0; i < mNLayers; i++){
    assert(mpWeights[i] != NULL);
    assert(mpBiases[i] != NULL);
    mpWeights[i]->Clear();
    mpBiases[i]->Clear();
  }
}

void SNet::Element::Add(Element *element, float c){
  assert(element != NULL);
  assert(this->mNLayers == element->mNLayers);
  for(int i=0; i < mNLayers; i++){
    assert(mpWeights[i] != NULL);
    assert(mpBiases[i] != NULL);
    mpWeights[i]->AddCMMul(c, *(element->mpWeights[i]));
    mpBiases[i]->AddCMMul(c, *(element->mpBiases[i]));
  }
}

void SNet::Element::SubByLearningRate(Element *element, float *c){
  assert(element != NULL);
  assert(this->mNLayers == element->mNLayers);
  for(int i=0; i < mNLayers; i++){
    assert(mpWeights[i] != NULL);
    assert(mpBiases[i] != NULL);
    mpWeights[i]->AddCMMul(-c[i], *(element->mpWeights[i]));
    mpBiases[i]->AddCMMul(-c[i], *(element->mpBiases[i]));
  }
}
