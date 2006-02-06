#include"nlayer.h"

SNet::NLayer::NLayer(Matrix<FLOAT>* weights, Matrix<FLOAT>* biases, int outFunc){
  mpWeights = weights;
  mpBiases = biases;
  mOutFunc = outFunc;
}
