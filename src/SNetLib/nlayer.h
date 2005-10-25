#ifndef NLAYER_H
#define NLAYER_H

#include"nlayer.h"

typedef float* Matrix;
typedef float* WindowMatrix;

class NLayer{
  private:
    int outFunc;

    Matrix weights;
    Matrix biases;
    Matrix changesWeights;
    Matrix changesBiases;
    Matrix *recWeights;
    Matrix *recBiases;
    Matrix in;
    Matrix out;
    Matrix error;
  public:

};

#endif
