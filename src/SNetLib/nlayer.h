#ifndef NLAYER_H
#define NLAYER_H

// :TODO: Include real matrixes from STK
typedef float* Matrix;
typedef float* WindowMatrix;

namespace SNet{

/// One neural net layer (Consists of linear transformation, biases and non-linearity)
class NLayer{
  private:
    int mOutFunc;            ///< Output function

    Matrix mWeights;         ///< Matrix of weights (MxN where N is size of next layer, transposed in memory)
    Matrix mBiases;          ///< Bias vector
    Matrix mChangesWeights;  ///< Changes matrix of weights computed during forward pass
    Matrix mChangesBiases;   ///< Changes bias vector computed during forward pass
    Matrix* mpRecWeights;    ///< Pointer to weights matrixes received from clients
    Matrix* mpRecBiases;     ///< Pointer to bias vectors received from clients
    Matrix mIn;              ///< Matrix of layer input vectors (first layer uses window of input cache)
    Matrix mOut;             ///< Matrix of layer output vectors (last layer uses window of input cache)
    Matrix mError;           ///< Layer error matrix (errors before non-linearity)
  public:

};

} // namespace

#endif
