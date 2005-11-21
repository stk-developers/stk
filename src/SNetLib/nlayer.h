#ifndef NLAYER_H
#define NLAYER_H

// :TODO: Include real matrixes from STK
typedef float* Matrix;
typedef float* WindowMatrix;

namespace SNet{
  /// One neural net layer (Consists of linear transformation, biases and non-linearity)
  class NLayer{
    private:
      int mOutFunc;              ///< Output function

      Matrix* mpWeights;         ///< Matrix of weights (MxN, N is size of next layer, transp. in memory)
      Matrix* mpBiases;          ///< Bias vector
      Matrix* mpChangesWeights;  ///< Changes matrix of weights computed during forward pass
      Matrix* mpChangesBiases;   ///< Changes bias vector computed during forward pass
      Matrix* mpRecWeights;      ///< Pointer to weights matrixes received from clients
      Matrix* mpRecBiases;       ///< Pointer to bias vectors received from clients
      Matrix* mpIn;           ///< Matrix of layer input vectors (first layer uses window of input cache)
      Matrix* mpOut;          ///< Matrix of layer output vectors (last layer uses window of input cache)
      Matrix* mpError;           ///< Layer error matrix (errors before non-linearity)
    public:

  };
} // namespace

#endif
