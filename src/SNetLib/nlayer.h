#ifndef NLAYER_H
#define NLAYER_H

#include "../STKLib/fileio.h"
#include "../STKLib/common.h"
#include "../STKLib/Models.h"
#include "../STKLib/Viterbi.h"
#include "../STKLib/Matrix.h"

namespace SNet{
  /// One neural net layer (Consists of linear transformation, biases and non-linearity)
  class NLayer{
    private:
      int mOutFunc;              ///< Output function

      Matrix<FLOAT>* mpWeights;         ///< Matrix of weights (MxN, N is size of next layer, transp. in memory)
      Matrix<FLOAT>* mpBiases;          ///< Bias vector
      Matrix<FLOAT>* mpChangesWeights;  ///< Changes matrix of weights computed during forward pass
      Matrix<FLOAT>* mpChangesBiases;   ///< Changes bias vector computed during forward pass
      Matrix<FLOAT>* mpRecWeights;      ///< Pointer to weights matrixes received from clients
      Matrix<FLOAT>* mpRecBiases;       ///< Pointer to bias vectors received from clients
      Matrix<FLOAT>* mpIn;              ///< Matrix of layer input vectors (first layer uses window of input cache)
      Matrix<FLOAT>* mpOut;             ///< Matrix of layer output vectors (last layer uses window of input cache)
      Matrix<FLOAT>* mpError;           ///< Layer error matrix (errors before non-linearity)
    public:
      NLayer(Matrix<FLOAT>* weights, Matrix<FLOAT>* biases, int outFunc);
      void In(Matrix<FLOAT>* in){mpIn = in;};
      void Out(Matrix<FLOAT>* out){mpOut = out;};
      Matrix<FLOAT>* Out(){return mpOut;};
      Matrix<FLOAT>* Weights(){return mpWeights;};
      Matrix<FLOAT>* Biases(){return mpBiases;};
      
      void BunchBias();
      void BunchLinear();
      void BunchNonLinear();
  };
} // namespace

#endif
