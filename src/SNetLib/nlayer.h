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
      int mOutFunc; ///< Output function number

      Matrix<FLOAT>* mpWeights;        ///< Matrix of weights (MxN, N is size of next layer, transp. in memory)
      Matrix<FLOAT>* mpBiases;         ///< Bias vector
      Matrix<FLOAT>* mpChangesWeights; ///< Changes matrix of weights computed during forward pass
      Matrix<FLOAT>* mpChangesBiases;  ///< Changes bias vector computed during forward pass
      Matrix<FLOAT>* mpRecWeights;     ///< Pointer to weights matrixes received from clients
      Matrix<FLOAT>* mpRecBiases;      ///< Pointer to bias vectors received from clients
      Matrix<FLOAT>* mpIn;             ///< Matrix of layer input vectors (first layer uses window of input cache)
      Matrix<FLOAT>* mpOut;            ///< Matrix of layer output vectors (last layer uses window of input cache)
      Matrix<FLOAT>* mpErr;            ///< Layer error matrix (errors before non-linearity)
      Matrix<FLOAT>* mpNextErr;        ///< Next layer errors - for error propagation
      Matrix<FLOAT>* mpNextWeights;    ///< Next layer weights - for error propagation
    public:
      NLayer(Matrix<FLOAT>* weights, Matrix<FLOAT>* biases, int outFunc); ///< Constructor
      ~NLayer();                                                          ///< Destructor

      void BunchBias();                         ///< Compute bias transformations
      void BunchLinear();                       ///< Compute linear transformations
      void BunchNonLinear();                    ///< Compute non-linear transformations
      void ChangeLayerWeights(FLOAT learnRate); ///< Updates weights and biases
      void ErrorPropagation();                  ///< Back-propagates errors in net
      void DerivateError();                     ///< Derivates error
      void ComputeLayerUpdates();               ///< Computes update matrixes
      
      // Accessors
      Matrix<FLOAT>* Out()     const {return mpOut;};
      Matrix<FLOAT>* Err()     const {return mpErr;};
      Matrix<FLOAT>* NextErr() const {return mpNextErr;};
      Matrix<FLOAT>* Weights() const {return mpWeights;};
      Matrix<FLOAT>* Biases()  const {return mpBiases;};      
      void In(Matrix<FLOAT>* in)   {mpIn = in;};
      void Out(Matrix<FLOAT>* out) {mpOut = out;};
      void Err(Matrix<FLOAT>* error) {mpErr = error;};
      void NextErr(Matrix<FLOAT>* perror) {mpNextErr = perror;};
      void NextWeights(Matrix<FLOAT>* pweights) {mpNextWeights = pweights;};
  };
} // namespace
#endif
