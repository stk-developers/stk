#ifndef NLAYER_H
#define NLAYER_H

#include "../STKLib/fileio.h"
#include "../STKLib/common.h"
#include "../STKLib/Models.h"
//#include "../STKLib/Decoder.h"
#include "../STKLib/Matrix.h"

namespace SNet{
  //! One neural net layer (Consists of linear transformation, biases and non-linearity)
  class NLayer{
    private:
      int mOutFunc; ///< Output function number

      STK::Matrix<FLOAT>* mpWeights;        ///< Matrix of weights (MxN, N is size of next layer, transp. in memory)
      STK::Matrix<FLOAT>* mpBiases;         ///< Bias vector
      STK::Matrix<FLOAT>* mpChangesWeights; ///< Changes matrix of weights computed during forward pass
      STK::Matrix<FLOAT>* mpChangesBiases;  ///< Changes bias vector computed during forward pass
      STK::Matrix<FLOAT>* mpRecWeights;     ///< Pointer to weights matrixes received from clients
      STK::Matrix<FLOAT>* mpRecBiases;      ///< Pointer to bias vectors received from clients
      
      STK::Matrix<FLOAT>* mpIn;             ///< Matrix of layer input vectors (first layer uses window of input cache)
      STK::Matrix<FLOAT>* mpOut;            ///< Matrix of layer output vectors (last layer uses window of input cache)
      // If you would like to make not full weight connections, you will need more matrixes, one for each block
      
      STK::Matrix<FLOAT>* mpErr;            ///< Layer error matrix (errors before non-linearity)
      STK::Matrix<FLOAT>* mpNextErr;        ///< Next layer errors - for error propagation
      STK::Matrix<FLOAT>* mpNextWeights;    ///< Next layer weights - for error propagation
    public:
      NLayer(STK::Matrix<FLOAT>* weights, STK::Matrix<FLOAT>* biases, int outFunc); ///< Constructor
      ~NLayer();                                                          ///< Destructor

      void BunchBias();                         ///< Compute bias transformations
      void BunchLinear();                       ///< Compute linear transformations
      void BunchNonLinear();                    ///< Compute non-linear transformations
      void ChangeLayerWeights(FLOAT learnRate); ///< Updates weights and biases
      void ErrorPropagation();                  ///< Back-propagates errors in net
      void DerivateError();                     ///< Derivates error
      void ComputeLayerUpdates();               ///< Computes update matrixes
      
      // Accessors
      STK::Matrix<FLOAT>* Out()            const {return mpOut;};
      STK::Matrix<FLOAT>* Err()            const {return mpErr;};
      STK::Matrix<FLOAT>* NextErr()        const {return mpNextErr;};
      STK::Matrix<FLOAT>* Weights()        const {return mpWeights;};
      STK::Matrix<FLOAT>* Biases()         const {return mpBiases;}; 
      STK::Matrix<FLOAT>* ChangesWeights() const {return mpChangesWeights;};
      STK::Matrix<FLOAT>* ChangesBiases()  const {return mpChangesBiases;};  
      STK::Matrix<FLOAT>* NextWeights()    const {return mpNextWeights;}; 
      void Weights(STK::Matrix<FLOAT>* weights)   {mpWeights = weights;};
      void Biases(STK::Matrix<FLOAT>* biases)   {mpBiases = biases;};
      void In(STK::Matrix<FLOAT>* in)   {mpIn = in;};
      void Out(STK::Matrix<FLOAT>* out) {mpOut = out;};
      void Err(STK::Matrix<FLOAT>* error) {mpErr = error;};
      void NextErr(STK::Matrix<FLOAT>* perror) {mpNextErr = perror;};
      void NextWeights(STK::Matrix<FLOAT>* pweights) {mpNextWeights = pweights;};
  };
} // namespace
#endif
