#ifndef ELEMENT_H
#define ELEMENT_H

#include "../STKLib/Matrix.h"

namespace SNet{
  class NNet;

  //! Object for sending / receiving NN weights and biases
  class Element{
    public: // everything public => it is only holding frame for send, receive
      int mNLayers;
      int mLast;
      int mFrom; // :TODO: can be deleted, not used right now
      
      STK::Matrix<FLOAT> **mpWeights;
      STK::Matrix<FLOAT> **mpBiases;
          
      Element(NNet *nnet, bool allocate); ///< Constructor, if not allocated, contains only pointers - Reference of ReferenceUpdate should be used
      ~Element();                        ///< Destructor
      
      void ReferenceUpdate(NNet *nnet);    ///< References element to update matrixes of NN
      void Reference(NNet *nnet);          ///< References element to main matrixes of NN
      void Clear();                        ///< Sets all matrixes to zeros
      void Add(Element *element, float c); ///< Adds element * constant
      void SubByLearningRate(Element *element, float *c);  ///< Subs element * constant[layer]
  };
} // namespace
#endif
