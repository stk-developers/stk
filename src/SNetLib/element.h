#ifndef ELEMENT_H
#define ELEMENT_H

#include "../STKLib/Error.h"
#include "nnet.h"

namespace SNet{
  // Object for sending / receiving
  class Element{
    public: // everything public => it is only holding frame for send, receive
      Matrix<FLOAT> **mpWeights;
      Matrix<FLOAT> **mpBiases;
      int mNLayers;
      int mLast;
      int mFrom;
    
      Element(NNet *nnet);
      ~Element();
  };
} // namespace
#endif
