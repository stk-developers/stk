#ifndef ELEMENT_H
#define ELEMENT_H

#include "../STKLib/Error.h"
#include "../STKLib/Matrix.h"

namespace SNet{
  class NNet;

  // Object for sending / receiving
  class Element{
    public: // everything public => it is only holding frame for send, receive
      Matrix<FLOAT> **mpWeights;
      Matrix<FLOAT> **mpBiases;
      int mNLayers;
      int mLast;
      int mFrom;
    
      Element(NNet *nnet, bool alocate);
      ~Element();
      void ReferenceUpdate(NNet *nnet);
      void Reference(NNet *nnet);
      void Clear();
      void Add(Element *element, float c);
  };
} // namespace
#endif
