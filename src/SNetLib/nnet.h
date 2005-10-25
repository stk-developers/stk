#ifndef NNET_H
#define NNET_H

#include"nlayer.h"

class NNet{
  private:
    int nOfLayers;
    int actualCache;
    int cacheSize;
    int bunchSize;
    int actualNOfBunch;
    int frames;
    int good;
    int dicarded;
    int learnRate;

    WindowMatrix inCache;
    WindowMatrix outCache;
    WindowMatrix compCache;
    Matrix error;

    NLayer *layers;
  public:
  
};

#endif
