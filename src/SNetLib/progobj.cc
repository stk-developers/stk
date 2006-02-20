#include"progobj.h"
   
SNet::ProgObj::ProgObj(XformInstance *NNetInstance, int cacheSize, int bunchSize, bool crossValidation, 
                       std::string version, float learningRate){
  mServer = 0;
  mClient = 0;
  mPort = 2020;
  mIp = "IP not set";
  mNoClients = 0;
  mpThreads = NULL; 
  mpMainSocket = NULL;
  mpClientSocket = NULL;
  
  if(NNetInstance->mpXform->mXformType !=  XT_COMPOSITE)
    Error("NN has to be CompositeXform");
  
  CompositeXform* nn = static_cast<CompositeXform*>(NNetInstance->mpXform);

  mpNNet = new NNet(nn, cacheSize, bunchSize, crossValidation, learningRate); // create NN
  
  std::cout << "===== SNET v" << version << " " << (crossValidation ? "CROSS-VALIDATION" : "TRAINING") << " STARTED ===== \n";  
}

SNet::ProgObj::~ProgObj(){
  delete mpNNet;
}

void SNet::ProgObj::NewVector(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize, bool last){
  mpNNet->AddToCache(inVector, outVector, inSize, outSize); // make copy of new vector to cache
  
  // There is full or last cache
  if(last || mpNNet->CacheFull()){
    if(!mpNNet->CrossValidation()){
      mpNNet->RandomizeCache();
    }
    mpNNet->ComputeCache(); // compute full or partial cache
  }
  if(last){
    mpNNet->PrintInfo(); // print numbers of vectors
  }
}

