#include"progobj.h"
   
SNet::ProgObj::ProgObj(XformInstance *NNet_instance, int cache_size, int bunch_size, bool cross_validation, std::string version){
  mServer = 0;
  mClient = 0;
  mPort = 2020;
  mIp = "IP not set";
  mNoClients = 0;
  mpThreads = NULL; 
  mpMainSocket = NULL;
  mpClientSocket = NULL;
  
  if(NNet_instance->mpXform->mXformType !=  XT_COMPOSITE)
    Error("NN has to be CompositeXform");
  
  CompositeXform* nn = static_cast<CompositeXform*>(NNet_instance->mpXform);

  mpNNet = new NNet(nn, cache_size, bunch_size, cross_validation);
  std::cout << "===== SNET v" << version << " " << (cross_validation ? "CROSS-VALIDATION" : "TRAINING") << " STARTED ===== \n";  
}

SNet::ProgObj::~ProgObj(){
  delete mpNNet;
}

void SNet::ProgObj::NewVector(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize, bool last){
  mpNNet->AddToCache(inVector, outVector, inSize, outSize);
  if(last || mpNNet->CacheFull()){
    if(!mpNNet->CrossValidation()) 
      mpNNet->RandomizeCache();
    mpNNet->ComputeCache();
  }
  if(last){
    mpNNet->PrintInfo();
  }
}

