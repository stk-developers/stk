#include"progobj.h"
   
SNet::ProgObj::ProgObj(XformInstance *NNetInstance, int cacheSize, int bunchSize, bool crossValidation, 
                       std::string version, float learningRate, int clients, char* ip){
  mPort = 2020;
  mNoClients = clients;
  
  if(clients == 0){
    mServer = false;
  }  
  else {
    mServer = true;
  }
  if(ip != NULL){
    mIp.assign(ip);
    mClient = true;
  }
  else{
    mIp = "IP not set";  
    mClient = false;
  }
  
  mpThreads = NULL; 
//  mpMainSocket = NULL;
//  mpClientSocket = NULL;
  mpServer = NULL;
  mpClient = NULL;  


  if(NNetInstance->mpXform->mXformType !=  XT_COMPOSITE)
    Error("NN has to be CompositeXform");
  
  CompositeXform* nn = static_cast<CompositeXform*>(NNetInstance->mpXform);

  mpNNet = new NNet(nn, cacheSize, bunchSize, crossValidation, learningRate); // create NN
  
  std::cout << "===== SNET v" << version << " " << (crossValidation ? "CROSS-VALIDATION" : "TRAINING") << " STARTED ===== \n";  
      
  srand48(GiveMeSeed());
}

SNet::ProgObj::~ProgObj(){
  delete mpNNet;
  if(mpThreads != NULL){
    delete [] mpThreads;
  }
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

int SNet::ProgObj::GiveMeSeed(){
  struct timeval tv;
  if (gettimeofday(&tv, 0) == -1) {
    assert(0 && "gettimeofday does not work.");
    exit(-1);
  }
  return (int)(tv.tv_sec) + (int)tv.tv_usec;
}

void* newServerThread(void *thdata){
  SNet::ThreadData *thread_data = (SNet::ThreadData*)thdata;
  thread_data->progObj->ServerReceivingThread(thread_data->number);

  return NULL;
}

void* newClientThread(void *thdata){
  SNet::ThreadData *thread_data = (SNet::ThreadData*)thdata;
  thread_data->progObj->ClientReceivingThread();

  return NULL;
}

void SNet::ProgObj::RunServer(){
  // Open server receiving threads
  mpThreads = new pthread_t[mNoClients]; // allocate thread pointers
  mpServer = new Socket::Server(mPort, mNoClients); // create server - will wait for clients
  
  // Create threads
  for(int i=0; i < mNoClients; i++){
    ThreadData *thread_data = new ThreadData;
    thread_data->number = i;
    thread_data->progObj = this;
    pthread_create(&(mpThreads[i]), NULL, newServerThread, (void*)&thread_data);
  }
  
    /// Main server thread
  // Get elements from queue
  // See if someone finished
  // Add to update
  // Make new element
  // Send to all clients
  // Loop unless all finished
}

void SNet::ProgObj::RunClient(){
  // Open client receiving thread
  mpThreads = new pthread_t[1];
  mpClient = new Socket::Client(mPort, (char*)mIp.c_str());

  // Create thread
  ThreadData *thread_data = new ThreadData;
  thread_data->progObj = this;
  pthread_create(&(mpThreads[0]), NULL, newClientThread, (void*)this);
  
  /// Main client thread
  // Do nothing, return to program
}

void SNet::ProgObj::ServerReceivingThread(int number){

}

void SNet::ProgObj::ClientReceivingThread(){

}
