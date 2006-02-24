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
  

  mpFreeMutex = new pthread_mutex_t;
  mpReceivedMutex = new pthread_mutex_t;    
  pthread_mutex_init(mpFreeMutex, NULL);
  pthread_mutex_init(mpReceivedMutex, NULL);
  mpBarrier = new barrier_t;
  barrier_init(mpBarrier, 2);
}

SNet::ProgObj::~ProgObj(){
  delete mpNNet;
  if(mpThreads != NULL){
    delete [] mpThreads;
  }
  delete mpFreeMutex;
  delete mpReceivedMutex;  
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
     
bool SNet::ProgObj::Finished(int noClients){
  int finished = 0; 
  for(int i=0; i<mNoClients; i++){ 
    if(mpClientFinished[i]) finished++;
  }
  return(finished >= noClients);
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

SNet::Element * SNet::ProgObj::GetOrCreate(std::queue<Element*> *que, NNet *nn, pthread_mutex_t *mutex){
  Element *element;
  int size = 0;
  pthread_mutex_lock(mutex);
  size = que->size();
  if(size > 0){
    element = que->front();
    que->pop();
    pthread_mutex_unlock(mutex);
    return element;
  }
  else{  
    pthread_mutex_unlock(mutex);
    return(new Element(nn, true));
  }
}

/// ********** Server **********
void SNet::ProgObj::RunServer(){
  // Open server receiving threads
  
  mpThreads = new pthread_t[mNoClients]; // allocate thread pointers
  mpServer = new Socket::Server(mPort, mNoClients); // create server - will wait for clients
  mpClientFinished = new bool[mNoClients];
  for(int i=0; i<mNoClients; i++) mpClientFinished[i] = false;
  
  // Create threads
  for(int i=0; i < mNoClients; i++){
    ThreadData *thread_data = new ThreadData;
    thread_data->number = i;
    thread_data->progObj = this;
    pthread_create(&(mpThreads[i]), NULL, newServerThread, (void*)thread_data);
  }
  
  /// Main server thread
  Element *element_nn = new Element(mpNNet, false);
  element_nn->Reference(mpNNet);
  Element *element_add = new Element(mpNNet, true);
  Element *element;
  int size = 0;
 
  // Until every client finishes work
  while(!Finished(mNoClients)){
    
    // Work with received elements
    pthread_mutex_lock(mpReceivedMutex);
     size = mReceivedElements.size();
    pthread_mutex_unlock(mpReceivedMutex);
    if(size >= mNoClients){
      std::cerr << "Have enough " << size << "\n";
      element_add->Clear();
      for(int i=0; i < mNoClients; i++){
      //while(mReceivedElements.size() >= 0){
        pthread_mutex_lock(mpReceivedMutex);
         element = mReceivedElements.front();
	 mReceivedElements.pop();
	pthread_mutex_unlock(mpReceivedMutex);
	element_add->Add(element, 1.0);
        pthread_mutex_lock(mpFreeMutex);
	 mFreeElements.push(element);
	pthread_mutex_unlock(mpFreeMutex);
      }
      element_nn->Add(element_add, -1.0*mpNNet->LearnRate());
      element_nn->mLast = 0;
      mpServer->SendElementBroad(element_nn);
      std::cerr << "New weights sent\n";
    }
    
  }
}

void SNet::ProgObj::ServerReceivingThread(int number){
  Element *element;
  
  while(!Finished(mNoClients)){
    element = GetOrCreate(&mFreeElements, mpNNet, mpFreeMutex);
    mpServer->ReceiveElement(element, number);
    pthread_mutex_lock(mpReceivedMutex);
     mReceivedElements.push(element);
     std::cerr << "THREAD: element received\n";
    pthread_mutex_unlock(mpReceivedMutex);
    if(element->mLast == 1){
      mpClientFinished[number] = true;
    }
  } 
}

/// ********** Client **********
void SNet::ProgObj::RunClient(){
  // Open client receiving thread
  mpThreads = new pthread_t;
  mpClient = new Socket::Client(mPort, (char*)mIp.c_str());

  // Create thread
  ThreadData *thread_data = new ThreadData;
  thread_data->progObj = this;
  
  pthread_create(&(mpThreads[0]), NULL, newClientThread, (void*)thread_data);
  
  /// Main client thread
  mpNNet->PrepareUpdateElement();
  mpNNet->Client(mpClient);
  mpNNet->ReceivedElements(&mReceivedElements);
  mpNNet->FreeElements(&mFreeElements);
  mpNNet->Mutexes(mpFreeMutex, mpReceivedMutex, mpBarrier);
 
  // Do nothing, return to program
  // ::TODO:: check where monitor is needed
  // ::TODO:: check last
}

void SNet::ProgObj::ClientReceivingThread(){
  Element *element;
  mClientShouldFinish = false;
  while(!mClientShouldFinish){
    element = GetOrCreate(&mFreeElements, mpNNet, mpFreeMutex);
    mpClient->ReceiveElement(element);
    pthread_mutex_lock(mpReceivedMutex);
     mReceivedElements.push(element);
     std::cerr << "THREAD: Received Weights\n";
    pthread_mutex_unlock(mpReceivedMutex);
    barrier_wait(mpBarrier);
    if(element->mLast == 1){
      mClientShouldFinish = true;
    }
  }  
}
