#include"progobj.h"
#include "barrier.h"

using namespace STK;
//******************************************************************************
//******************************************************************************
SNet::ProgObj::
ProgObj(STK::XformInstance* NNetInstance, int cacheSize, int bunchSize, 
    bool  crossValidation, std::string version, float learningRate, int clients, 
    char* ip, bool randomize, bool sync, int port, int seed, int ncumrbu, 
    char* learning_rate_list)
{
                       
  mpLearningRateList = NULL;
  int ls = (static_cast<CompositeXform*>(NNetInstance->mpXform))->mNLayers / 3;
  //mpLearningRateList = (float*)malloc(sizeof(float) * ls);
  mpLearningRateList = new float[ls];

  if(learning_rate_list != NULL)
  {
    // Parse learning_rate_list
    char* pch;
    pch = strtok(learning_rate_list, ",");
    int i=0;
    while (pch != NULL)
    {
      sscanf(pch, "%f", &mpLearningRateList[i]);
      pch = strtok (NULL, ",");
      mpLearningRateList[i] *= learningRate; // To have absolute learningRate multiplied by list
      i++;
    }
  }
  else
  {
    for(int i=0; i<ls; i++) 
    {
      mpLearningRateList[i] = learningRate;
    }
  }
                       
  mPort = port;
  mNoClients = clients;
  mRandomize = randomize;
  mSync = sync;
  
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
  mpServer = NULL;
  mpClient = NULL;  
  mNcumrbu = ncumrbu;

  if(NNetInstance->mpXform->mXformType !=  XT_COMPOSITE)
    Error("NN has to be CompositeXform");
  
  CompositeXform* nn = static_cast<CompositeXform*>(NNetInstance->mpXform);

  // Make neural network
  mpNNet = new NNet(nn, cacheSize, bunchSize, crossValidation, mpLearningRateList); // create NN
  
  std::cout << "===== SNET v" << version << " " << (crossValidation ? "CROSS-VALIDATION" : "TRAINING") << " STARTED ===== \n";  
  std::cout << "Learning rate list: ";
  for(int i=0; i<ls; i++) {
    std::cout << mpLearningRateList[i];
    if(i != ls-1)
      std::cout << ", ";
    else
      std::cout << "\n";
  }
      
  srand48((seed == 0) ? GiveMeSeed() : seed);
  
  // Create sync mechanisms
  mpFreeMutex       = new pthread_mutex_t;
  mpReceivedMutex   = new pthread_mutex_t;    
  pthread_mutex_init(mpFreeMutex, NULL);
  pthread_mutex_init(mpReceivedMutex, NULL);
  mpBarrier         = new barrier_t;
  mpEndBarrier      = new barrier_t;

  if(mClient) 
  {
    barrier_init(mpBarrier, 2); // receiving thread + main
    barrier_init(mpEndBarrier, 2); 
  }

  if(mServer) 
    barrier_init(mpBarrier, mNoClients+1); // all client threads + main
}


//******************************************************************************
//******************************************************************************
SNet::ProgObj::
~ProgObj()
{
  delete mpNNet;
  
  if(mpThreads != NULL)
    delete [] mpThreads;
  
  delete mpFreeMutex;
  delete mpReceivedMutex;  
  delete mpBarrier;
  delete mpEndBarrier;

  //free(mpLearningRateList);
  delete [] mpLearningRateList;
}


//******************************************************************************
//******************************************************************************
void 
SNet::ProgObj::
NewVector(FLOAT* pInVector, FLOAT* pOutVector, int inSize, int outSize, bool last)
{
  mpNNet->AddToCache(pInVector, pOutVector, inSize, outSize); // make copy of new vector to cache
  
  // There is full or last cache
  if(last || mpNNet->CacheFull()){
    if(!mpNNet->CrossValidation() && mRandomize){
      mpNNet->RandomizeCache();
    }
    mpNNet->ComputeCache(last); // compute full or partial cache
  }
  if(last){
    if(mClient){
      TimersGet()->End(0);
      if(DEBUG_PROG) std::cerr << "END OF CLIENT, WAIT ON END-BARRIER\n";
      barrier_wait(mpEndBarrier);
      
      mpClient->SendInt(mpNNet->Vectors());
      mpClient->SendInt(mpNNet->Good());
      mpClient->SendInt(mpNNet->Discarded());
      if(DEBUG_PROG) std::cerr << "Numbers sent!\n";
    }  
    if(!mClient && !mServer){
      TimersGet()->End(0); // 1 CPU 
    }
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
  for(int i=0; i < mNoClients; i++){ 
    if(mpClientFinished[i]) finished++;
  }
  return(finished >= noClients);
}

bool SNet::ProgObj::LastSent(int noClients){
  int finished = 0; 
  for(int i=0; i < mNoClients; i++){ 
    if(mpLastSent[i]) finished++;
  }
  return(finished >= noClients);
}

int 
SNet::ProgObj::
ActiveClients()
{
  int finished = 0; 
  
  for(int i=0; i < mNoClients; i++)
  { 
    if(mpClientFinished[i]) finished++;
  }
  return(mNoClients - finished);
}

void* 
newServerThread(void* thdata)
{
  //Just run receiving thread
  SNet::ThreadData* thread_data = (SNet::ThreadData*)thdata;
  thread_data->progObj->ServerReceivingThread(thread_data->number);

  return NULL;
}


void* 
newClientThread(void* thdata)
{
  //Just run receiving thread
  SNet::ThreadData* thread_data = (SNet::ThreadData*)thdata;
  thread_data->progObj->ClientReceivingThread();

  return NULL;
}


SNet::Element* 
SNet::ProgObj::
GetOrCreate(std::queue<Element*> *que, NNet *nn, pthread_mutex_t *mutex)
{
  Element*    element;
  int         size = 0;
  pthread_mutex_lock(mutex);
  size =      que->size();

  if(size > 0)
  {
    element = que->front();
    que->pop();
    pthread_mutex_unlock(mutex);
    return element;
  }
  else
  {  
    pthread_mutex_unlock(mutex);
    return(new Element(nn, true));
  }
}


//******************************************************************************
//******************************************************************************
//********** Server **********
void 
SNet::ProgObj::
RunServer()
{
  // Start server
  mpServer = new Socket::Server(mPort, mNoClients); // create server - will wait for clients
  
  // Set flags to default
  mpClientFinished = new bool[mNoClients];
  for(int i=0; i<mNoClients; i++) mpClientFinished[i] = false;
  mpLastSent = new bool[mNoClients];
  for(int i=0; i<mNoClients; i++) mpLastSent[i] = false;
  
  // Create and run threads
  mpThreads = new pthread_t[mNoClients]; // allocate thread pointers
  for(int i=0; i < mNoClients; i++){
    ThreadData *thread_data = new ThreadData;
    thread_data->number = i;
    thread_data->progObj = this;
    pthread_create(&(mpThreads[i]), NULL, newServerThread, (void*)thread_data);
  }
  
  /// Main server thread
  mpRecVectors = new int[mNoClients]; 
  mpRecGood = new int[mNoClients]; 
  mpRecDiscarded = new int[mNoClients];  
  Element *element_nn = new Element(mpNNet, false);
  element_nn->Reference(mpNNet);
  Element *element_add = new Element(mpNNet, true); // element for adding update matrixes
  Element *element;
  int size = 0;
  
  for(int i=0; i < mNoClients; i++){
    if(DEBUG_PROG) std::cerr << "Sending first weights " << i << "\n";
    mpServer->SendElement(element_nn, i); // at start send First Weights
  }
  
  TimersGet()->Start(0);
  
  // While there is active client
  while(!LastSent(mNoClients)){
    pthread_mutex_lock(mpReceivedMutex);
     size = mReceivedElements.size();
    pthread_mutex_unlock(mpReceivedMutex);
    // if( (mSync) ?(size >= ActiveClients()) :(size >= mNoClients || Finished(mNoClients))){ 
    // Changed for mNcumrbu
    if( (mSync) ?(size >= ActiveClients()) :(size >= mNcumrbu || Finished(mNoClients))){ 
      // if there is enough received elements to make update
      if(DEBUG_PROG) std::cerr << "Have enough " << size << "\n";
      element_add->Clear();
      for(int i=0; i < size; i++){
        pthread_mutex_lock(mpReceivedMutex);
         element = mReceivedElements.front();
         mReceivedElements.pop();
        pthread_mutex_unlock(mpReceivedMutex);
        element_add->Add(element, 1.0); // add update matrixes
        pthread_mutex_lock(mpFreeMutex);
         mFreeElements.push(element);
        pthread_mutex_unlock(mpFreeMutex);
      }
      // element_nn->Add(element_add, -1.0*mpNNet->LearnRate()); // make new weights
      element_nn->SubByLearningRate(element_add, mpNNet->LearningRate()); // make new weights
      
      // For all clients, send element if client is active
      for(int i=0; i<mNoClients; i++){
        if(!mpLastSent[i]){ // only if last element is not sent
          if(mpClientFinished[i]){
            element_nn->mLast = 1;
            mpLastSent[i] = true;
          }
          else{
            element_nn->mLast = 0;
          }
          mpServer->SendElement(element_nn, i); // send new weights
          if(DEBUG_PROG) {
            if(element_nn->mLast == 0) std::cerr << "New weights (client "<<i<<") sent\n";
            else std::cerr << "New weights (client "<<i<<") sent LAST\n";  
          }
        }
      }
      
    }  
  }
  
  TimersGet()->End(0);
  
  barrier_wait(mpBarrier); // waiting for all clients to send info
  int vectors = 0;
  int good = 0;
  int discarded = 0;
  for(int i=0; i < mNoClients; i++){
     vectors += mpRecVectors[i];
     good += mpRecGood[i];
     discarded += mpRecDiscarded[i];
  }
  mpNNet->Vectors(vectors);
  mpNNet->Good(good);
  mpNNet->Discarded(discarded);
  mpNNet->PrintInfo(); // print numbers of vectors
  
  delete[] mpRecVectors; 
  delete[] mpRecGood; 
  delete[] mpRecDiscarded; 
}

void SNet::ProgObj::ServerReceivingThread(int number){
  Element *element;
  while(!mpClientFinished[number]){
    element = GetOrCreate(&mFreeElements, mpNNet, mpFreeMutex);
    if(DEBUG_PROG) std::cerr << "THREAD " << number << " : waiting for element\n";
    mpServer->ReceiveElement(element, number);
    pthread_mutex_lock(mpReceivedMutex);
     mReceivedElements.push(element);
     if(DEBUG_PROG) {
       if(element->mLast == 1) std::cerr << "THREAD " << number << " : element received LAST\n";
       else std::cerr << "THREAD " << number << " : element received\n"; 
     }
    pthread_mutex_unlock(mpReceivedMutex);
    if(element->mLast == 1){ // last indicates that this client is done
      mpClientFinished[number] = true;
    }
  }
  if(DEBUG_PROG) std::cerr << "THREAD " << number << " waiting for numbers!\n";
  mpRecVectors[number] = mpServer->ReceiveInt(number);
  mpRecGood[number] = mpServer->ReceiveInt(number);
  mpRecDiscarded[number] = mpServer->ReceiveInt(number);
  if(DEBUG_PROG) std::cerr << "THREAD " << number << " receivec numbers!\n";
  barrier_wait(mpBarrier); // waiting for all clients to send info
}

/// ********** Client **********
void 
SNet::ProgObj::
RunClient()
{
  // Open client receiving thread
  mpThreads = new pthread_t[1]; // because of delete[]
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
  mpNNet->Mutexes(mpFreeMutex, mpReceivedMutex, mpBarrier, &mSync);
  
  mpNNet->WaitForStartingWeights(); // wait for first weights
  
  if(DEBUG_PROG) if(mSync) std::cerr << "THREAD: Waiting on barrier\n";
  if(mSync) barrier_wait(mpBarrier); // to destroy barrier in receiving thread
  
  TimersGet()->Start(0);
  // Do nothing, return to program
}


//******************************************************************************
//******************************************************************************
void 
SNet::ProgObj::
ClientReceivingThread()
{
  Element *element;
  mClientShouldFinish = false;

  while(!mClientShouldFinish)
  {
    element = GetOrCreate(&mFreeElements, mpNNet, mpFreeMutex);
    mpClient->ReceiveElement(element);
    pthread_mutex_lock(mpReceivedMutex);
    mReceivedElements.push(element);
    mpNNet->TimersGet()->Count(0);

    if(DEBUG_PROG) 
    {
      if(element->mLast == 1)  std::cerr << "THREAD: Received Weights " << mpNNet->TimersGet()->Counter(0) << " LAST\n";
      else std::cerr << "THREAD: Received Weights " << mpNNet->TimersGet()->Counter(0) << "\n";
    }

    pthread_mutex_unlock(mpReceivedMutex);
    if(DEBUG_PROG) if(mSync)  std::cerr << "THREAD: Waiting on barrier\n";
    if(mSync) barrier_wait(mpBarrier);
    if(element->mLast == 1){ // waiting for last element from server even if useless
      mClientShouldFinish = true;
    } 
  }

  if(DEBUG_PROG) std::cerr << "THREAD: END OF THREAD WAITING ON END-BARRIER\n";

  barrier_wait(mpEndBarrier);
}
