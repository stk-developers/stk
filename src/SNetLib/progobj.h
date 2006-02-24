#ifndef PROGOBJ_H
#define PROGOBJ_H

#include "nnet.h"
#include <iostream>
#include <string>
#include <queue>


namespace SNet{
  class ProgObj;

  typedef struct ThreadData {
    ProgObj* progObj;
    int number;
  } ThreadData;

  // Core SNet program class
  class ProgObj{
    private:
      bool mServer;    ///< If program is server true
      bool mClient;    ///< If program is client true
      int mPort;       ///< Used network port number
      std::string mIp; ///< Client only - IP adress where client should join
      int mNoClients;  ///< Server only - Number of clients
      bool mClientShouldFinish;
      
      std::queue<Element*> mFreeElements;
      std::queue<Element*> mReceivedElements;
      
      pthread_mutex_t *mpFreeMutex;
      pthread_mutex_t *mpReceivedMutex;
      barrier_t *mpBarrier;
      
      bool *mpClientFinished;
      
      pthread_t* mpThreads;      ///< Pointer to threads
      //SocketObj* mpMainSocket;   ///< Main socket
      //SocketObj* mpClientSocket; ///< Server only - Pointer to sockets to clients
      Socket::Server *mpServer;
      Socket::Client *mpClient;
      NNet *mpNNet;              ///< Main NN structure
      
      int GiveMeSeed(); ///< Generates random seed
        
      bool Finished(int noClients);
      Element *GetOrCreate(std::queue<Element*> *que, NNet *nn, pthread_mutex_t *mutex);
    public:
      ProgObj(XformInstance *NNetInstance, int cacheSize, int bunchSize, bool crossValidation, ///< Constructor
              std::string version, float learningRate, int clients, char* ip); 
      ~ProgObj();                                                                              ///< Destructor
      
      void NewVector(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize, bool last); ///< Work with new input vector
      void RunServer();
      void RunClient();
      void ServerReceivingThread(int number);
      void ClientReceivingThread();
      
      // Accessors
      bool Server() const {return mServer;};
      bool Client() const {return mClient;};
  };
} // namespace
#endif
