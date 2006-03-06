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

  //! Core SNet program class
  class ProgObj{
    private:
      bool mServer;                           ///< True if program is server
      bool mClient;                           ///< True if program is client
      bool mRandomize;                        ///< True if randomization switched on
      int mPort;                              ///< Used network port number
      std::string mIp;                        ///< Client only - IP adress where client should join
      int mNoClients;                         ///< Server only - Number of clients
      bool mClientShouldFinish;               ///< Client only - 
      bool mSync;                             ///< Client only - True if parallel synchronzation, false for asynchronous parallelization
      std::queue<Element*> mFreeElements;     ///< Queue of free elements
      std::queue<Element*> mReceivedElements; ///< Queue of received elements
      
      NNet *mpNNet;                     ///< Main NN structure
      pthread_mutex_t *mpFreeMutex;     ///< Mutex for queue of free elements
      pthread_mutex_t *mpReceivedMutex; ///< Mutex for queue of received elements
      barrier_t *mpBarrier;             ///< Barrier for synchronized parallelism
      barrier_t *mpEndBarrier;          ///< End of program (client) barrier
      bool *mpClientFinished;           ///< Array of flags if client finished
      bool *mpLastSent;                 ///< Array of flags if sent last data to client - really finished client
      pthread_t* mpThreads;             ///< Pointer to threads
      Socket::Server *mpServer;         ///< Server only - TCP/IP Server object
      Socket::Client *mpClient;         ///< Client only - TCP/IP Client object
      int *mpRecVectors;                ///< Array for received number of vectors from clients
      int *mpRecGood;                   ///< Array for received number of good vectors from clients 
      int *mpRecDiscarded;              ///< Array for received number of discarded vectors from clients 
      
      int GiveMeSeed();                                       ///< Generates random seed
      bool Finished(int noClients);                           ///< Returns if noClients clients finished
      bool LastSent(int noClients);                           ///< Returns if noClients clients really finished from server view
      int ActiveClients();                                    ///< Returns number of finished clients
      Element *GetOrCreate(std::queue<Element*> *que, 
                           NNet *nn, pthread_mutex_t *mutex); ///< Gets element from queue (using mutex) or makes new one (using NN)
    public:
      ProgObj(XformInstance *NNetInstance, int cacheSize, int bunchSize, bool crossValidation, ///< Constructor
              std::string version, float learningRate, int clients, char* ip, bool randomize, 
	      bool sync); 
      ~ProgObj();                                                                              ///< Destructor
      
      void NewVector(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize, bool last);   ///< Work with new input vector
      void RunServer();                                                                        ///< Runs server threads 
      void RunClient();                                                                        ///< Runs client threads 
      void ServerReceivingThread(int number);                                                  ///< Receiving thread for server
      void ClientReceivingThread();                                                            ///< Receiving thread for client 
      
      // Accessors
      bool Server() const {return mServer;};
      bool Client() const {return mClient;};
  };
} // namespace
#endif
