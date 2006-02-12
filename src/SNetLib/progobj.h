#ifndef PROGOBJ_H
#define PROGOBJ_H

#include"socketobj.h"
#include"nnet.h"
#include<iostream>
#include<string>


namespace SNet{
  // Core SNet program class
  class ProgObj{
    private:
      int mServer;               ///< If program is server 1, else 0
      int mClient;               ///< If program is client 1, else 0
      int mPort;                 ///< Used network port number
      std::string mIp;           ///< Client only - IP adress where client should join
      int mNoClients;            ///< Server only - Number of clients
      pthread_t* mpThreads;      ///< Pointer to threads

      SocketObj* mpMainSocket;   ///< Main socket
      SocketObj* mpClientSocket; ///< Server only - Pointer to sockets to clients
      
      NNet *mpNNet;
    public:
      ProgObj(XformInstance *NNet_instance, int cache_size, int bunch_size, bool cross_validation, std::string version, float learning_rate);
      ~ProgObj();
      void NewVector(FLOAT *inVector, FLOAT *outVector, int inSize, int outSize, bool last);
  };
} // namespace

#endif
