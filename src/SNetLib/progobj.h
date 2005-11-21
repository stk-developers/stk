#ifndef PROGOBJ_H
#define PROGOBJ_H

#include"socketobj.h"
#include<iostream>
#include<string>


namespace SNet{
  // Core SNet program class
  class ProgObj{
    private:
      int mServer;               ///< If program is server 1, else 0
      int mClient;               ///< If program is client 1, else 0
      int mCrossValidation;      ///< If program is server 1, else 0
      int mPort;                 ///< Used network port number
      std::string mIp;           ///< Client only - IP adress where client should join
      int mNoClients;            ///< Server only - Number of clients
      pthread_t* mpThreads;      ///< Pointer to threads
      std::string mVersion;
      std::string mHelp;

      SocketObj* mpMainSocket;   ///< Main socket
      SocketObj* mpclientSocket; ///< Server only - Pointer to sockets to clients
    public:
      ProgObj();
      void PrintCommandLine(int argc, char *argv[]);
      void PrintVersion();
      void PrintHelp();
  };
} // namespace

#endif
