#ifndef PROGOBJ_H
#define PROGOBJ_H

#include"socketobj.h"

namespace SNet{

// Core SNet program class
class Progobj{
  private:
    int mServer;               ///< If program is server 1, else 0
    int mClient;               ///< If program is client 1, else 0
    int mCrossValidation;      ///< If program is server 1, else 0
    int mPort;                 ///< Used network port number
    char mIp[17];              ///< Client only - IP adress where client should join
    int mNoClients;            ///< Server only - Number of clients
    pthread_t* mpThreads;      ///< Pointer to threads

    Socketobj mMainSocket;     ///< Main socket
    Socketobj* mpclientSocket; ///< Server only - Pointer to sockets to clients
  public:

};

} // namespace

#endif
