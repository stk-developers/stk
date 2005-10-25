#ifndef PROGOBJ_H
#define PROGOBJ_H

#include"socketobj.h"

class Progobj{
  private:
    int server;
    int client;
    int crossValidation;
    int port;
    int ip;
    int nOfClients;
    int *threads;

    Socketobj mainSocket;
    Socketobj *clientSocket;
  public:

};

#endif
