#ifndef SOCKETOBJ_H
#define SOCKETOBJ_H

// STK Common include
#include "element.h"
#include "barrier.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <unistd.h>


#define DEBUG 0

typedef struct sockaddr_in Tsockaddr_in;
typedef struct sockaddr Tsockaddr;
typedef struct in_addr Tin_addr;
typedef struct hostent Thostent;

namespace Socket{
  //! Object for sockets
  class SocketObj{
    protected:
      int mSocket;  ///< Socket identifier
    public:
      int ReceiveInt();                     ///< Receives one integer from socket
      void SendInt(int data);               ///< Sends one integer to socket
      void ReceiveData(char* data, int n);  ///< Receives N bytes to array DATA using socket
      void SendData(char* data, int n);     ///< Sends N bytes from array DATA using socket
  };
  
  //! Object for server
  class Server : public SocketObj{
    private:
      int *pmClientSocket; ///< Pointer of client sockets
      int mNOfClients;     ///< Number of clients
    public:
      Server(int port, int nOfClients); ///< Constructor for server
      ~Server();                        ///< Destructor for server
      
      void SendInt(int data, int who);
      int ReceiveInt(int who);
      void Send(char *data, int n, int who);
      void SendElement(SNet::Element *element, int who);
      void Receive(char *data, int n, int who);
      void ReceiveElement(SNet::Element *element, int who);
      void SendIntBroad(int data);
      void SendBroad(char *data, int n);
      void SendElementBroad(SNet::Element *element);
  }; 
  
  //! Object for client
  class Client : public SocketObj{
    private:
    public:
      Client(int port, char *ip); ///< Constructor for client
      ~Client();                  ///< Destructor for client
      
      void SendInt(int data);
      int ReceiveInt();
      void Send(char *data, int n);
      void SendElement(SNet::Element *element);
      void Receive(char *data, int n);
      void ReceiveElement(SNet::Element *element);
  }; 
} // namespace
#endif
