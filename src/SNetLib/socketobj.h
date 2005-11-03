#ifndef SOCKETOBJ_H
#define SOCKETOBJ_H

#include <sys/socket.h>

namespace SNet{

  class Socketobj{
    private:
      int mSocket;                            ///< Socket identifier
    public:
      int receiveInt();                       ///< Receives one integer from socket
      void sendInt(int data);                 ///< Sends one integer to socket
      void receiveData(char* data, int n);    ///< Receives N bytes to array DATA using socket
      void sendData(char* data, int n);       ///< Sends N bytes from array DATA using socket
  };

} // namespace

#endif
