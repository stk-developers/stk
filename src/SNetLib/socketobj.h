#ifndef SOCKETOBJ_H
#define SOCKETOBJ_H

#include <sys/socket.h>

namespace SNet{
  /// Object for sockets
  class SocketObj{
    private:
      int mSocket;                            ///< Socket identifier
    public:
      int ReceiveInt();                       ///< Receives one integer from socket
      void SendInt(int data);                 ///< Sends one integer to socket
      void ReceiveData(char* data, int n);    ///< Receives N bytes to array DATA using socket
      void SendData(char* data, int n);       ///< Sends N bytes from array DATA using socket
  };
} // namespace

#endif
