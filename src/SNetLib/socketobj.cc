#include"socketobj.h"
#include "../STKLib/Error.h"

using namespace Socket;

Server::Server(int port, int nOfClients){
  mNOfClients = nOfClients;
  int socket_num;
  if((socket_num = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == -1){
    STK::Error("Cannot create socket");
  }
  Tsockaddr_in sock_name;
  sock_name.sin_family = AF_INET;
  sock_name.sin_port = htons(port);
  sock_name.sin_addr.s_addr = INADDR_ANY;
  int opt = 1;
  setsockopt(socket_num, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof (opt));
  if(bind(socket_num, (Tsockaddr*)&sock_name, sizeof(sock_name)) == -1){
    STK::Error("Cannot name socket");
  }
  printf("SERVER is running\n");
  mSocket = socket_num;
  
  if(listen(mSocket, 10) == -1){
    STK::Error("Cannot make listen queue");
  }
  Tsockaddr_in clientInfo;
  socklen_t addrlen = sizeof(clientInfo);
  int i;
  
  pmClientSocket = new int[nOfClients];
  
  // Create all connections
  for(i=0; i<mNOfClients; i++){ 
    int client = accept(mSocket, (Tsockaddr*)&clientInfo, &addrlen);
    if(client == -1){
      STK::Error("Problems with client join");
    }
    printf("SERVER accepted client from %s\n",inet_ntoa((Tin_addr)clientInfo.sin_addr));
    pmClientSocket[i] = client;
  }
  
}

Server::~Server(){
  for(int i=0; i<mNOfClients; i++){
    shutdown(pmClientSocket[i], SHUT_RDWR);
  }
  delete [] pmClientSocket;
  shutdown(mSocket, SHUT_RDWR);
}

void Server::SendInt(int data, int who){
  int size;
  char *pom = (char*)&data;
  if(DEBUG) printf("Sending int to client %d (%d)\n", who, data);
  size = send(pmClientSocket[who], pom, sizeof(int), 0);
  if(DEBUG) printf("Sent %d bytes to client\n", size);
  if(size == -1){
    STK::Error("Server can not send int to client %d - problem with sending data", who);
  }
}

int Server::ReceiveInt(int who){
  int data, size;
  int totalSize = 0;
  char *pom = (char*)&data;
  if(DEBUG) printf("Waiting for int from client %d\n", who);
  while(totalSize != sizeof(int)){
    size = recv(pmClientSocket[who], pom, sizeof(int) - totalSize, 0);
    if(size == -1){
      STK::Error("Server can not receive int from client %d - problem with receiving data", who);
    }
    pom += size;
    totalSize += size;
  }
  if(DEBUG) printf("Received from client %d (%d bytes)\n", who, totalSize);
  if(DEBUG) printf("Received int from client %d (%d)\n", who, data);
  return data;
}

void Server::Send(char *data, int n, int who){
  int size;
  char *pom = data;
  if(DEBUG) printf("Sending data to client %d (%d bytes)\n", who, n);
  size = send(pmClientSocket[who], pom, n, 0);
  if(DEBUG) printf("Sent data to client %d (%d bytes)\n", who, size);
  if(size == -1){
    STK::Error("Server can not send data to client %d - problem with sending data", who);
  }
}

void Server::SendElement(SNet::Element *element, int who){
   assert(element != NULL);
   SendInt(element->mLast, who);
   SendInt(element->mNLayers, who);
   for(int i=0; i < element->mNLayers; i++){
     Send((char*)(*(element->mpWeights[i]))[0], element->mpWeights[i]->MSize(), who);
     Send((char*)(*(element->mpBiases[i]))[0], element->mpBiases[i]->MSize(), who);     
   }
}

void Server::Receive(char *data, int n, int who){
  int totalSize, size;
  char *pom = data;
  totalSize = 0;
  if(DEBUG) printf("Receiveing data from client %d (wait for %d bytes)\n", who, n);
  while(totalSize != n){
    size = recv(pmClientSocket[who], pom, n - totalSize, 0);
    if(size == -1){
      STK::Error("Server can not receive data from client %d - problem with receiving data", who);
    }
    pom += size;
    totalSize += size;
  }
  if(DEBUG) printf("Received data from client %d (%d bytes)\n", who, totalSize);
}

void Server::ReceiveElement(SNet::Element *element, int who){
   element->mLast = ReceiveInt(who);
   element->mNLayers = ReceiveInt(who);   
   for(int i=0; i < element->mNLayers; i++){
     Receive((char*)(*(element->mpWeights[i]))[0], element->mpWeights[i]->MSize(), who);
     Receive((char*)(*(element->mpBiases[i]))[0], element->mpBiases[i]->MSize(), who);     
   }
}

void Server::SendIntBroad(int data){
  for(int i=0; i < mNOfClients; i++){
    SendInt(data, i);
  }
}

void Server::SendBroad(char *data, int n){
  for(int i=0; i < mNOfClients; i++){
    Send(data, n, i);
  }
}

void Server::SendElementBroad(SNet::Element *element){
  for(int i=0; i < mNOfClients; i++){
    SendElement(element, i);
  }
}

Client::Client(int port, char *ip){
  int socket_num;
  Thostent *host;
  if((socket_num = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == -1){
    STK::Error("Cannot create socket");
  }
  if((host = gethostbyname(ip)) == NULL){
    STK::Error("Cannot resolve host name");
  }
  Tsockaddr_in sock_name;
  sock_name.sin_family = AF_INET;
  sock_name.sin_port = htons(port);
  memcpy(&(sock_name.sin_addr), host->h_addr, host->h_length);
  while(connect(socket_num, (Tsockaddr*)&sock_name, sizeof(sock_name)) == -1){
    fprintf(stderr, "WAITING - cannot connect to server\n");
    sleep(1);
  }
  printf("CLIENT joined server %s\n", ip);
  mSocket = socket_num;
}

Client::~Client(){
  shutdown(mSocket, SHUT_RDWR);
}

void Client::SendInt(int data){
  int size;
  char *pom = (char *)&data;
  if(DEBUG) printf("Sending int to server (%d)\n", data);
  size = send(mSocket, pom, sizeof(int), 0);
  if(DEBUG) printf("Sent %d bytes to server\n", size);
  if(size == -1){
    STK::Error("Client can not send int to server - problem with sending data");
  }
}

int Client::ReceiveInt(){
  int data, size;
  int totalSize = 0;
  char *pom = (char*)&data;
  if(DEBUG) printf("Waiting for int from server\n");
  while(totalSize != sizeof(int)){
    size = recv(mSocket, pom, sizeof(int) - totalSize, 0);
    if(size == -1){
      STK::Error("Client can not receive int from server %d - problem with receiving data");
    }
    pom += size;
    totalSize += size;
  }
  if(DEBUG) printf("Received from server (%d bytes)\n", totalSize);
  if(DEBUG) printf("Received int from server (%d)\n", data);
  return data;
}

void Client::Send(char *data, int n){
  int size;
  char *pom = data;
  if(DEBUG) printf("Sending data to server (%d bytes)\n", n);
  size = send(mSocket, pom, n, 0);
  if(DEBUG) printf("Sent data to server (%d bytes)\n", size);
  if(size == -1){
    STK::Error("Server can not send data to server - problem with sending data");
  }
}

void Client::SendElement(SNet::Element *element){
   assert(element != NULL);
   SendInt(element->mLast);
   SendInt(element->mNLayers);
   for(int i=0; i < element->mNLayers; i++){
     Send((char*)(*(element->mpWeights[i]))[0], element->mpWeights[i]->MSize());
     Send((char*)(*(element->mpBiases[i]))[0], element->mpBiases[i]->MSize());     
   }
}

void Client::Receive(char *data, int n){
  int totalSize, size;
  char *pom = data;
  totalSize = 0;
  if(DEBUG) printf("Receiveing data from server (wait for %d bytes)\n", n);
  while(totalSize != n){
    size = recv(mSocket, pom, n - totalSize, 0);
    if(size == -1){
      STK::Error("Client can not receive data from server - problem with receiving data");
    }
    pom += size;
    totalSize += size;
  }
  if(DEBUG) printf("Received data from server (%d bytes)\n", totalSize);
}

void Client::ReceiveElement(SNet::Element *element){
   assert(element != NULL);
   element->mLast = ReceiveInt();
   element->mNLayers = ReceiveInt();
   for(int i=0; i < element->mNLayers; i++){
     Receive((char*)(*(element->mpWeights[i]))[0], element->mpWeights[i]->MSize());
     Receive((char*)(*(element->mpBiases[i]))[0], element->mpBiases[i]->MSize());     
   }
}
