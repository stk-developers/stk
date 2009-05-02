
/***************************************************************************
 *   copyright            : (C) 2004 by Lukas Burget,UPGM,FIT,VUT,Brno     *
 *   email                : burget@fit.vutbr.cz                            *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#define SVN_DATE       "$Date: 2009-02-20 15:02:14 +0100 (Fri, 20 Feb 2009) $"
#define SVN_AUTHOR     "$Author: burget $"
#define SVN_REVISION   "$Revision: 272 $"
#define SVN_ID         "$Id: SERest.cc 272 2009-02-20 14:02:14Z burget $"

#define MODULE_VERSION "2.0.7 "__TIME__" "__DATE__" "SVN_ID  


#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
//#include "STKLib/Decoder.h"
//#include "STKLib/labels.h"
#include "STKLib/stkstream.h"
//#include "STKLib/Features.h"
//#include "STKLib/MlfStream.h"

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <malloc.h>
#include <assert.h>

// Use internal version of getopt function if not available
#ifndef HAVE_UNISTD_H
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include <iostream>
#include <string>
#include <sstream>
#include <queue>


#ifdef WIN32
#include <winsock2.h>
#include <process.h>
#include <conio.h>
#define socklen_t int
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
#define SOCKET int
#define WORD int
#define WSADATA int
#define MAKEWORD(a, b) 0
#define WSAStartup(a, b) 0
pthread_t thread;
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define SOCKET_ERROR   (-1)
#define INVALID_SOCKET (-1)
#define CoInitialize(a)
#define HOSTENT struct hostent
#define closesocket close
#define WSAGetLastError() errno
#define WSAENOTSOCK  ENOTSOCK
#endif


using namespace std;
using namespace STK;

//******************************************************************************
//******************************************************************************
void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\n%s version " MODULE_VERSION "\n"
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                                     Default\n\n"
" -d s       Dir to find hmm definitions                     Current\n"
" -o s       Extension for new hmm files                     As src\n"
" -s s       Print statistics to file s                      Off\n"
" -x s       Extension for hmm files                         None\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
" -H mmf     Load HMM macro file mmf                         \n"
" -M dir     Dir to write HMM macro files                    Current\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname, progname);
  exit(-1);
}

#ifdef WIN32
void process_sock(void *);
void accept_clients(void *);
void process_socktEx(void *params)
{
  try {
    process_sock(params);
  } catch(...){
  }
}
void accept_clientsEx(void *params)
{
  try {
    accept_clients(params);
  } catch(...){
  }
}
#else
void *process_sock(void *);
void *accept_clients(void *);
#endif


#define SNAME "STCPACCDUMP"
const char *optionStr =
" -d r   SOURCEMODELDIR"
" -x r   SOURCEMODELEXT"
" -D n   PRINTCONFIG=TRUE"
" -H l   SOURCEMMF"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE";


//" -o r   TARGETMODELEXT"
//" -s r   SAVESTATS"
//" -M r   TARGETMODELDIR"

SOCKET sock; // Listen sockets
pthread_cond_t cond   = PTHREAD_COND_INITIALIZER;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


template <typename t_data>
class volatile_queue {
  struct queue_elem {
    t_data data;
    queue_elem *next;
  };
  
  queue_elem *first;
  queue_elem *last;

public:
  volatile_queue() : first(NULL) {}
  
  void push(t_data &d)  volatile
  {
    queue_elem *tmp = new queue_elem;
    tmp->data = d;
    tmp->next = NULL;
    
    if(first == NULL) {
      first = tmp;
    } else {
      last->next = tmp;
    }
    last = tmp;
  }
  
  void pop()  volatile
  {
    assert(first != NULL);
    queue_elem *tmp = first;
    first = tmp->next;
    delete tmp;
  }

  t_data &front()  volatile
  { 
    assert(first != NULL);
    return first->data; 
  }
  
  bool empty() volatile
  { 
    return first == NULL; 
  }
};

volatile volatile_queue<pair<SOCKET, int> > sock_queue;
volatile volatile_queue<FILE *> file_queue;
char const *error_message = NULL;

//******************************************************************************
//******************************************************************************
int main(int argc, char* argv[]) 
{
  try
  {
    ModelSet                        hset;
    int                             i;
    int                             fcnt = 0;
    MyHSearchData                   cfgHash;
    FLOAT                           totLogLike      = 0;
    int                             totFrames       = 0;
    const char*                     src_hmm_list;
    const char*                     src_hmm_dir;
    const char*                     src_hmm_ext;
          char*                     src_mmf;
//    const char*                     trg_hmm_dir;
//    const char*                     trg_hmm_ext;
//    const char*                     trg_mmf;
    const char*                     xformList;
    const char*                     when_ready;
//    const char*                     stat_file;
//    const char*                     mix_occup_file;
    int                             trace_flag;
    int                             tcp_port;
    std::list<std::string> read_acc_names;

    WORD wVersionRequested;
    WSADATA wsaData;
    int err;
    struct sockaddr_in addr;

    // show help if no arguments specified
    if (argc == 1) {
      usage(argv[0]);
    }
    
    if (!my_hcreate_r(100,  &cfgHash)) {
      Error("Insufficient memory");
    }

        
    i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);

    gpFilterWldcrd=GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD","$");
    xformList    = GetParamStr(&cfgHash, SNAME":XFORMLIST",       NULL);
    gpHListFilter= GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
    gpMmfFilter  = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
    trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
//    stat_file    = GetParamStr(&cfgHash, SNAME":SAVESTATS",       NULL);
//    mix_occup_file=GetParamStr(&cfgHash, SNAME":SAVEMIXOCCUPS",   NULL);
    src_hmm_list = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
    src_hmm_dir  = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
    src_hmm_ext  = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
    src_mmf=(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
//    trg_hmm_dir  = GetParamStr(&cfgHash, SNAME":TARGETMODELDIR",  NULL);
//    trg_hmm_ext  = GetParamStr(&cfgHash, SNAME":TARGETMODELEXT",  NULL);
//    trg_mmf      = GetParamStr(&cfgHash, SNAME":TARGETMMF",       NULL);
    tcp_port     = GetParamInt(&cfgHash, SNAME":TCPPORT",         0);
    when_ready   = GetParamStr(&cfgHash, SNAME":WHENREADY",       NULL);
    
    // initialize basic ModelSet
    hset.Init(MODEL_SET_WITH_ACCUM);

    bool print_all_options = GetParamBool(&cfgHash,SNAME":PRINTALLOPTIONS", false);
    
    if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", false)) {
      PrintConfig(&cfgHash);
    }
    if (GetParamBool(&cfgHash, SNAME":PRINTVERSION", false)) {
      puts("Version: "MODULE_VERSION"\n");
    }
    if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", false)) {
      CheckCommandLineParamUse(&cfgHash);
    }
    if (print_all_options) {
      print_registered_parameters();
    } 
    
    if (NULL != src_mmf) {
      for (src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) {
        bool read_only = false;

        if (src_mmf[0] == '@') {
          src_mmf++;
        }
        std::cout << "Adding " << src_mmf << std::endl;
        hset.ParseMmf(src_mmf, NULL, read_only);
      }
    }

    if (src_hmm_list) {
      hset.ReadHMMList(src_hmm_list, src_hmm_dir ? src_hmm_dir : "", 
                                     src_hmm_ext ? src_hmm_ext : "");
    }
    
    hset.ResetAccums();
    
    if (xformList != NULL) {
      hset.ReadXformList(xformList);
    }
    hset.AllocateAccumulatorsForXformStats();
    
    wVersionRequested = MAKEWORD(2, 2);
 
    err = WSAStartup(wVersionRequested, &wsaData);
    if(err != 0) {
      Error("WSAStartup() failed");
    }
    if((sock = socket(AF_INET, SOCK_STREAM, 0)) == INVALID_SOCKET) {
      err = WSAGetLastError();
      Error("socket() failed %d", err);
    }
    memset(&addr, 0, sizeof(addr));
    addr.sin_family      = AF_INET;
    addr.sin_port        = htons(tcp_port);
    addr.sin_addr.s_addr = INADDR_ANY;
    
    if(bind(sock, (struct sockaddr *) &addr, sizeof(addr)) == SOCKET_ERROR) {
      err = WSAGetLastError();
      Error("bind() failed %d", err);
    }
    if((err = listen(sock, SOMAXCONN)) == SOCKET_ERROR ) {
      err = WSAGetLastError();
      Error("listen() failed %d", err);
    }
    
    socklen_t namelen = sizeof(addr);
    if((err = getsockname(sock, (struct sockaddr *) &addr, &namelen)) == SOCKET_ERROR ) {
      err = WSAGetLastError();
      Error("getsockname() failed %d", err);
    }
    tcp_port = ntohs(addr.sin_port);

    if(when_ready) {  
//      char host_name[1024] = {0};
//      std::ostringstream host_and_port;    
//      if(gethostname(host_name, sizeof(host_name)-1) != 0) {
//        err = WSAGetLastError();
//        Error("gethostkname() failed %d", err);
//      }
//      host_and_port << host_name << '/' << tcp_port;
      std::ostringstream ports;
      ports << tcp_port;
      char *cmd = expandFilterCommand(when_ready, ports.str().c_str());
      if(system(cmd) != 0) {
        Error("system() failed %d");
      }
      free(cmd);
    }

#ifdef WIN32
    _beginthread(process_sockEx, 0, NULL);
    _beginthread(accept_clientsEx, 0, NULL);
#else
    pthread_create(&thread, NULL,  process_sock, NULL);
    pthread_create(&thread, NULL,  accept_clients, NULL);
#endif
    TraceLog("Listening on port: %d", tcp_port);

    for(;;) {
      FILE *fp;
      char line_buff[1024];
      int  buff_len;
      
      char *chptr;
      FLOAT P;
      FLOAT P2;
      long  S;

      pthread_mutex_lock(&mutex);

      while(file_queue.empty() && NULL == error_message)
        pthread_cond_wait(&cond, &mutex);

      if(error_message)
        Error(error_message);
      
      fp = file_queue.front();
      file_queue.pop();
      pthread_mutex_unlock(&mutex);

      if(NULL == fp
      || NULL == fgets(line_buff, sizeof(line_buff), fp)) {
        pthread_mutex_lock(&mutex);
        if(error_message) {
          Error(error_message);
        } else {
          Error("Cannt read the input TCP stream");
        } 
      }
      if(line_buff[0] == 'q') {
        break;
      }
      if(line_buff[0] != 'r' && line_buff[0] != 'w') {
          Error("Invalid command in the TCP stream. First line must start with 'r', 'w' or 'q'.");
      }
      buff_len = strlen(line_buff);
      if(line_buff[buff_len-1] != '\n') {
        Error("Invalid command in the TCP stream. First line is too long");
      }
      line_buff[buff_len-1] = '\0';
      
      if(line_buff[0] == 'r') {
        read_acc_names.push_back(string(line_buff+1));
        
        if(trace_flag & 1) {
          TraceLog("Processing accumulator %d '%s'", ++fcnt, line_buff+1);
        }
        hset.ReadAccums(fp, line_buff+1, 1.0, &S, &P, UT_ML);
        totFrames  += S;
        totLogLike += P;

        if (trace_flag & 1) {
          TraceLog("[%d frames] %f", S, P/S);       
        }
      } else {
        if (trace_flag & 1) {
          TraceLog("Writing accumulator: '%s'", line_buff+1);
        }
        hset.WriteAccums(line_buff+1, NULL, totFrames, totLogLike);
        
        std::string acc_list_file = string(line_buff+1) + ".lst";
        OStkStream  output_stream(acc_list_file.c_str());
        
        if (!output_stream.good()) { 
          Error("Cannot open output accumulator list file: '%s'", acc_list_file.c_str());
        }
        for(std::list<std::string>::iterator it = read_acc_names.begin(); it != read_acc_names.end(); it++) {
          output_stream << *it << std::endl;
        }
        output_stream.close();
        
        if (trace_flag & 1) {
          TraceLog("Total number of frames: %d\nTotal log likelihood: %e", totFrames, totLogLike);
          TraceLog("Average log likelihood per frame: %e", totLogLike/totFrames);
        }
      }
      fclose(fp);

    }
/*    
    if (stat_file)
        hset.WriteHMMStats(stat_file);

    if (mix_occup_file)
      hset.WriteMixtureOccups(mix_occup_file);
*/
    closesocket(sock);  
    hset.Release();
    
    for (size_t i = 0; i < cfgHash.mNEntries; i++) 
      free(cfgHash.mpEntry[i]->data);
          
    my_hdestroy_r(&cfgHash, 1);
  } // try
  catch (std::exception& rExc) {
    std::cerr << "Exception thrown" << std::endl;
    std::cerr << rExc.what() << std::endl;
    return 1;
  }
  return 0;
}


#ifdef WIN32
void process_sock(void *params)
#else
void *process_sock(void *params)
#endif
{
  SOCKET sck;
  int fd;
  char buffer[256*1024];
  int buf_len = sizeof(buffer);
  
  for(;;) {
    pthread_mutex_lock(&mutex);

    while(sock_queue.empty())
      pthread_cond_wait(&cond, &mutex);
      
    sck = sock_queue.front().first;
    fd  = sock_queue.front().second;
    sock_queue.pop();
    pthread_mutex_unlock(&mutex);
    
//    socklen_t optlen = sizeof(buf_len); 
//    if(getsockopt(sock, SOL_SOCKET, SO_RCVBUF,(char *)&buf_len, &optlen) != 0)
//      buf_len = 4096;

    for(;;) {
      ssize_t i = recv(sck, buffer, buf_len, 0);
      
      if(0 == i)
        break;

      if(SOCKET_ERROR == i  
      || i != write(fd,  buffer, i)) 
      {
//        pthread_mutex_lock(&mutex);
//        error_message = "recv or write failed\n");
        closesocket(sck);
        close(fd);
        closesocket(sock);
//        pthread_cond_broadcast(&cond);
//        pthread_mutex_unlock(&mutex);
        return NULL;
      }
    }
    closesocket(sck);
    close(fd);
  }
}

#ifdef WIN32
void accept_clients(void *params)
#else
void *accept_clients(void *params)
#endif
{
  for(;;) {
    int filedes[2];
    FILE *fp = NULL;

    SOCKET socka = accept(sock, NULL, 0 );
//    printf("\n######### accept() %d #########\n", socka);
        
    if(socka == INVALID_SOCKET ) {  
      if(WSAGetLastError() == WSAENOTSOCK) {
        error_message = "recv or write failed";
      } else {
        error_message = "accept failed";
      }
    } else if(0 != pipe(filedes)) {
      error_message = "pipe failed";
    } else if(NULL == (fp = fdopen(filedes[0], "rb"))) {
      error_message = "fdopen failed";
    }
    pthread_mutex_lock(&mutex);

    if(NULL == error_message) {
      pair<SOCKET, int> tp(socka, filedes[1]);
      sock_queue.push(tp);
      file_queue.push(fp);
    }
    pthread_cond_broadcast(&cond);
    pthread_mutex_unlock(&mutex);

    if(NULL != error_message) {
      return NULL;
    } 
  }
}
