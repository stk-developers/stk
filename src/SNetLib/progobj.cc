#include"progobj.h"

using namespace SNet;
using std::cerr;
      
ProgObj::ProgObj():
  mServer(0),
  mClient(0),
  mCrossValidation(0),
  mPort(2020),
  mIp("IP not set"),
  mNoClients(0),
  mpThreads(NULL), 
  mVersion("2.0.0"),
  mHelp( 
    "There should be \n"
    "a very long help. \n"
  ), 
  mpMainSocket(NULL),
  mpclientSocket(NULL)
{}

void ProgObj::PrintCommandLine(int argc, char *argv[]){
  for (int i = 0; i<argc; i++){
    cerr << argv[i] << " ";
  }
  cerr << "\n";
}

void ProgObj::PrintVersion(){
  cerr << "Program version " << mVersion << "\n";
}

void ProgObj::PrintHelp(){
  cerr << mHelp;
}
