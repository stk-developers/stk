#include"timers.h"

double SNet::Timers::GiveTime(){
  struct timeval tv;
  if (gettimeofday(&tv, 0) == -1) {
    assert(0 && "gettimeofday does not work.");
    exit(-1);
  }
  return (double)(tv.tv_sec) + (double)tv.tv_usec/1000000.;
};

SNet::Timers::Timers(int nTimers, int nCounters){
  mpTimer = new double[nTimers];
  mpTimerAdd = new double[nTimers];
  mpCounter = new int[nCounters]; 
  for(int i=0; i < nCounters; i++)
    mpCounter[i] = 0;
  for(int i=0; i < nTimers; i++)
    mpTimerAdd[i] = 0;
}

SNet::Timers::~Timers(){
  delete[] mpTimer;
  delete[] mpTimerAdd;
  delete[] mpCounter; 
}

void SNet::Timers::Start(int number){
  mpTimer[number] = GiveTime();
}

void SNet::Timers::End(int number){
  mpTimerAdd[number] += GiveTime() - mpTimer[number];
}

double SNet::Timers::Timer(int number){
  return mpTimerAdd[number];
}

void SNet::Timers::Count(int number){
  mpCounter[number]++;
}

int SNet::Timers::Counter(int number){
  return mpCounter[number];
}
