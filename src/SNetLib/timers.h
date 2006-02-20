#ifndef TIMERS_H
#define TIMERS_H

#include <assert.h>
#include <sys/time.h>
#include <stdlib.h>

namespace SNet{
  // Core SNet program class
  class Timers{
   private: 
     double *mpTimer;
     double *mpTimerAdd;
     int *mpCounter;  
     
     double GiveTime();  
   public:
     Timers(int nTimers, int nCounters);
     ~Timers();
     void Start(int number);
     void End(int number);
     double Timer(int number);
     void Count(int number);
     int Counter(int number);
   };
}   
#endif
