/**
 * @author      : Evan Gerritz (evan.gerritz@yale.edu)
 * @file        : timing
 * @created     : Saturday Dec 17, 2022 21:07:45 EST
 */

#include "timing.h"

void timing_(double* wcTime, double* cpuTime)
{
   timing(wcTime, cpuTime);
}

void timing(double* wcTime, double* cpuTime)
{
   struct timeval tp;
   struct rusage ruse;

   gettimeofday(&tp, NULL);
   *wcTime=(double) (tp.tv_sec*1000.0 + tp.tv_usec/1000.0); 
  
   if (cpuTime) {
       getrusage(RUSAGE_SELF, &ruse);
       *cpuTime=(double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0);
   }
}
