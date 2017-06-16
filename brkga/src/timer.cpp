#include "timer.h"

/*
 *  The virtual time of day and the real time of day are calculated and
 *  stored for future use.  The future use consists of subtracting these
 *  values from similar values obtained at a later time to allow the user
 *  to get the amount of time used by the algorithm.
 */
timer::timer(void) {
  getrusage( RUSAGE_SELF, &res );
  virtual_time = (double) res.ru_utime.tv_sec +
    (double) res.ru_stime.tv_sec +
    (double) res.ru_utime.tv_usec * 1.0E-6 +
    (double) res.ru_stime.tv_usec * 1.0E-6;
  
  gettimeofday( &tp, 0 );
  real_time =    (double) tp.tv_sec + (double) tp.tv_usec * 1.0E-6;
}

void
timer::resetTime() {
  getrusage( RUSAGE_SELF, &res );
  virtual_time = (double) res.ru_utime.tv_sec +
    (double) res.ru_stime.tv_sec +
    (double) res.ru_utime.tv_usec * 1.0E-6 +
    (double) res.ru_stime.tv_usec * 1.0E-6;
  
  gettimeofday( &tp, 0 );
  real_time =    (double) tp.tv_sec + (double) tp.tv_usec * 1.0E-6;
}

/*
 *  Stop the stopwatch and return the time used in seconds (either
 *  REAL or VIRTUAL time, depending on ``type'').
 */
double timer::elapsedTime(const TYPE& type) {
  if (type == REAL) {
    gettimeofday( &tp, 0 );
    return( (double) tp.tv_sec + (double) tp.tv_usec * 1.0E-6 - real_time );
  }
  else {
    getrusage( RUSAGE_SELF, &res );
    return( (double) res.ru_utime.tv_sec +
	    (double) res.ru_stime.tv_sec +
	    (double) res.ru_utime.tv_usec * 1.0E-6 +
	    (double) res.ru_stime.tv_usec * 1.0E-6
	    - virtual_time );
  }
}

