#ifndef timer_H
#define timer_H
#include <sys/time.h>
#include <sys/resource.h>

class timer {
private:
  struct rusage res;
  struct timeval tp;
  double virtual_time, real_time;

public:
  enum TYPE {REAL, VIRTUAL};
  timer(void);
  void resetTime();
  double elapsedTime(const TYPE& type);
};
#endif

