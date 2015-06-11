#include <Timer.h>

#include <sys/time.h>

#include <time.h>

Timer::Timer() : fStart(0), fCycles(0) {

}

Timer::~Timer() {

}

void Timer::Start() {

  struct timeval t;
  gettimeofday(&t, 0);
  fStart = (int64)t.tv_sec*1000 + (int64)t.tv_usec/1000;
 
}

void Timer::Stop() {
  
  struct timeval t;
  gettimeofday(&t, 0);
  
  int64 end = (int64)t.tv_sec*1000 + (int64)t.tv_usec/1000;

  end -= fStart;
  fStart = 0;
  fCycles += end;

}

void Timer::Reset() {

  fCycles = 0;

}

double Timer::GetElapsedTime() {

  return double(fCycles)/1e3;

}
