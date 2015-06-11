// Really simple timer -- can be made more complicated if desired!

#ifndef _SimpleTimer_h_
#define _SimpleTimer_h_

typedef long long int64;

class Timer {

 public:

  Timer();
  ~Timer();

  void Start();
  void Stop();
  void Reset();

  double GetElapsedTime();

 private:

  int64 fStart, fCycles;

};

#endif
