#ifndef _TIMER_H_ 
#define _TIMER_H_

/**
*    Class to get timing results (basic stopwatch).
*    Reports time in seconds, accurate up to microseconds
*	   
*/

#include <string>
#include "config.h"

#ifdef HAS_BOOST
#include "boost/timer/timer.hpp"

class  Timer
{
public:
  Timer();
  ~Timer();

  //initialize timer(set the elapsed time to 0)
  void Init();

  //start the timer(also set the elapsed time to 0)
  void Start();

  //Restart the tiemr(do not reset the elapsed time)
  void ReStart();

  //stop the timer 
  void Stop();

  //get the elapsed time in seconds or milliseconds
  double GetElapsedTimeSecond();
  double GetElapsedTimeMilliSecond();

  //get the elapsed time in seconds as a string
  std::string GetElapsedTimeAsString();

private:

  boost::timer::cpu_timer ter;

  boost::timer::nanosecond_type  mElapsed_time;//in u seconds;
};
#else

class  Timer
{
public:
  void Start();
  void Stop();
  void ReStart();
  double GetElapsedTimeSecond();
  double GetElapsedTimeMilliSecond();
  std::string GetElapsedTimeAsString();
private:
  double ter;
  double mElapsed_time;
};

#endif


#endif