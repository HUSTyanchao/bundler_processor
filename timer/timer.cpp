
#include <sstream>
#include <iostream>

#include "timer.h"

#ifdef HAS_BOOST

using namespace boost;

Timer::Timer():mElapsed_time(0){
	
}

Timer::~Timer(){}

//start the timer(also set the elapsed time to 0)
void Timer::Start(){
	mElapsed_time = 0;
}

//Restart the timer(reset the elapsed time)
void Timer::ReStart(){
  mElapsed_time = ter.elapsed().wall;
}

//stop the timer 
void Timer::Stop(){
	//since when restart we also add the original time
  mElapsed_time = ter.elapsed( ).wall - mElapsed_time;
}

//get the elapsed time in seconds 
double Timer::GetElapsedTimeSecond(){
	return (mElapsed_time / 1e9);
}
double Timer::GetElapsedTimeMilliSecond(){
  return (mElapsed_time / 1e6);
}

//get the elapsed time in seconds in a string
std::string Timer::GetElapsedTimeAsString(){
	// parse elapsed time into string stream
	std::ostringstream s;

	double elapsed_time_min = (mElapsed_time / 1e9);
	// get the elapsed minutes and seconds
	double elapsed_minutes = (elapsed_time_min < 60) ? 0.0 : floor(elapsed_time_min / 60);
	double elapsed_seconds = elapsed_time_min - floor(elapsed_time_min / 60) * 60;
	s << elapsed_minutes << " minutes " << elapsed_seconds << " seconds ";
	return s.str();
}
#else

void Timer::Start(){}
void Timer::Stop( ){}
void Timer::ReStart( ){}
double Timer::GetElapsedTimeSecond( ){ return  0; }
double Timer::GetElapsedTimeMilliSecond( ){ return 0; }
std::string Timer::GetElapsedTimeAsString(){ return std::string( "no timer" ); }


#endif // HAS_BOOST