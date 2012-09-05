#include "utilities.h"

using namespace std;

static time_t time_begin;
static time_t time_last;

time_t Initial_Time()
{
	time_begin = time(NULL);
	time_last = time_begin;
	return time_begin;
};

//time used during the past step
time_t Cal_StepTime()
{
	time_t tused = time(NULL)-time_last;
	time_last = time(NULL);
	return tused;
};

//total time exhaust
time_t Cal_AllTime()
{
	return time(NULL)-time_begin;
};

//current time on string format
char * Curr_Time()
{
	time_t t=time(NULL);
	return ctime(&t);
}
