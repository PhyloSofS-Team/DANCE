#include "cTimer.hpp"
#include <stdio.h>


#if defined(_WIN64) || defined(_WIN32)
#define exp7           10000000i64     //1E+7     //C-file part
#define exp9         1000000000i64     //1E+9
#define w2ux 116444736000000000i64     //1.jan1601 to 1.jan1970

void unix_time(struct timespec *spec)
{
	__int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
	wintime -= w2ux;  spec->tv_sec = wintime / exp7;
	spec->tv_nsec = wintime % exp7 * 100;
}
int clock_gettime(int, timespec *spec)
{
	static  struct timespec startspec; static double ticks2nano;
	static __int64 startticks, tps = 0;    __int64 tmp, curticks;
	QueryPerformanceFrequency((LARGE_INTEGER*)&tmp); //some strange system can
	if (tps != tmp) {
		tps = tmp; //init ~~ONCE         //possibly change freq ?
		QueryPerformanceCounter((LARGE_INTEGER*)&startticks);
		unix_time(&startspec); ticks2nano = (double)exp9 / tps;
	}
	QueryPerformanceCounter((LARGE_INTEGER*)&curticks); curticks -= startticks;
	spec->tv_sec = startspec.tv_sec + (curticks / tps);
	spec->tv_nsec = startspec.tv_nsec + (double)(curticks % tps) * ticks2nano;
	if (!(spec->tv_nsec < exp9)) { spec->tv_sec++; spec->tv_nsec -= exp9; }
	return 0;
}
#endif // defined(_WIN64) || defined(_WIN32)


void	cTimer::printTimeInterval(const char* info) {

	double dt = this->getTimeInterval();
	printf("%s \t\t%g s\n", info, dt);

}

void	cTimer::printTime(const char* info) {

	double dt = this->getTime();
	printf("%s \t\t%g s\n", info, dt);

}

cSequentialTimer::cSequentialTimer() {

#if defined(__APPLE__) // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	baseSequential.tv_sec = mts.tv_sec;
	baseSequential.tv_nsec = mts.tv_nsec;
#elif defined(__linux__)
	clock_gettime(CLOCK_REALTIME, &baseSequential);
#elif defined(_WIN64) || defined(_WIN32)
	clock_gettime(0, &baseSequential);
#endif
	currentSequential = baseSequential;
}

void cSequentialTimer::beginMeasurement() {
#if defined(__APPLE__) // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	baseSequential.tv_sec = mts.tv_sec;
	baseSequential.tv_nsec = mts.tv_nsec;
#elif defined(__linux__)
	clock_gettime(CLOCK_REALTIME, &baseSequential);
#elif defined(_WIN64) || defined(_WIN32)
	clock_gettime(0, &baseSequential);
#endif
	currentSequential = baseSequential;
}

double cSequentialTimer::getTimeInterval() {

	struct timespec cur;
#ifdef __APPLE__ // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	cur.tv_sec = mts.tv_sec;
	cur.tv_nsec = mts.tv_nsec;
#elif defined(__linux__)
	clock_gettime(CLOCK_REALTIME, &cur);
#elif defined(_WIN64) || defined(_WIN32)
	clock_gettime(0, &cur);
#endif
	double dt = (cur.tv_sec - currentSequential.tv_sec)
		+ (double)(cur.tv_nsec - currentSequential.tv_nsec)
		/ (double)1000000000L;
	times.push_back(dt);

	currentSequential = cur;
	return dt;
}


double cSequentialTimer::getTime() {

	struct timespec cur;
#ifdef __APPLE__ // OS X does not have clock_gettime, use clock_get_time
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	cur.tv_sec = mts.tv_sec;
	cur.tv_nsec = mts.tv_nsec;
#elif defined(__linux__)
	clock_gettime(CLOCK_REALTIME, &cur);
#elif defined(_WIN64) || defined(_WIN32)
	clock_gettime(0, &cur);
#endif
	double dt = (cur.tv_sec - baseSequential.tv_sec)
		+ (double)(cur.tv_nsec - baseSequential.tv_nsec)
		/ (double)1000000000L;
	return dt;
}

#ifdef MULTICORE

cParallelTimer::cParallelTimer() {

	baseTime = omp_get_wtime();
	currentTime = baseTime;

}

double cParallelTimer::getTimeInterval() {

	double cur = omp_get_wtime();
	double dt = cur - currentTime;
	currentTime = cur;
	times.push_back(dt);
	return dt;
}

double cParallelTimer::getTime() {

	double cur = omp_get_wtime();
	return cur - baseTime;
}
#endif
