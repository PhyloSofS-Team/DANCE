#pragma once

#include <ctime>
#include <vector>

#if defined(_WIN64) || defined(_WIN32)
#include <windows.h>
#endif

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

class cTimer {
protected:
	struct timespec baseSequential, currentSequential;
	double baseTime, currentTime;
public:
	std::vector<double> times;
	virtual double getTimeInterval() = 0;
	virtual double getTime() = 0;
	void	printTimeInterval(const char*);
	void	printTime(const char*);

};

class cSequentialTimer : public cTimer {

public:
	cSequentialTimer();
	virtual void   beginMeasurement();
	virtual double getTimeInterval();
	virtual double getTime();

};

#ifdef MULTICORE

class cParallelTimer : public cTimer {

public:
	cParallelTimer();
	virtual double getTimeInterval();
	virtual double getTime();
};
#endif

//class cProgressBar {
//private:
//    int width;
//    const int barWidth = 70;
//    int displayedPos;
//public:
//
//    cProgressBar(int _w) {
//        width = _w;
//        displayedPos = -1;
//    };
//
//    ~ cProgressBar() {
//        std::cout << std::endl;
//    }
//    void print(int currentPos) {
//        float progress = float (currentPos+1)/width;
//        int pos = barWidth * progress;
//        if (pos == displayedPos) {
//            return;
//        } else {
//            displayedPos = pos;
//        }
//        std::cout << "[";
//        for (int i = 0; i < barWidth; ++i) {
//            if (i < pos) std::cout << "=";
//            else if (i == pos) std::cout << ">";
//            else std::cout << " ";
//        }
//        std::cout << "] " << int(progress * 100.0) << " %\r";
//        std::cout.flush();
//        
//    }
//};
