#pragma once

#include <string>
#include <vector>
#include <set>

class cInfo;

struct sChangeLog {

    std::string version;
    std::string log;
    std::string date;
    
};

struct cParameters {

    // output
    cInfo   *info;

	// names of files, global variables
	std::string inputFileName;
    std::string outputFileName;

	int nCPUs;

    static    std::vector<sChangeLog> changeLog;

    static    void    setChangeLog();
    std::string    getChangeLog();

	void parseInputArguments(int argc, char** argv);

    cParameters(cInfo *i) {info = i;}

    ////////////////////////////////
    //Local machine
    ///////////////////////////////
    char hostname[128];
    std::string startingTime;


private:
    std::string arguments;
    std::string logo;

    void    init();
    void    setLogo();
    void    setArguments(int argc, const char * const * argv);

};
