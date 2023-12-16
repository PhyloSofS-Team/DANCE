#include "cParameters.hpp"
#include "tclap/CmdLine.h"
#include "cInfo.hpp"
#include <locale>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif // _WIN32

std::vector<sChangeLog> cParameters::changeLog;

#include "tclap/Visitor.h"

class LogoVisitor : public TCLAP::Visitor {
protected:
    std::string _logo;
public:
    LogoVisitor(const std::string& name ) : Visitor(), _logo(name) {} ;
    void visit() { std::cout << _logo << std::endl;  exit(0); };
};

void cParameters::setChangeLog() {
    sChangeLog tmp;

    tmp.version = "0.1";
    tmp.date = "Nov 2021";
    tmp.log = "Initial release.";
    changeLog.push_back(tmp);
}

void cParameters::parseInputArguments(int argc, char** argv) {

    info->printInfo("Parsing the Command Line");

    cParameters::setChangeLog(); // calling static function
    cParameters::setLogo(); // calling static function

	// parse the arguments
	try {
        TCLAP::CmdLine cmd("", ' ', cParameters::changeLog.back().version+", "+cParameters::changeLog.back().date, true, getChangeLog());
		
		TCLAP::UnlabeledValueArg<std::string> _input("input","PDB input file for the complex",true,"","input PDB");
		cmd.add( _input );

		
		TCLAP::ValueArg<std::string> _out("o","output","Output path (by default, standard output will be used)",false,"","output file");
		cmd.add( _out );


        TCLAP::SwitchArg _logo("", "logo", "prints logo and exits", false,
                               new LogoVisitor(logo));
        cmd.add(_logo);

		cmd.parse( argc, argv );
		
		//check the arguments...
		
		inputFileName = _input.getValue();
        outputFileName = _out.getValue();

//        json = _json.getValue();
//        verbose = _verbose.getValue();

        if (_logo.getValue()) {
            std::cout << logo;
            exit(0);
        }

	}
	catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
	catch (...) {
		std::cerr << "unknown error" << std::endl;
		exit(0);
	}

    init();
    setArguments(argc, argv);
#ifndef _WIN32
	info->printInfo("Host name", hostname);

#endif // !_WIN32
    info->printInfo("Started on",startingTime);
    info->printInfo("Command-line arguments",arguments);

}

bool is_number(const std::string& s)
{
	std::locale loc;
    std::string::const_iterator it = s.begin();
	while (it != s.end() && std::isdigit(*it, loc)) ++it;
    return !s.empty() && it == s.end();
}

void cParameters::init() {

}

std::string cParameters::getChangeLog() {

    std::stringstream out;

    for (auto cl = changeLog.begin(); cl != changeLog.end(); ++ cl) {

        out <<"Version "<<cl->version<<" from "<<cl->date<<":\n";
        out <<"\t"<<cl->log<<"\n";
    }

    return out.str();
    
}

void cParameters::setArguments(int argc, const char * const * argv) {

    for (int i = 0; i < argc; i++) {
        arguments +=argv[i];
        arguments +=" ";
    }
    
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    startingTime = ctime(&now);
    startingTime.pop_back();
#ifndef _WIN32
    gethostname(hostname, sizeof(hostname));
#endif // !_WIN32
}

void cParameters::setLogo() {
        logo =
"                                    *               .*,                 \n                    /*            .(#/            .(,                   \n                     *#           (/#(          .(/,                    \n                      /%,        *//*#. ,.     **,.   .                 ...,,*****/(###%%%(/&&&&&%%#(/****,,,..                \n                    .....,,,,,,,,,,,,,,,,,,,,,......                    \n";
}
