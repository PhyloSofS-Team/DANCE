#pragma once

#include <string>
#include <vector>
#include <set>
#include  <sstream>

    class cInfo {

    public:
        enum formats {NONE=0, JSON, XML, YAML};

    private:
        int outLength;
        char buf[100];
        char buf1[100];

        double totalTime;

        std::vector<std::pair<std::string, double> > records;
        void	printHeader(const char *, int opt=0);

    public:
        formats format;
        std::stringstream formatStr;
        int	formatLevel;
        int prevFormatLevel;
        void indent();
//        void printFormat(const char *str);


        cInfo(int Length);

        void printInfo(const char *info);

        void printInfo(const char *info, const char *arg, const char *val =NULL);
        void printInfo(const char *info, const std::string &arg, const char *val =NULL);
        void printInfo(const char *info, int arg, const char *val =NULL);
        void printInfo(const char *info, int arg, double val);
        void printInfo(const char *info, float arg, const char *val =NULL);
        void printInfo(const char *info, double arg, const char *val =NULL);
        void printInfo(const char *info, const double * vec, int len, const char *val);
        void printInfo(const char *info, const std::set<std::string> &arg);


        void saveTimeInterval(double);
        void saveTimeInterval(const char *, double);
        void setTotalTime (double);
        void printTimings();

        void setFormat (formats f);
        void printLog (const std::string &fileName);
//        template <class In, class Arg, class Unit>
//        void printFormat(In *info, Arg arg, Unit *unit = NULL);

    public:
				inline void enableVerbose() { m_quiet = false; }
				inline void disableVerbose() { m_quiet = true; }
private:
				bool m_quiet;
        
    };

//template <class In, class Arg, class Unit>
//void cInfo::printFormat(In *info, Arg str, Unit *val ) {
//
//    switch (format) {
//        case JSON:
//            if (prevFormatLevel != formatLevel) {
//            } else {
//                formatStr << ",\n";
//            }
//
//            indent();
//            formatStr << "\""<<info<<"\": {\n";
//
//            formatLevel++;
//            indent();
//            formatStr << "\"value\": ";
//            formatStr << "\""<<str<<"\"";
//            if (val) {
//                formatStr << ", \"unit\": ";
//                formatStr << "\""<< val << "\"";
//            }
//            
//            formatStr << "\n";
//            formatLevel--;
//            prevFormatLevel = formatLevel;
//            indent();
//            formatStr << "}";
//            
//            break;
//    }
//}
//
