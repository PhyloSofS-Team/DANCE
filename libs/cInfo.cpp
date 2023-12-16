#include "cInfo.hpp"
#include <string.h>
#include <stdio.h>
#include <set>

cInfo::cInfo (int l) : outLength(l)
{
	strcpy(buf,"..................................................................................................\0");
	strcpy(buf1,"==================================================================================================\0");
	m_quiet = false;
    format = NONE;
}

void cInfo::printInfo(const char *info) { if (!m_quiet) {
	int i1 = strlen(info);
	int i2 = (outLength - i1)/2;
	int i3 = i2;
	if (i2*2+i1 != outLength)
		i3++;
	printf("%.*s%s%.*s\n",i2,buf1,info,i3,buf1);
	records.push_back(std::pair<std::string, double>(info,0.0));
//    if (format) {printFormat(info);}
}}

void cInfo::printInfo(const char *info, const char *arg, const char *val) { if (!m_quiet) {
	int i1 = strlen(info);
	int i3 = outLength - i1;
	printf("%s%.*s : %s",info,i3,buf,arg);
	if (val) printf(" %s",val);
	printf("\n");
//    if (format) {printFormat(info, arg, val);}
}}

void cInfo::printInfo(const char *info, const std::string &arg, const char *val) { if (!m_quiet) {
	int i1 = strlen(info);
	int i3 = outLength - i1;
	printf("%s%.*s : %s",info,i3,buf,arg.c_str());
	if (val) printf(" %s",val);
	printf("\n");
//    if (format) {printFormat(info, arg, val);}
}}

void cInfo::printInfo(const char *info, int arg, const char *val) { if (!m_quiet) {
	int i1 = strlen(info);
	int i3 = outLength - i1;
	printf("%s%.*s : %d",info,i3,buf,arg);
	if (val) printf(" %s",val);
	printf("\n");
}}

void cInfo::printInfo(const char *info, int arg, double val) { if (!m_quiet) {
    int i1 = strlen(info);
    int i3 = outLength - i1;
    printf("%s%.*s : %d",info,i3,buf,arg);
    printf(", %f",val);
    printf("\n");
}}

void cInfo::printInfo(const char *info, float arg, const char *val) { if (!m_quiet) {
	int i1 = strlen(info);
	int i3 = outLength - i1;
	printf("%s%.*s : %g",info,i3,buf,arg);
	if (val) printf(" %s",val);
	printf("\n");
//    if (format) {printFormat(info, arg, val);}
}}

void cInfo::printInfo(const char *info, double arg, const char *val) { if (!m_quiet) {
	int i1 = strlen(info);
	int i3 = outLength - i1;
	printf("%s%.*s : %g",info,i3,buf,arg);
	if (val) printf(" %s",val);
	printf("\n");
//    if (format) {printFormat(info, arg, val);}
}}

void cInfo::printInfo(const char *info, const double * vec, int len, const char *val) { if (!m_quiet) {
    int i1 = strlen(info);
    int i3 = outLength - i1;
    printf("%s :\n",info);

    const char *whiteBuf = "                                                              ";
    for (int i=0; i<len; ++i) {

        int i2 = 2+2+1;
        printf("%.*s[%2d]%.*s : %g",i1+1,whiteBuf,i+1, i3-i2,buf,vec[i]);

        if (val) printf(" %s",val);
        printf("\n");
    }
//    if (format) {
//        printFormat(info, arg, val);
//    }
}}

void cInfo::saveTimeInterval(double dt) {
	records.back().second = dt;
}

void cInfo::saveTimeInterval(const char *info, double dt) {
	records.push_back(std::pair<std::string, double>(info, dt));
}

void cInfo::setTotalTime (double t) {
	totalTime = t;
}

void cInfo::printHeader(const char *info, int opt) {
	
	if (opt == 0) printf("%.*s\n",outLength,buf1);
	else if (opt == 1) printf("%.*s\n",outLength,buf);
	
	int i1 = strlen(info);
	int i2 = (outLength - i1)/2;
	int i3 = i2;
	if (i2*2+i1 != outLength)
		i3++;
	if (opt == 0)
		printf("%.*s%s%.*s\n",i2,buf1,info,i3,buf1);
	else if (opt == 1)
		printf("%s%.*s : %g s\n",info,i2+i3,buf,totalTime);
	printf("%.*s\n",outLength,buf1);
}

void cInfo::printTimings() {

	printHeader(" Timing : ");

	for (std::size_t i =0; i < records.size(); ++i) {
		std::string &str = records[i].first;
		double dt = records[i].second;
		
		int i1 = str.size();
		int i3 = outLength - i1;
		printf("%s%.*s : %g s\n",str.c_str(),i3,buf,dt);

	}
	printHeader("Total time : ",1);
}

void cInfo::setFormat(cInfo::formats f) {

    format = f;
    formatLevel = 0;

    switch (format) {
        case JSON:
            formatStr << "{\n";
            formatLevel++;
            indent();
            formatStr << "\"AnAnaS\": {\n";
            formatLevel++;
            break;
        case XML:
            formatStr << "<AnAnaS>\n";
            formatLevel++;
            break;

        default:
            break;
    }
}

void cInfo::indent() {
    for (int i=0; i<4*formatLevel; i++)
        formatStr <<" ";
}

//void cInfo::printFormat(const char *str) {
//
//    switch (format) {
//        case JSON:
//
//            if (formatLevel > 2) {
//                formatStr << "\n";
//                formatLevel--;
//                indent();
//                formatStr << "}";
//            }
//
//            if (prevFormatLevel >= formatLevel) {
//                formatStr << ",\n";
//            }
//            indent();
//
//            formatStr << "\""<<str<<"\": {\n";
//            prevFormatLevel = formatLevel;
//            formatLevel++;
//            break;
//    }
//
//}

void cInfo::printLog (const std::string &fileName) {

    for(formatLevel--; formatLevel>=0; formatLevel--) {
        formatStr<<"\n";
        indent();
        formatStr<<"}";
    }
    //	std::cout <<formatStr.str()<<std::endl;

    FILE* outFile = fopen(fileName.c_str(), "w");

    if (!outFile) {
        printInfo("Cannot open file for writing",fileName);
    } else {
        fputs (formatStr.str().c_str(),outFile);
        fclose(outFile);
    }
    
}

void cInfo::printInfo(const char *info, const std::set<std::string> &arg) { if (!m_quiet) {
    int i1 = strlen(info);
    int i3 = outLength - i1;

    printf("%s%.*s : ",info,i3,buf);
    bool first = 1;
    for (auto item : arg) {
        if (!first) {
            printf(",");
        }
        first = 0;
        printf("%s",item.c_str());
    }
    printf("\n");
    //    if (format) {printFormat(info, arg, val);}
}}
