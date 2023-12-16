// -*- Mode: c++; c-basic-offset: 4; tab-width: 4; -*-

/****************************************************************************** 
 * 
 *  file:  StdOutput.h
 * 
 *  Copyright (c) 2004, Michael E. Smoot
 *  All rights reverved.
 * 
 *  See the file COPYING in the top directory of this distribution for
 *  more information.
 *  
 *  THE SOFTWARE IS PROVIDED _AS IS_, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 *  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 *  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 *  DEALINGS IN THE SOFTWARE.  
 *  
 *****************************************************************************/ 

#ifndef TCLAP_STDCMDLINEOUTPUT_H
#define TCLAP_STDCMDLINEOUTPUT_H

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

#include <tclap/CmdLineInterface.h>
#include <tclap/CmdLineOutput.h>
#include <tclap/XorHandler.h>
#include <tclap/Arg.h>

#ifdef _WIN32
#include <Windows.h>
#else
//#include "utf8.h"
#endif

namespace TCLAP {

/**
 * A class that isolates any output from the CmdLine object so that it
 * may be easily modified.
 */
class StdOutput : public CmdLineOutput
{

	public:

		/**
		 * Prints the usage to stdout.  Can be overridden to 
		 * produce alternative behavior.
		 * \param c - The CmdLine object the output is generated for. 
		 */
		virtual void usage(CmdLineInterface& c);

		/**
		 * Prints the version to stdout. Can be overridden 
		 * to produce alternative behavior.
		 * \param c - The CmdLine object the output is generated for. 
		 */
		virtual void version(CmdLineInterface& c);

    /**
     * Prints the version to stdout. Can be overridden
     * to produce alternative behavior.
     * \param c - The CmdLine object the output is generated for.
     */
    virtual void changelog(CmdLineInterface& c);

    /**
		 * Prints (to stderr) an error message, short usage 
		 * Can be overridden to produce alternative behavior.
		 * \param c - The CmdLine object the output is generated for. 
		 * \param e - The ArgException that caused the failure. 
		 */
		virtual void failure(CmdLineInterface& c, 
				     ArgException& e );

	protected:

        /**
         * Writes a brief usage message with short args.
		 * \param c - The CmdLine object the output is generated for. 
         * \param os - The stream to write the message to.
         */
        void _shortUsage( CmdLineInterface& c, std::ostream& os ) const;

        /**
		 * Writes a longer usage message with long and short args, 
		 * provides descriptions and prints message.
		 * \param c - The CmdLine object the output is generated for. 
		 * \param os - The stream to write the message to.
		 */
		void _longUsage( CmdLineInterface& c, std::ostream& os ) const;

		/**
		 * This function inserts line breaks and indents long strings 
		 * according the  params input. It will only break lines at spaces, 
		 * commas and pipes.
		 * \param os - The stream to be printed to.
		 * \param s - The string to be printed.
		 * \param maxWidth - The maxWidth allowed for the output line. 
		 * \param indentSpaces - The number of spaces to indent the first line. 
		 * \param secondLineOffset - The number of spaces to indent the second
		 * and all subsequent lines in addition to indentSpaces.
		 */
		void spacePrint( std::ostream& os, 
						 const std::string& s, 
						 int maxWidth, 
						 int indentSpaces, 
						 int secondLineOffset ) const;

		void spacePrint(std::ostream& os,
			const std::wstring& s,
			int maxWidth,
			int indentSpaces,
			int secondLineOffset) const;
};


inline void StdOutput::version(CmdLineInterface& _cmd) 
{
	std::string progName = _cmd.getProgramName();
	std::string xversion = _cmd.getVersion();

	std::cout << std::endl << progName << "  version: " 
			  << xversion << std::endl << std::endl;
}

inline void StdOutput::changelog(CmdLineInterface& _cmd)
    {
    std::string progName = _cmd.getProgramName();
    std::string xversion = _cmd.getVersion();
    std::string xchangelog = _cmd.getChangeLog();

    std::cout << std::endl << progName << "  ChangeLog:\n"
    << xchangelog << std::endl;
}

inline void StdOutput::usage(CmdLineInterface& _cmd ) 
{
	std::cout << std::endl << "USAGE: " << std::endl << std::endl; 

	_shortUsage( _cmd, std::cout );

	std::cout << std::endl << std::endl << "Where: " << std::endl << std::endl;

	_longUsage( _cmd, std::cout );

	//std::cout << std::endl;

}

inline void StdOutput::failure( CmdLineInterface& _cmd,
								ArgException& e ) 
{
	std::string progName = _cmd.getProgramName();

	std::cerr << "PARSE ERROR: " << e.argId() << std::endl
		      << "             " << e.error() << std::endl << std::endl;

	if ( _cmd.hasHelpAndVersion() )
		{
			std::cerr << "Brief USAGE: " << std::endl;

			_shortUsage( _cmd, std::cerr );	

			std::cerr << std::endl << "For complete USAGE and HELP type: " 
					  << std::endl << "   " << progName << " --help" 
					  << std::endl << std::endl;
		}
	else
		usage(_cmd);

	throw ExitException(1);
}

inline void 
StdOutput::_shortUsage( CmdLineInterface& _cmd, 
						std::ostream& os ) const
{
	std::list<Arg*> argList = _cmd.getArgList();
	std::string progName = _cmd.getProgramName();
	XorHandler xorHandler = _cmd.getXorHandler();
	std::vector< std::vector<Arg*> > xorList = xorHandler.getXorList();

	std::string s = progName + " ";

	// first the xor
	for ( int i = 0; static_cast<unsigned int>(i) < xorList.size(); i++ )
		{
			s += " {";
			for ( ArgVectorIterator it = xorList[i].begin(); 
				  it != xorList[i].end(); it++ )
				s += (*it)->shortID() + "|";

			s[s.length()-1] = '}';
		}

	// then the rest
	for (ArgListIterator it = argList.begin(); it != argList.end(); it++)
		if ( !xorHandler.contains( (*it) ) && (*it)->getName() != "version" && (*it)->getName() != "help" && (*it)->getName() != "log" && !(*it)->isHidden())
			s += " " + (*it)->shortID();

	// then version and help
	for (ArgListIterator it = argList.begin(); it != argList.end(); it++)
		if (  (*it)->getName() == "version" || (*it)->getName() == "help" || (*it)->getName() == "log")
			s += " " + (*it)->shortID();
	
	// if the program name is too long, then adjust the second line offset
	int secondLineOffset = static_cast<int>(progName.length()) + 2;
	if ( secondLineOffset > 75/2 )
		secondLineOffset = static_cast<int>(75/2);

	spacePrint( os, s, 75, 3, secondLineOffset );
}

inline void 
StdOutput::_longUsage( CmdLineInterface& _cmd, 
					   std::ostream& os ) const
{
	std::list<Arg*> argList = _cmd.getArgList();
	std::string message = _cmd.getMessage();
	XorHandler xorHandler = _cmd.getXorHandler();
	std::vector< std::vector<Arg*> > xorList = xorHandler.getXorList();

	// first the xor 
	for ( int i = 0; static_cast<unsigned int>(i) < xorList.size(); i++ )
		{
			for ( ArgVectorIterator it = xorList[i].begin(); 
				  it != xorList[i].end(); 
				  it++ )
				{
					spacePrint( os, (*it)->longID(), 75, 3, 3 );
					spacePrint( os, (*it)->getDescription(), 75, 5, 0 );

					if ( it+1 != xorList[i].end() )
						spacePrint(os, "-- OR --", 75, 9, 0);
				}
			os << std::endl << std::endl;
		}

	// then the rest
	for (ArgListIterator it = argList.begin(); it != argList.end(); it++)
		if ( !xorHandler.contains( (*it) )  && (*it)->getName() != "version" && (*it)->getName() != "help" && (*it)->getName() != "log" && !(*it)->isHidden())
			{
				spacePrint(os, (*it)->longID(), 75, 3, 3);
				if ((*it)->_description != "") {
					spacePrint(os, (*it)->getDescription(), 75, 5, 0);
					os << std::endl;
				}
				if ((*it)->_wdescription != L"") {
					spacePrint(os, (*it)->getWDescription(), 75, 5, 0);
					os << std::endl;
				}
			}

	// then version and help
	for (ArgListIterator it = argList.begin(); it != argList.end(); it++)
		if (  (*it)->getName() == "version" || (*it)->getName() == "help" || (*it)->getName() == "log")
			{
			spacePrint(os, (*it)->longID(), 75, 3, 3);
			spacePrint(os, (*it)->getDescription(), 75, 5, 0);
			os << std::endl;
		}
	
	//os << std::endl;

	if (message!="")
		spacePrint( os, message, 75, 3, 0 );
}

inline void StdOutput::spacePrint( std::ostream& os, 
						           const std::string& s, 
						           int maxWidth, 
						           int indentSpaces, 
						           int secondLineOffset ) const
{
	int len = static_cast<int>(s.length());

	if ( (len + indentSpaces > maxWidth) && maxWidth > 0 )
		{
			int allowedLen = maxWidth - indentSpaces;
			int start = 0;
			while ( start < len )
				{
					// find the substring length
					// int stringLen = std::min<int>( len - start, allowedLen );
					// doing it this way to support a VisualC++ 2005 bug 
					using namespace std; 
					int stringLen = min<int>( len - start, allowedLen );

					// trim the length so it doesn't end in middle of a word
					if ( stringLen == allowedLen )
						while ( stringLen >= 0 &&
								s[stringLen+start] != ' ' && 
								s[stringLen+start] != ',' &&
								s[stringLen+start] != '|' ) 
							stringLen--;
	
					// ok, the word is longer than the line, so just split 
					// wherever the line ends
					if ( stringLen <= 0 )
						stringLen = allowedLen;

					// check for newlines
					for ( int i = 0; i < stringLen; i++ )
						if ( s[start+i] == '\n' )
							stringLen = i+1;

					// print the indent	
					for ( int i = 0; i < indentSpaces; i++ )
						os << " ";

					if ( start == 0 )
						{
							// handle second line offsets
							indentSpaces += secondLineOffset;

							// adjust allowed len
							allowedLen -= secondLineOffset;
						}

					os << s.substr(start,stringLen) << std::endl;

					// so we don't start a line with a space
					while ( s[stringLen+start] == ' ' && start < len )
						start++;
			
					start += stringLen;
				}
		}
	else
		{
			for ( int i = 0; i < indentSpaces; i++ )
				os << " ";
			os << s << std::endl;
		}
}

inline void StdOutput::spacePrint(std::ostream& os,
	const std::wstring& s,
	int maxWidth,
	int indentSpaces,
	int secondLineOffset) const
{
	int len = static_cast<int>(s.length());

	if ((len + indentSpaces > maxWidth) && maxWidth > 0)
	{
		int allowedLen = maxWidth - indentSpaces;
		int start = 0;
		while (start < len)
		{
			// find the substring length
			// int stringLen = std::min<int>( len - start, allowedLen );
			// doing it this way to support a VisualC++ 2005 bug 
			using namespace std;
			int stringLen = min<int>(len - start, allowedLen);

			// trim the length so it doesn't end in middle of a word
			if (stringLen == allowedLen)
			while (stringLen >= 0 &&
				s[stringLen + start] != L' ' &&
				s[stringLen + start] != L',' &&
				s[stringLen + start] != L'|')
				stringLen--;

			// ok, the word is longer than the line, so just split 
			// wherever the line ends
			if (stringLen <= 0)
				stringLen = allowedLen;

			// check for newlines
			for (int i = 0; i < stringLen; i++)
			if (s[start + i] == L'\n')
				stringLen = i+1;

			// print the indent	
			for (int i = 0; i < indentSpaces; i++)
				os<<" ";

			if (start == 0)
			{
				// handle second line offsets
				indentSpaces += secondLineOffset;

				// adjust allowed len
				allowedLen -= secondLineOffset;
			}
			std::wstring sub = s.substr(start, stringLen);
			const wchar_t *s2 =sub.c_str();
#ifdef _WIN32
			int bufferSize2 = WideCharToMultiByte(CP_UTF8, 0,
				s2, -1, NULL, 0, NULL, NULL);
			char* m2 = new char[bufferSize2];
			WideCharToMultiByte(CP_UTF8, 0, s2, -1, m2,
				bufferSize2, NULL, NULL);

			wprintf(L"%ls", m2);
			delete m2;
#else
			wprintf(s2);

#endif
			os << std::endl;

			//os << s.substr(start, stringLen) << std::endl;

			// so we don't start a line with a space
			while (s[stringLen + start] == L' ' && start < len)
				start++;

			start += stringLen;
		}
	}
	else
	{
		for (int i = 0; i < indentSpaces; i++)
			os << " ";
	const wchar_t *s2 = s.c_str();
#ifdef _WIN32
		int bufferSize2 = WideCharToMultiByte(CP_UTF8, 0,
			s2, -1, NULL, 0, NULL, NULL);
		char* m2 = new char[bufferSize2];
		WideCharToMultiByte(CP_UTF8, 0, s2, -1, m2,
			bufferSize2, NULL, NULL);

		wprintf(L"%ls", m2);
		delete m2;
#else
	wprintf(s2);
#endif
	
		os << std::endl;
	}
	
}

} //namespace TCLAP
#endif 
