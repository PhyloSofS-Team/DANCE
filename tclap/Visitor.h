
/****************************************************************************** 
 * 
 *  file:  Visitor.h
 * 
 *  Copyright (c) 2003, Michael E. Smoot .
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


#ifndef TCLAP_VISITOR_H
#define TCLAP_VISITOR_H

#include <tclap/CmdLineInterface.h>
//#include <tclap/CmdLineOutput.h>

namespace TCLAP {

/**
 * A base class that defines the interface for visitors.
 */
class Visitor
{
	public:

		/**
		 * Constructor. Does nothing.
		 */
		Visitor() { }

		/**
		 * Destructor. Does nothing.
		 */
		virtual ~Visitor() { }

		/**
		 * Does nothing. Should be overridden by child.
		 */
		virtual void visit() { }
};


    class EmptyVisitor: public Visitor
        {
        private:
        /**
         * Prevent accidental copying.
         */
        EmptyVisitor(const EmptyVisitor& rhs);
        EmptyVisitor& operator=(const EmptyVisitor& rhs);

        public:

        /**
         * Constructor.
         * \param cmd - The CmdLine the output will be generated for.
         * \param out - The type of output.
         */
        EmptyVisitor()
        : Visitor() { }
        
        /**
         * Calls the usage method of the CmdLineOutput for the 
         * specified CmdLine.
         */
        void visit() {
            return;
        }
        
        };
        

}

#endif
