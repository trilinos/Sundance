/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#ifndef SUNDANCE_TESTFUNCELEMENT_H
#define SUNDANCE_TESTFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFuncElement.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    class TestFunctionStub;
  }

  namespace Internal
  {
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * TestFuncElement represents a scalar-valued element of a (possibly)
     * list-valued TestFunction
     */
    class TestFuncElement : public SymbolicFuncElement
    {
    public:
      /** */
      TestFuncElement(const TestFunctionStub* master,
                      const string& name,
                      const string& suffix,
                      int myIndex);

      /** virtual destructor */
      virtual ~TestFuncElement() {;}

      /** Get the master test function 
       * of which this object is an element */
      const TestFunctionStub* master() const {return master_;}


      /** Test whether all terms have test functions. 
       * I'm a test function, so return true */
      virtual bool allTermsHaveTestFunctions() const {return true;}

      /** Test whether this expr contains a test function. 
       * I'm a test function, so return true. */
      virtual bool hasTestFunctions() const {return true;}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
      const TestFunctionStub* master_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
