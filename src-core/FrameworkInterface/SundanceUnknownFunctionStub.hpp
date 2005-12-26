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

#ifndef SUNDANCE_UNKNOWNFUNCTIONSTUB_H
#define SUNDANCE_UNKNOWNFUNCTIONSTUB_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "SundanceUnknownFuncDataStub.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;

  namespace Internal
  {
    using std::string;
    using std::ostream;

    /** 
   * UnknownFunctionStub is the base class for unknown functions. 
   * Each framework will need to implement its own subclass of
   * UnknownFunctionStub. 
   *
   * The interface is left very minimal so as to not place
   * any constraints on how a framework might specify the basis.
   * When a framework needs any information about the
   * unknown function, it will have to get it by downcasting
   * to the appropriate framework-specific subclass.
   *
   * <h4> Writing a UnknownFunctionStub subclass </h4>
   *
   * For purposes of interaction with the Sundance core, no 
   * additional methods are required.
   * However, most frameworks will require extensions to 
   * UnknownFunctionStub that can supply the framework with information
   * on the basis used by the unknown func. See the
   * demo and standard frameworks for information on how to do this.
   */
    class UnknownFunctionStub : public SymbolicFunc
    {
    public:
      /** */
      UnknownFunctionStub(const string& name, int nElems=1,
                          const RefCountPtr<const UnknownFuncDataStub>& data=RefCountPtr<const UnknownFuncDataStub>());

      /** virtual destructor */
      virtual ~UnknownFunctionStub() {;}

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}

    protected:

      /** */
      const RefCountPtr<const UnknownFuncDataStub>& dataStub() const {return data_;}

    private:
      RefCountPtr<const UnknownFuncDataStub> data_;


    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
