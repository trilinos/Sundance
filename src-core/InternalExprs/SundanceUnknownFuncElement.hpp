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

#ifndef SUNDANCE_UNKNOWNFUNCELEMENT_H
#define SUNDANCE_UNKNOWNFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceUnknownFuncDataStub.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;


  namespace Internal
  {
    using namespace Teuchos;

    using std::string;
    using std::ostream;

    /** 
     * UnknownFuncElement represents a scalar-valued element of a (possibly)
     * list-valued UnknownFunction
     */
    class UnknownFuncElement : public SymbolicFuncElement
    {
    public:
      /** */
      UnknownFuncElement(const RefCountPtr<const UnknownFuncDataStub>& data,
                         const string& name,
                         const string& suffix,
                         int myIndex);

      /** virtual destructor */
      virtual ~UnknownFuncElement() {;}

      /** Get the data associated with the vector-valued function 
       * that contains this function element. */
      const UnknownFuncDataStub* commonData() const {return commonData_.get();}


      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual RefCountPtr<Internal::ExprBase> getRcp() {return rcp(this);}
      
    private:
       const RefCountPtr<const UnknownFuncDataStub> commonData_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
