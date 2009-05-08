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

#ifndef SUNDANCE_QUADRATUREFAMILYSTUB_H
#define SUNDANCE_QUADRATUREFAMILYSTUB_H

#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_XMLObject.hpp"

namespace SundanceCore
{
using namespace SundanceUtils;
using namespace Teuchos;

class QuadratureFamilyStub 
  : public TSFExtended::Handleable<QuadratureFamilyStub>,
    public TSFExtended::Printable,
    public TSFExtended::Describable,
    public Noncopyable
{
public:
  /** Empty ctor */
  QuadratureFamilyStub(int order);

  /** virtual dtor */
  virtual ~QuadratureFamilyStub(){;}

  /** Return the order of the quadrature rule */
  int order() const {return order_;}

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** Ordering for storage in STL maps */
  virtual bool lessThan(const QuadratureFamilyStub* other) const 
    {return order() < other->order();}

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const {os << toXML();}
  //@}

  /** \name Describable interface */
  //@{
  /** Print to a stream */
  virtual string description() const 
    {return "QuadratureFamilyStub[order=" + Teuchos::toString(order()) 
        +  "]";}
  //@}

  /** */
  virtual RefCountPtr<QuadratureFamilyStub> getRcp() {return rcp(this);}

  /** */
  static RefCountPtr<QuadratureFamilyStub>& defaultQuadrature()
    {static RefCountPtr<QuadratureFamilyStub> rtn
        = rcp(new QuadratureFamilyStub(1)); return rtn;}
      
private:
  int order_;
};
}


#endif
