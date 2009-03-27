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

#ifndef SUNDANCE_CELLFILTERSTUB_H
#define SUNDANCE_CELLFILTERSTUB_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceObjectWithInstanceID.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "Teuchos_XMLObject.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;
  using namespace Internal;
  
  namespace Internal
  {
    /** 
     * Stub class for cell filter objects, i.e., objects that can 
     * select a subset of mesh cells on which an integral or 
     * BC is to be applied.
     *
     * <h4> Notes for framework interface implementors </h4>
     *
     *  
     */
    class CellFilterStub : public TSFExtended::Handleable<CellFilterStub>,
                           public TSFExtended::Printable,
                           public Teuchos::Describable,
                           public Noncopyable,
                           public ObjectWithInstanceID<CellFilterStub>,
                           public TSFExtended::ObjectWithVerbosity<CellFilterStub>
    {
    public:
      /** Empty ctor */
      CellFilterStub();

      /** virtual dtor */
      virtual ~CellFilterStub(){;}

      /** Write to XML */
      virtual XMLObject toXML() const ;

      /** Ordering for storage in STL maps */
      virtual bool lessThan(const CellFilterStub* other) const 
      {return id() < other->id();}

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(std::ostream& os) const {os << toXML();}
      //@}

      /** \name Describable interface */
      //@{
      /** Print to a stream */
      virtual string description() const 
      {return "CellFilterStub[id=" + Teuchos::toString(id()) + "]";}
      //@}

      /** */
      virtual RefCountPtr<CellFilterStub> makeNullRegion() const ;

      /* */
      GET_RCP(CellFilterStub);

      /** */
      bool isNullRegion() const ;

      /** */
      bool operator!=(const CellFilterStub& other) const 
        {
          return this->lessThan(&other) || other.lessThan(this);
        }

      /** */
      bool operator==(const CellFilterStub& other) const 
        {
          return !(*this != other);
        }

    };
  }
}


#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
