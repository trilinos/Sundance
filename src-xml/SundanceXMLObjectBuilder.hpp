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

#ifndef SUNDANCE_XMLOBJECTBUILDER_H
#define SUNDANCE_XMLOBJECTBUILDER_H

#include "Sundance.hpp"


namespace SundanceXML
{
  template <class T>
  class XMLObjectBuilder
  {
  public:
    /** */
    XMLObjectBuilder() : varMap_() {;}
    
    /** */
    bool hasName(const string& name) const
    {
      return varMap_.containsKey(name);
    }

    /** */
    const T& get(const string& name) const 
    {
      return varMap_.get(name);
    }

    /** */
    void checkOptionValidity(const string& option, 
                             const Set<string>& validChoices) const 
    {
      TEST_FOR_EXCEPTION(!valid.contains(key), RuntimeError,
                         "option [" << option << "] not found among valid "
                         "choices " << validChoices);
    }

    /** */
    void checkTag(const XMLObject& xml, const string& expectedTag) const 
    {
      TEST_FOR_EXCEPTION(xml.getTag() != expectedTag, RuntimeError,
                         "found tag [" << xml.getTag() << "] where ["
                         << expectedTag << "] was expected");
    }

  protected:
    
    /** */
    void addToMap(const string& name, const T& value) const
    {
      varMap_.put(name, value);
    }
  private:
    
    mutable Map<string, T> varMap_;

  };

}
#endif
