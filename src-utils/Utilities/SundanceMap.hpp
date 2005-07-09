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

#ifndef SUNDANCE_MAP_H
#define SUNDANCE_MAP_H

#include "SundanceDefs.hpp"
#include <map>

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  using namespace Teuchos;

  /** */
  template<class Key, class Value, class Compare = less<Key> >
    class Map : public std::map<Key, Value, Compare>
    {
    public:
      /** */
      Map() : std::map<Key,Value,Compare>() {;}

      /** */
      inline bool containsKey(const Key& key) const {return this->find(key) != this->end();}

      /** */
      inline void put(const Key& key, const Value& value)
        {operator[](key) = value;}

      /** */
      inline const Value& get(const Key& key) const
        {return (*(this->find)(key)).second;}

      /** */
      inline Value& get(const Key& key) 
        {return (*(this->find)(key)).second;}
    };

}

namespace std
{
   /** \relates Map */
  template<class Key, class Value, class Compare> inline
    std::ostream& operator<<(std::ostream& os, const SundanceUtils::Map<Key, Value, Compare>& m)
    {
      typename SundanceUtils::Map<Key, Value, Compare>::const_iterator iter;

      os << "Map[";
      unsigned int k = 0 ;
      for (iter=m.begin(); iter != m.end(); iter++, k++)
        {
          os << "{" << (*iter).first << ", " << (*iter).second << "}";
          if (k < m.size()-1) os << ", ";
        }
      os << "]";
      return os;
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
