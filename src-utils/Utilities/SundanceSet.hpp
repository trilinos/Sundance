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

#ifndef SUNDANCE_SET_H
#define SUNDANCE_SET_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include <set>

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceUtils
{
  using namespace Teuchos;

  /** */
  template<class Key, class Compare = less<Key> >
  class Set : public std::set<Key, Compare>
  {
  public:
    /** */
    Set() : std::set<Key, Compare>() {;}

    /** */
    bool contains(const Key& key) const {return this->find(key) != this->end();}

    /** */
    void put(const Key& key) {insert(key);}

    /** */
    Array<Key> elements() const ;

    /** */
    void elements(Array<Key>& keys) const ;

    /** */
    void merge(const Set<Key, Compare>& other);

    /** */
    ostream& toStream(ostream& os) const ;

    /** */
    string toString() const ;
  };


  template<class Key, class Compare> inline
  Array<Key> Set<Key, Compare>::elements() const
  {
    Array<Key> rtn;

    typename Set<Key, Compare>::const_iterator iter;

    for (iter=this->begin(); iter != this->end(); iter++)
      {
        rtn.append(*iter);
      }
    return rtn;
  }


  template<class Key, class Compare> inline
  void Set<Key, Compare>::elements(Array<Key>& rtn) const
  {
    rtn.resize(0);
    typename Set<Key, Compare>::const_iterator iter;

    for (iter=this->begin(); iter != this->end(); iter++)
      {
        rtn.append(*iter);
      }
  }

  template<class Key, class Compare> inline
  void Set<Key, Compare>::merge(const Set<Key, Compare>& other)
  {
    typename Set<Key, Compare>::const_iterator iter;

    for (iter=other.begin(); iter != other.end(); iter++)
      {
        put(*iter);
      }
  }

  template<class Key, class Compare> inline
  ostream& Set<Key, Compare>::toStream(ostream& os) const
  {
    typename Set<Key, Compare>::const_iterator iter;

    unsigned int k = 0;
    os << "{";
    for (iter=this->begin(); iter != this->end(); iter++, k++)
      {
        os << *iter;
        if (k<(this->size()-1)) os << ", ";
      }
    os << "}";

    return os;
  }

  template<class Key, class Compare> inline
  string Set<Key, Compare>::toString() const
  {
    ostringstream os;
    os << *this;
    return os.str();
  }

}

namespace std
{
  /** \relates Set */
  template<class Key, class Compare> inline
  ostream& operator<<(ostream& os, const SundanceUtils::Set<Key, Compare>& m)
  {return m.toStream(os);}
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
