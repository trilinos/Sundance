/* @HEADER@ */
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
    bool contains(const Key& key) const {return find(key) != end();}

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

    for (iter=begin(); iter != end(); iter++)
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

    for (iter=begin(); iter != end(); iter++)
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
    for (iter=begin(); iter != end(); iter++, k++)
      {
        os << *iter;
        if (k<(size()-1)) os << ", ";
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
