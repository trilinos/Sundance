/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MULTISET_H
#define SUNDANCE_MULTISET_H

#include "SundanceDefs.hpp"
#include <set>

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceUtils
{
  using namespace Teuchos;

  /** */
  template<class Key>
    class MultiSet : public std::multiset<Key>
    {
    public:
      /** */
      MultiSet() : std::multiset<Key>() {;}

      /** */
      bool contains(const Key& key) const {return find(key) != end();}

      /** */
      void put(const Key& key) {insert(key);}

      /** */
      Array<Key> elements() const ;

      /** */
      ostream& toStream(ostream& os) const ;

      /** */
      string toString() const ;
    };


  template<class Key> inline
    Array<Key> MultiSet<Key>::elements() const
    {
      Array<Key> rtn;

      typename MultiSet<Key>::const_iterator iter;

      for (iter=begin(); iter != end(); iter++)
        {
          rtn.append(*iter);
        }
      return rtn;
    }

  template<class Key> inline
    ostream& MultiSet<Key>::toStream(ostream& os) const
    {
      typename MultiSet<Key>::const_iterator iter;

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

  template<class Key> inline
    string MultiSet<Key>::toString() const
    {
      ostringstream os;
      os << *this;
      return os.str();
    }

}

namespace std
{
/** \relates MultiSet */
  template<class Key> inline
    ostream& operator<<(ostream& os, const SundanceUtils::MultiSet<Key>& m)
    {return m.toStream(os);}
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif

