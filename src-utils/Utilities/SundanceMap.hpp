/* @HEADER@ */
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
      inline bool containsKey(const Key& key) const {return find(key) != end();}

      /** */
      inline void put(const Key& key, const Value& value)
        {operator[](key) = value;}

      /** */
      inline const Value& get(const Key& key) const
        {return (*find(key)).second;}

      /** */
      inline Value& get(const Key& key) 
        {return (*find(key)).second;}
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
