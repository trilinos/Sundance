/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXPLICITCELLITERATOR_H
#define SUNDANCE_EXPLICITCELLITERATOR_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    /**
     *
     */
    class ExplicitCellIterator
    {
    public:
      /** */
      ExplicitCellIterator(const Set<int>::const_iterator& iter);
      
      /** Dereferencing operator */
      const int& operator*() const {return *iter_;}
      
      /** Postfix increment: advances iterator and returns previous value  */
      CellIterator operator++(int) 
      {
        CellIterator old = *this;
        iter_++;
        return old;
      }
      

      /** Prefix increment: advances iterator, returning new value */
      CellIterator& operator++()
      {
        iter_++;
        return *this;
      }

      /** */
      bool operator==(const CellIterator& other) const 
      {
        return iter_ == other.iter_;
      }

      /** */
      bool operator!=(const CellIterator& other) const 
      {
        return currentLID_ == other.currentLID_;
      }

      
    private:

      /** The LID to which this iterator is currently pointing. */
      int currentLID_;

      /** Unmanaged pointer to the cell set through which this iterator
       * is iterating */
      CellSetBase* cellSet_;

      
      
    }
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
