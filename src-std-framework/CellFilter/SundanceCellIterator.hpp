/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLITERATOR_H
#define SUNDANCE_CELLITERATOR_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceCellReordererImplemBase.hpp"
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
     * CellIterator is an iterator for walking through cell sets.
     * It satisfies the requirements for an input iterator in STL.
     *
     * This class design violates the usual rules of good OO style:
     * it has polymorphic behavior packed into a single class. The
     * justification for this decision is to avoid the expense
     * of the clone() operations that would be required by the
     * copy ctor for iterators were 
     * a polymorphic class heirarchy used. 
     *
     * Two cell set types exist: explicit, where the member cells LIDs
     * are enumerated in a physical Set<int> object, and implicit,
     * where no physical Set is made, rather, the sequence of cell LIDs
     * is obtained through some scheme of walking the mesh. 
     *
     * The only cell sets that can represented implicitly are
     * the set of all cells of a given dimension.
     * 
     * \see CellSet, CellFilter, CellIteratorPos
     */
    class CellIterator : public std::iterator<input_iterator_tag, int>
    {
    public:

      
      /** 
       * CellIteratorPos is used to specify whether a new CellIterator
       * is positioned at the beginning or end of a set.
       */
      enum CellIteratorPos {Begin, End};

      /** Empty ctor */
      CellIterator();

      /** Construct an implicit iterator for walking all cells of a given
       * dimension on the given mesh. */
      CellIterator(const Mesh& mesh, int cellDim, CellIteratorPos pos);

      /** Construct an explicit iterator for walking an explicitly
       * enumerated set of cells. */
      CellIterator(const Set<int>* cells, CellIteratorPos pos);
      
      /** Dereferencing operator */
      const int& operator*() const 
      {
        if (isImplicit_) return currentLID_;
        else return *iter_;
      }
      
      /** Postfix increment: advances iterator and returns previous value  */
      CellIterator operator++(int) 
      {
        CellIterator old = *this;
        advance();
        return old;
      }
      

      /** Prefix increment: advances iterator, returning new value */
      CellIterator& operator++()
      {
        advance();
        return *this;
      }

      /** */
      bool operator==(const CellIterator& other) const 
      {
        if (isImplicit_)
          {
            return currentLID_ == other.currentLID_;
          }
        else
          {
            return iter_ == other.iter_;
          }
      }

      /** */
      bool operator!=(const CellIterator& other) const 
      {
        return !(*this == other);
      }

      
    private:

      /** Advance the iterator */
      void advance()
      {
        CellIterator old = *this;
        if (isImplicit_) 
          {
            if (reorderer_ != 0) 
              {
                currentLID_ = reorderer_->advance(currentLID_);
              }
            else currentLID_++;
          }
        else iter_++;
      }
      
      /** Flag indicating whether this iterator is implicit */
      bool isImplicit_;
      
      /** The LID to which this iterator is currently pointing.
       * Used only for implicit iterators. */
      int currentLID_;

      /** Unmanaged pointer to the reorderer used for walking 
       * implicit cell sets. Used only for implicit iterators. */
      const CellReordererImplemBase* reorderer_; 

      /** iterator for enumerated cells.
       * Used only for explicit iterators. */
      Set<int>::const_iterator iter_;
    };
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
