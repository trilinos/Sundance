/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLFILTERBASE_H
#define SUNDANCE_CELLFILTERBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceCellSet.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;
  
  namespace Internal
  {
    /** 
     * Base class for CellFilter objects.
     *
     * <h4> Notes for subclass implementors </h4>
     * 
     * Derived classes must implement the methods
     * <ul>
     * <li> internalGetCells() -- returns the set of cells that 
     * pass through this filter
     * <li> dimension() -- returns the dimension of the cells that
     * will pass through this filter
     * <li> toXML() -- writes an XML description of the filter
     * <li> lessThan() -- compares to another cell filter. Used to store
     * cell filters in STL containers. 
     * <li> typeName() -- returns the name of the subclass. Used in ordering.
     * </ul>
     */
    class CellFilterBase : public CellFilterStub
    {
    public:
      /** Empty ctor */
      CellFilterBase();

      /** virtual dtor */
      virtual ~CellFilterBase(){;}

      /** Find the cells passing this filter on the given mesh. This
       * method will cache the cell sets it computes for each mesh  */
      CellSet getCells(const Mesh& mesh) const ;

      /** Return the dimension of the cells that will be identified
       * by this filter when acting on the given mesh */
      virtual int dimension(const Mesh& mesh) const = 0 ;

    protected:

      /** */
      virtual CellSet internalGetCells(const Mesh& mesh) const = 0 ;

    private:
      /** cache of previously computed cell sets */
      mutable CellSet cellSetCache_;
      //      mutable Map<int, CellSet> cellSetCache_;

    };
  }
}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
