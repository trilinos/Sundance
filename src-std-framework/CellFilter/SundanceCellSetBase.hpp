/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLSETBASE_H
#define SUNDANCE_CELLSETBASE_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceCellIterator.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"
#include "SundanceNoncopyable.hpp"
#include "TSFHandleable.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * CellSetBase is the base class for cell sets. There are two cell
     * set subtypes: ExplicitCellSet and ImplicitCellSet.
     *
     * @see CellFilter
     **/
    class CellSetBase : public TSFExtended::ObjectWithVerbosity<CellSetBase>,
                        public TSFExtended::Printable,
                        public Noncopyable,
                        public TSFExtended::Handleable<CellSetBase>
    {
    public:
      /** Construct, initializing to an empty set */
      CellSetBase(const Mesh& mesh, int cellDim,
                  const CellType& cellType);

      /** Return an iterator pointing to the first element in the set */
      virtual CellIterator begin() const = 0 ;

      /** Return an iterator containing the past-the-end value */
      virtual CellIterator end() const = 0 ;

      /** Return the type of cells in this set */
      const CellType& cellType(const CellType& cellType) const
      {return cellType_;}

      /** The ID number of the mesh in which these cells exist */
      int meshID() const {return mesh_.id();}

      /** The dimension of the cells contained in this set */
      int dimension() const {return dim_;}
      
      /** The mesh in which these cells exist */
      const Mesh& mesh() const {return mesh_;}

      /** The type of the cells contained in this set */
      const CellType& cellType() const {return cellType_;}

    private:

      /** the mesh in which the set exists */
      Mesh mesh_;

      /** the type of cell in the set */
      CellType cellType_;

      /** the dimension of the cells in the set */
      int dim_;
      
    };
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
