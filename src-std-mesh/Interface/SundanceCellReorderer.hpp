/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLREORDERER_H
#define SUNDANCE_CELLREORDERER_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellReordererBase.hpp"
#include "TSFHandle.hpp"


namespace Sundance
{
  using namespace Teuchos;
  using namespace Internal;
  using namespace TSFExtended;
  class Mesh;
  /**
   * User-level handle class for abstract specification of
   * cell reordering algorithms.
   *
   * <h4> Examples </h4>
   *
   * \code
   * CellReorderer bfs = new BreadthFirstReorderer();
   * mesh.setReorderer(bfs);
   * \endcode
   */
  class CellReorderer
    : public TSFExtended::Handle<CellReordererFactoryBase>
  {
  public:
    HANDLE_CTORS(CellReorderer, CellReordererFactoryBase);
      
    /** */
    RefCountPtr<CellReordererImplemBase> 
    createInstance(const MeshBase* mesh) const ;
  };
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
