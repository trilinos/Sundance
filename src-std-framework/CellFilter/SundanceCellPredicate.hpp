/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLPREDICATE_H
#define SUNDANCE_CELLPREDICATE_H



#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellPredicateBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;
  using namespace Internal;
  
    /** 
     * User-level handle for predicates (deriving from CellPredicateBase)
     * used to decide whether
     * a given cell passes a CellFilter.
     */
  class CellPredicate : public TSFExtended::Handle<CellPredicateBase>
  {
    public:
    
    /* handle boilerplate */
    HANDLE_CTORS(CellPredicate, CellPredicateBase);

    /** write to XML */
    XMLObject toXML() const {return ptr()->toXML();}

#ifndef DOXYGEN_DEVELOPER_ONLY

    /** set the mesh on which cells are to be tested */
    void setMesh(const Mesh& mesh, int cellDim) const 
    {ptr()->setMesh(mesh, cellDim);}

    /** compare to another predicate, used for placement in STL containers */
    bool operator<(const CellPredicate& other) const ;


#endif  /* DOXYGEN_DEVELOPER_ONLY */    
    };

}

namespace std
{
  inline ostream& operator<<(ostream& os, const SundanceStdFwk::CellPredicate& pred)
  {
    os << pred.toXML() << endl;
    return os;
  }
}


#endif
