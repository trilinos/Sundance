/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SIDECELLFILTER_H
#define SUNDANCE_SIDECELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilterBase.hpp"



namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal {}
  using namespace Internal;
  using namespace Teuchos;

  /** 
   *
   **/
  class SideCellFilter : public CellFilterBase 
  {
  public:
    /** Empty ctor */
    SideCellFilter();

#ifndef DOXYGEN_DEVELOPER_ONLY

    /** */
    virtual ~SideCellFilter(){;}

    /** Return the dimension of the cells that will be identified
     * by this filter when acting on the given mesh */
    virtual int dimension(const Mesh& mesh) const ;

    /** Write to XML */
    virtual XMLObject toXML() const {return XMLObject(typeName());}

    /** Return the type name */
    virtual string typeName() const {return "SideCellFilter";}

    /** Compare to another object */
    virtual bool lessThan(const CellFilterStub* other) const ;

    /* Handleable boilerplate */
    GET_RCP(CellFilterStub);
    

  protected:
    /** */
    virtual CellSet internalGetCells(const Mesh& mesh) const ;

#endif  /* DOXYGEN_DEVELOPER_ONLY */
  };

}



#endif
