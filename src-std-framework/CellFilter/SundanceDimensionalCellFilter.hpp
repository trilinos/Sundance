/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DIMENSIONALCELLFILTER_H
#define SUNDANCE_DIMENSIONALCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilterBase.hpp"


namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;
  using namespace Internal;


  /**
   * DimensionalCellFilter is a filter that identifies all cells of a 
   * specified dimension. 
   *
   * <h4> Example: </h4> get all faces in a 3D mesh
   *
   * \code
   * Mesh myMesh = myReader.getMesh();
   * CellFilter faceFilter = new DimensionalCellFilter(2);
   * CellSet faces = faceFilter.getCells(myMesh);
   * \endcode
   */
  class DimensionalCellFilter : public CellFilterBase 
  {
  public:
    /** */
    DimensionalCellFilter(int dim);

#ifndef DOXYGEN_DEVELOPER_ONLY

    /** */
    virtual ~DimensionalCellFilter(){;}

    /** */
    virtual int dimension(const Mesh& mesh) const {return dim_;}

    /** */
    virtual XMLObject toXML() const ;

    /** */
    virtual string typeName() const {return "DimensionalCellFilter";}

    /** */
    virtual bool lessThan(const CellFilterBase* other) const ;

    /* */
    GET_RCP(CellFilterBase);

  protected:
    /** get the cells */
    virtual CellSet internalGetCells(const Mesh& mesh) const ;
    
  private:
    
    int dim_;

#endif  /* DOXYGEN_DEVELOPER_ONLY */

  };

}



#endif
