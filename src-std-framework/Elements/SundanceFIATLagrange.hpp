/* @HEADER@ */
/* @HEADER@ */


// <-------- CHANGED INCLUDE GUARDS TO FIATLAGRANGE_H
#ifndef SUNDANCE_FIATLAGRANGE_H
#define SUNDANCE_FIATLAGRANGE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceScalarBasis.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace Teuchos;

  /** 
   * Lagrange basis 
   */
  class FIATLagrange : public ScalarBasis
  {
  public:
    /** */
    FIATLagrange(int order);

    /** */
    virtual ~FIATLagrange(){;}

    /** */
    virtual void print(ostream& os) const ;

    /** */
    virtual int order() const {return order_;}

    /** return the number of nodes for this basis on the given cell type */
    virtual int nNodes(const CellType& cellType) const ;

    /** */
    virtual void getLocalDOFs(const CellType& cellType,
                              Array<Array<Array<int> > >& dofs) const ;

    /** */
    virtual void refEval(const CellType& cellType,
                         const Array<Point>& pts,
                         const MultiIndex& deriv,
                         Array<Array<double> >& result) const ;

    /* Handleable boilerplate */
    GET_RCP(BasisFamilyBase);

  private:
    int order_;

    static const int MAXDEGREE_ = 8;

    static Array<int> makeRange(int low, int high);

    /** evaluate on a triangle cell  */
    void evalOnTriangle(const Point& pt,
                        const MultiIndex& deriv,
                        Array<double>& result) const ;
    
    static const int MAX_DEGREE_ = 1;
    Array<Array<double> > VDM_;
    Array<Array<Array<double> > > derivMats_;

  };
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif

