/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_BASISFAMILYBASE_H
#define SUNDANCE_BASISFAMILYBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellType.hpp"
#include "SundancePoint.hpp"
#include "SundanceMultiIndex.hpp"
#include "TSFHandleable.hpp"
#include "TSFPrintable.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "Teuchos_XMLObject.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace Teuchos;
  using namespace TSFExtended;
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  
  namespace Internal
  {
    /** 
     * 
     */
    class BasisFamilyBase : public Handleable<BasisFamilyBase>,
                            public TSFExtended::Printable,
                            public TSFExtended::ObjectWithVerbosity<BasisFamilyBase>
    {
    public:
      /** */
      BasisFamilyBase();

      /** */
      virtual ~BasisFamilyBase(){;}

      /** 
       * Get a description of the DOF numbering scheme for this
       * basis function on the given cell type. The DOF description
       * is returned by reference argument, and is to be interpreted
       * as follows: The outer dimension of the description array 
       * is D, where D is the spatial dimension of the cell.
       * The DOFs attached to
       * facets are stored in array elements 0 through D-1, those
       * associated with the cell body are stored in array element D. 
       * For dofs attached to facets, the 
       * <tt>dofIndex'</tt>th dof associated
       * with the <tt>facetIndex'</tt>th facet of dimension  <tt>dim</tt>
       * is given by the array entry:
       * \code
       * dofs[facetDim][cellIndex][dofIndex]; 
       * \endcode
       * For dofs attached to the cell body, the <tt>dofIndex'</tt>
       * dof is given by 
       * \code
       * dofs[cellDim][0][dofIndex];
       * \endcode 
       * For example, the Lagrange basis functions of order 0 through
       * 3 on triangles would have the following dof arrays:
       * \code
       * Order 0:
       * { {}, {}, {{0}} }
       *
       * Order 1:
       * { {{0}, {1}, {2}}, {}, {} }
       * 
       * Order 2:
       * { {{0}, {1}, {2}}, {{3}, {4}, {5}}, {} }
       * 
       * Order 3:
       * { {{0}, {1}, {2}}, {{3,4}, {5,6}, {7,8}}, {9} }
       * \endcode
       * We have used the ordering given in Hughes' textbook.
       *
       * @param cellType specification of the cell topology
       * @param dofs array of dof numbering information, to be filled
       * in during the call
       */
      virtual void getLocalDOFs(const CellType& cellType,
                                Array<Array<Array<int> > >& dofs) const = 0 ;

      /** Return the number of nodes associated with this basis on 
       * the given cell type */
      virtual int nNodes(const CellType& cellType) const = 0 ;

      /** 
       * Evaluate the basis function on an array of points on the 
       * reference element for the given cell type.
       *
       * @param cellType specification of the type of cell on which
       * to evaluate the basis function in this call
       * @param pts array of quadrature points at which the basis function
       * is to be evaluated
       * @param deriv multiindex specifying the derivative to be evaluated
       * @param result vector of basis function values
       */
      virtual void refEval(const CellType& cellType,
                           const Array<Point>& pts,
                           const MultiIndex& deriv,
                           Array<Array<double> >& result) const = 0 ;


      /** Return the polynomial order of the basis */
      virtual int order() const = 0 ;

      /** Return the spatial dimension of the basis */
      virtual int dim() const = 0 ;

      /** Ordering for storage in STL maps. No need to check type equality 
       * here, because it's already done at the handle level. */
      bool lessThan(const BasisFamilyBase* other) const 
      {return order() < other->order();}
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
