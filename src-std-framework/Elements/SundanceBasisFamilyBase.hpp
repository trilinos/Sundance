/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
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

namespace SundanceStdFwk {

using namespace Teuchos;
using namespace TSFExtended;
using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;

namespace Internal {

/** Abstract interface for a basis family for defining discretized functions
 * on a mesh cell.
 *
 * This interface allows the definition of a basis representation of a
 * function on a single 1D, 2D, 3D (or higher) cell in a mesh.  Currently,
 * only an enumerated set of cell types are supported (see <tt>CellType</tt>).
 *
 * A function \f$g(x)\f$ defined on a cell is represented as:

 \f[
  g(x) = \sum_{i=0}^{N-1} \bar{g}_i \psi_i(x)
 \f]

 * where \f$x\in\Re^{D}\f$ is the spatial spatial coordinate for spatical dim
 * \f$D\f$ that lives in the cell's domain, \f$\psi_i(x)\f$ is the \f$i\f$th
 * basis function (of order = <tt>order()</tt>), \f$\bar{g}_i\f$ is the
 * (constant) coefficient for the \f$i\f$th basis function, and \f$N\f$ is the
 * number of basis functions.  Therefore, given the coefficients of the basis
 * function on the cell \f$\bar{g}_i\f$, one can compute the value of the
 * function \f$g(x)\f$ on the cell at any point \f$x\f$ in the domain of the
 * cell given the above summation formula.  This interface refers to the
 * coefficients \f$\bar{g}_i\f$ as degrees of freedom (<b>DOFs</b>).
 *
 * This interface allows the specification basis functions and basis
 * coefficients (i.e. DOFs) to be associated with any of the facets of a cell,
 * including the cell itself without restriction.  See the function
 * <tt>getLocalDOFs()</tt> for how this mapping of DOFs to facets for a single
 * function defined on the cell.
 *
 * It is important to note that each cell type, i.e. the values of the enum
 * <tt>CellType</tt>, has a "agreed upon" geometry for the "reference cell"
 * and this geometry must be known by the client of the basis family and the
 * concrete implementations of the basis family.
 *
 * \section Sundance_BasisFamily_ToDo_sec Possible Refactorings
 *
 * <ul>
 *
 * <li>Remove the constructor
 *
 * <li>Factor out a TopoBasisFamilyBase interface that does not include the
 * refEval() function.  This allows us to get around the problem of the
 * definition of the "reference cell".
 *
 * <li>We need to define a query function that determines if particular basis
 * family implementation supports all cell types.  We should not expect every
 * concrete basis family implementation to support every cell type do we?
 * Therefore, we need to add a bool function like
 * <tt>supportsCellType(cellType)</tt> that returns true if the cell type
 * supported and false otherwise.
 *
 * </ul>
 */
class BasisFamilyBase
  : public Handleable<BasisFamilyBase>,
    public TSFExtended::Printable,
    public TSFExtended::ObjectWithVerbosity<BasisFamilyBase>
{
public:

  /** \brief .
   *
   * ToDo: Remove this constructor for the base interface.
   */
  BasisFamilyBase();

  /** \brief .
   *
   * ToDo: Remove this function since it is not needed.
   */
  virtual ~BasisFamilyBase(){;}

  /** \brief Get a description of the DOF numbering and distribution scheme
   * for this basis function on the given cell type.
   *
   * \param  cellType
   *           [in] Specification of the cell topology
   * \param  dofs
   *           [out] Array of dof numbering information, to be filled in
   *           during the call.  On output,
   *           <tt>dofs.size()==dimension(cellType)</tt>.  See description of
   *           <tt>dofs</tt> below for more details.
   *
   * The DOF description is returned in the nested array <tt>dofs</tt>, and is
   * to be interpreted as follows: The outer dimension of the description
   * array <tt>dofs.size()</tt> is <tt>cellDim</tt>, where <tt>cellDim</tt> is
   * the spatial dimension of the cell.  The DOFs attached to facets are
   * stored in array entries <tt>dofs[s]</tt> where <tt>s=0...cellDim-1</tt>.
   * Those associated with the cell body are stored in
   * <tt>dofs[cellDim-1]</tt>. For cell dofs attached to facets, the dof
   * <tt>facetDofIndex</tt> associated with the facet <tt>facetIndex</tt> of
   * facet dimension <tt>facetDim</tt> is given by:

   \code
     dofs[facetDim][facetIndex][faceDofIndex] 
   \endcode
   
   * For dofs attached to the cell body, the local DOF within the entire cell
   * is given by dof is given by

   \code
   dofs[cellDim][0][dofIndex]
   \endcode 

   * More specifically:<ul>
   *
   * <li><tt>dof[facetDim].size()</tt> gives the number of facets of the facet
   * dimension <tt>facetDim</tt>, where <tt>0 <= facetDim <= cellDim</tt>
   *
   * <li><tt>dof[facetDim][facetIndex].size()</tt> gives the number of degrees
   * of freedom (DOFs) on the facet <tt>facetIndex</tt> with facet dimension
   * <tt>facetDim</tt>, where <tt>0 <= facetDim <= cellDim</tt> and <tt>0 <=
   * facetIndex < numFacets(cellType,facetDim)</tt>.
   *
   * </ul>
   *
   * For example, the Lagrange basis functions of order 0 through 3 on 2D
   * triangles would have the following dof arrays:

   \verbatim

    Order 0:

    { {}, {}, {{0}} }
   
    Order 1:

    { { {0}, {1}, {2} }, {}, {} }
    
    Order 2:

    { { {0}, {1}, {2} }, { {3}, {4}, {5} }, {} }
    
    Order 3:

    { { {0}, {1}, {2} }, { {3,4}, {5,6}, {7,8} }, {9} }

   \endverbatim

   * Above, we have used the ordering given in Hughes' textbook.
   */
  virtual void getLocalDOFs(
    const CellType& cellType,
    Array<Array<Array<int> > >& dofs
    ) const = 0 ;

  /** \brief Return the number of basis coefficients nodes associated with
   * this basis on the given cell type.
   *
   * ToDo: Change the name of this function to something like
   * <tt>nDOFsPerCell()</tt>.
   *
   * ToDo: Remove the spatialDim argument since it is not needed and not used
   * and the dimension of the cell is intrisic in the cellType.
   */
  virtual int nNodes(
    int spatialDim, const CellType& cellType
    ) const = 0;

  /** \brief Evaluate the basis functions (or some mixed spatial derivative of
   * the basis functions) for an array of points on the "reference cell" for a
   * given cell type.
   *
   * \param  spatialDim
   *           [in] Dimension of the embedding space
   * \param  cellType
   *           [in] The type of cell
   * \param  pts
   *           [in] Array of points on the reference cell (or master cell)
   *           where the basis functions are to be computed.
   * \param  deriv
   *           [in] Specifies which derivatives in each spatical direction
   *           to evaluate the basis functions derivatives for.
   * \param  result
   *           [out] On output, gives a doublely nested array which contain
   *           the basis functions (or the requested basis function
   *           derivatives) evaluated at the given points <tt>pts</tt>.  The
   *           size of the outer array <tt>results</tt> is \code
   *           results.size() == pts.size() *
   *           this->nNodes(spatialDim,cellType) \endcode.  More specifically,
   *           \code results[pointIndex][basisIndex] \endcode gives the value
   *           of the spatial derivative
   *           \f[\frac{\partial^{d_x+d_y+d_z}}{\partial x^{d_x} \partial
   *           y^{d_y} \partial z^{d_z}}\psi_i(x,y,z)\f] (where \f$d_x\f$ =
   *           <tt>deriv[0]</tt>, \f$d_y\f$ = <tt>deriv[1]</tt>, and \f$d_Z\f$
   *           = <tt>deriv[2]</tt>) at the point <tt>pointIndex</tt> (where
   *           <tt>0 <= pointIndex < pts.size()</tt>) for the basis function
   *           \f$i\f$ = <tt>basisIndex</tt> (where <tt>0 <= basisIndex <
   *           mapStructure.numBasisChunks()</tt>).  If the cell is just a 2D
   *           cell then <tt>deriv[2]</tt> is ignored and if it is just a 1D
   *           cell then <tt>deriv[1]</tt> and <tt>deriv[2]</tt> are ignored.
   *           Obviously, if <tt>deriv[k]==0</tt> for <tt>k=0,1,2</tt> then
   *           just the value of the basis function at the points is evaluated
   *           and returned.
   *
   * ToDo: Remove the spatialDim argument since it is not needed and not used
   * and the dimension of the cell is intrisic in cellType.
   */
  virtual void refEval(int spatialDim,
                       const CellType& cellType,
                       const Array<Point>& pts,
                       const MultiIndex& deriv,
                       Array<Array<double> >& result) const = 0 ;


  /** \brief Return the polynomial order of the basis.
   *
   * The order of the basis functions determines how many non-zero spatial
   * derivatives that each basis function has.  For example, a basis family of
   * order 2 gives a quadratic approximation of the function over the cell.
   */
  virtual int order() const = 0 ;

  /** \brief Return the spatial dimension of the basis.
   *
   * For scalar basis types <tt>dim()==1</tt>.  For vector basis functions the
   * dimension would typically be equal to the cell dimension.
   *
   * ToDo: I think this function needs to be changed to <tt>dim(cellType)</tt>
   * since the dimension of the basis will be determined by the spatical
   * dimension of the cell if this is indeed a vector basis.
   */
  virtual int dim() const = 0 ;

  /** \brief Ordering rule for basis families.
   *
   * ToDo: This function is not virtual and does not need to be in this
   * interface.  In fact, it is not clear to me how you can even implement
   * this function in general since you are really asking any arbitrary basis
   * family implementation to know how to compare itself to any other basis
   * family implementation.  Suggestion: Why not just order the basis
   * functions the combination of <tt>typeid(*this)</tt>,
   * <tt>this->order()</tt>, and <tt>this->dim()</tt>.  This should give a unique
   * way to order basis families.
   */
  bool lessThan(const BasisFamilyBase* other) const 
    {return order() < other->order();}

};

} // namespace Internal

} // namespace SundanceStdFwk

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
