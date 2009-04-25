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

#ifndef SUNDANCE_BASISREFERENCEEVALUATIONBASE_H
#define SUNDANCE_BASISREFERENCEEVALUATIONBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellType.hpp"
#include "Teuchos_Array.hpp"


namespace SundanceUtils {class Point;}
namespace SundanceCore {namespace Internal {class MultiIndex;}}

namespace SundanceStdFwk {


using Teuchos::Array;
using SundanceUtils::Point;
using SundanceStdMesh::CellType;
using SundanceCore::Internal::MultiIndex;


/** 
 * Abstract interface for evaluation of basis functions and their
 * spatial derivatives on reference cells.
 */
class BasisReferenceEvaluationBase
{
public:

  /** \brief Evaluate the basis functions (or some mixed spatial derivative of
   * the basis functions) for an array of points on the "reference cell" for a
   * given cell type.
   *
   * \param  maximalCellType
   *           [in] Cell type for the maximal-dimension reference cell
   *           in which the current cell is embedded. This information 
   *           is required for certain basis families where basis members
   *           from throughout the cell have nonzero support on facets
   *           to which their associated DOFs are not directly attached. 
   *           If the current cell
   *           is already of maximal dimension, the value of 
   *           <tt>maximalCellType</tt> is simply equal to <tt>cellType.</tt>
   * \param  cellType
   *           [in] The type of cell on which the basis is currently being
   *           evaluated. 
   *           Precondition: <tt>dimension(cellType) 
   *           <= dimension(maximalCellType).</tt>
   * \param  pts
   *           [in] Array of points on the reference cell (or master cell)
   *           where the basis functions are to be computed. 
   * \param  deriv
   *           [in] Specifies which derivatives in each spatial direction
   *           to evaluate the basis functions derivatives for. 
   *           Preconditions: <tt>deriv[d]==0</tt> for any <tt>d</tt> such
   *           that <tt>dimension(maximalCellType) <= d.</tt>
   * \param  result
   *           [out] On output, gives a triply nested array which contain
   *           the basis functions (or the requested basis function
   *           derivatives) evaluated at the given points <tt>pts</tt>.  The
   *           size of the outer array <tt>results</tt> is either zero
   *           or spatialDim, depending on whether this is a scalar or
   *           vector basis, respectively. The size of the next
   *           array level is equal to the number of evaluation points. 
   *           Finally, the size of the innermost array level is equal to
   *           the number of DOFs visible from the given cell type.
x   *           Specifically,
   *           \code 
   *           results[k][pointIndex][basisIndex] 
   *           \endcode gives the value
   *           of the spatial derivative of the \f$k\f$-th component of
   *           \f[\frac{\partial^{d_x+d_y+d_z}}{\partial x^{d_x} \partial
   *           y^{d_y} \partial z^{d_z}}\psi_i(x,y,z)\f],
   *           where \f$d_x\f$ =
   *           <tt>deriv[0]</tt>, \f$d_y\f$ = <tt>deriv[1]</tt> (in 2D or 3D)
   *           and \f$d_Z\f$
   *           = <tt>deriv[2]</tt> (in 3D) at the point <tt>pointIndex</tt> 
   *           (where
   *           <tt>0 <= pointIndex < pts.size()</tt>) for the basis function
   *           \f$i\f$ = <tt>basisIndex</tt> (where <tt>0 <= basisIndex <
   *           mapStructure.numBasisChunks()</tt>). 
   */
  virtual void refEval(
    const CellType& maximalCellType,
    const CellType& cellType,
    const Array<Point>& pts,
    const MultiIndex& deriv,
    Array<Array<Array<double> > >& result
    ) const = 0 ;  

};

}


#endif
