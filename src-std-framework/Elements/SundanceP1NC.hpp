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

#ifndef SUNDANCE_P1NC_H
#define SUNDANCE_P1NC_H

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

  /** 
   * P1 Nonconforming basis 
   */
  class P1NC : public ScalarBasis
  {
  public:
    /** */
    P1NC(){;}

    /** */
    virtual ~P1NC(){;}

    /** */
    virtual void print(ostream& os) const ;

    /** */
    virtual int order() const {return 1;}

    /** return the number of nodes for this basis on the given cell type */
    virtual int nNodes(int spatialDim, 
                       const CellType& cellType) const ;

    /** */
    virtual void getLocalDOFs(const CellType& cellType,
                              Array<Array<Array<int> > >& dofs) const ;

    /** */
    virtual void refEval(int spatialDim, 
                         const CellType& cellType,
                         const Array<Point>& pts,
                         const MultiIndex& deriv,
                         Array<Array<double> >& result) const ;

    /* Handleable boilerplate */
    GET_RCP(BasisFamilyBase);

  private:
    /** evaluate on a line cell  */
    void evalOnLine(const Point& pt,
                    const MultiIndex& deriv,
                    Array<double>& result) const ;
    
    /** evaluate on a triangle cell  */
    void evalOnTriangle(const Point& pt,
                        const MultiIndex& deriv,
                        Array<double>& result) const ;
  };
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
