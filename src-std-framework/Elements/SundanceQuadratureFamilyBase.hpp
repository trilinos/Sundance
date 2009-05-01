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

#ifndef SUNDANCE_QUADRATUREFAMILYBASE_H
#define SUNDANCE_QUADRATUREFAMILYBASE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellType.hpp"
#include "SundancePoint.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace Teuchos;

  namespace Internal
  {
    /** 
     * QuadratureFamilyBase extends QuadratureFamilyStub to provide
     * an interface for getting quadrature points for a given cell type. 
     */
    class QuadratureFamilyBase : public QuadratureFamilyStub 
    {
    public:
      /** */
      QuadratureFamilyBase(int order) : QuadratureFamilyStub(order) {;}

      /** */
      virtual ~QuadratureFamilyBase(){;}

      /** Gets number of points associated with a particular cell type:
          WARNING: this is slow.  Call it once and store the result. 
          TODO: make it pure virtual and override with queries in
          the derived classes, making them supply the information.  */
      virtual int getNumPoints( const CellType &cellType ) const 
      {
	Array<Point> qp;
	Array<double> qw;
	this->getPoints(cellType,qp,qw);
	return qp.size();
      }

      /** Get the quadrature points and weights for the given cell type */
      virtual void getPoints(const CellType& cellType, 
                     Array<Point>& quadPoints,
                     Array<double>& quadWeights) const ;
      
    protected:

      /** compute a rule for the reference line cell */
      virtual void getLineRule(Array<Point>& quadPoints,
                               Array<double>& quadWeights) const ;

      /** compute a rule for the reference triangle cell */
      virtual void getTriangleRule(Array<Point>& quadPoints,
                                   Array<double>& quadWeights) const ;

      /** compute a rule for the reference quad cell */
      virtual void getQuadRule(Array<Point>& quadPoints,
                               Array<double>& quadWeights) const ;

      /** compute a rule for the reference tet cell */
      virtual void getTetRule(Array<Point>& quadPoints,
                              Array<double>& quadWeights) const ;


      /** compute a rule for the reference brick cell */
      virtual void getBrickRule(Array<Point>& quadPoints,
                                Array<double>& quadWeights) const ;
    private:
    };
  }
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif

