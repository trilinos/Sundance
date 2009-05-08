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

#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureFamily.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;


unsigned int iPow(unsigned int base, unsigned int exponent)
{
  int rtn = 1;
  for (unsigned int i=0; i<exponent; i++)
    {
      rtn *= base;
    }
  return rtn;
}

namespace SundanceStdFwk {
void checkbasis( BasisFamily &b1 , BasisFamily &b2 )
{
  unsigned int maxDim=3;
  double tol = 1.0e-13;
  int maxDiffOrder = 0;
  int numErrors = 0;
  QuadratureFamily quad = new GaussianQuadrature(4);
  
  for (unsigned int spatialDim=1; spatialDim<=maxDim; spatialDim++) {
    cerr << "\t" << "spatial dimension =" << spatialDim << endl;
    for (unsigned int cellDim=0; cellDim<=spatialDim; cellDim++) { 
      cerr << "\t\t" << "cell dimension =" << cellDim << endl;
      CellType cellType;
      if (cellDim==0) cellType=PointCell;
      if (cellDim==1) cellType=LineCell;
      if (cellDim==2) cellType=TriangleCell;
      if (cellDim==3) cellType=TetCell;
      
      Array<Point> qPts;
      Array<double> qWts;
      quad.getPoints(cellType, qPts, qWts);
      
      for (int d=0; d<=maxDiffOrder; d++) {
	if (cellDim==0 && d>0) continue;
	cerr << "\t\t\t" << "differentiation order = " << d << endl;
	for (unsigned int dir=0; dir<iPow(cellDim, d); dir++) {
	  cerr << "\t\t\t\t" << "direction = " << dir << endl;
	  MultiIndex mi;
	  mi[dir]=d;
	  Array<Array<double> > values1;
	  Array<Array<double> > values2;
	  cerr << "\t\t\t\t" << "computing basis1...";
	  b1.ptr()->refEval(spatialDim, cellType, qPts, mi, values1);
	  cerr << "done" << endl << "\t\t\t\t" << "computing basis2...";
	  b2.ptr()->refEval(spatialDim, cellType, qPts, mi, values2);
	  cerr << "done" << endl;
	  int nNodes1 = b1.ptr()->nNodes(spatialDim, cellType);
	  int nNodes2 = b2.ptr()->nNodes(spatialDim, cellType);
	  cerr << "\t\t\t\t" << "num nodes: basis1=" << nNodes1
	       << " basis2=" << nNodes2 << endl;
	  if (nNodes1 != nNodes2) {	
	    cerr << "******** ERROR: node counts should be equal" << endl;
	    numErrors++;
	    continue;
	  }
	  if (values1.size() != values2.size()) {
	    cerr << "******** ERROR: value array outer sizes should be equal" << endl;
	    numErrors++;
	    continue;
	  }
	  if (values1.size() != qPts.size()) {
	    cerr << "******** ERROR: value array outer size should be equal to number of quad points" << endl;
	    numErrors++;
	    continue;
	  }
	  for (int q=0; q<qPts.length(); q++) {
	    if (values1[q].length() != nNodes1) {
	      cerr << "******** ERROR: value array inner size should be equal to number of nodes" << endl;
	      numErrors++;
	      continue;
	    }
	    cerr << "\t\t\t\t\t" << "quad point q=" << q << " pt=" << qPts[q]
		 << endl;
	    for (int n=0; n<nNodes1; n++) {
	      cerr << "\t\t\t\t\t\t" << "node n=" << n << " phi1="
		   << values1[q][n]
		   << " phi2=" << values2[q][n]
		   << " |phi1-phi2|=" << fabs(values1[q][n]-values2[q][n])
		   << endl;
	      if (fabs(values1[q][n]-values2[q][n]) > tol) {
		cout << "ERROR" << endl; numErrors++;
	      }
	    }
	  }
	}
      }
    }
  }    
  cerr << endl << endl << "Summary: detected " << numErrors << " errors " << endl;
}
}
