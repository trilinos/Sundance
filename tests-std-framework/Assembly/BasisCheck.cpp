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

#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
#include "SundanceFIATLagrange.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceEssentialBC.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceAssembler.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceRefIntegral.hpp"
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

unsigned int iPow(unsigned int base, unsigned int exponent)
{
  int rtn = 1;
  for (unsigned int i=0; i<exponent; i++)
  {
    rtn *= base;
  }
  return rtn;
}

int main(int argc, char** argv)
{
  int stat = 0;

  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());
    Tabs tab0;

    unsigned int pMax = 3;
    unsigned int maxDim=2;
    double tol = 1.0e-13;
    int maxDiffOrder = 0;
    int numErrors = 0;
      
    QuadratureFamily quad = new GaussianQuadrature(4);

    for (unsigned int p=1; p<=pMax; p++)
    {
      Tabs tab1;
      cerr << tab1 << "Polynomial order p=" << p << endl;
      for (unsigned int spatialDim=0; spatialDim<=maxDim; spatialDim++)
      {
        Tabs tab2;
        cerr << tab2 << "spatial dimension =" << spatialDim << endl;
        for (unsigned int cellDim=0; cellDim<=spatialDim; cellDim++)
        {   
          Tabs tab3;
          cerr << tab3 << "cell dimension =" << cellDim << endl;
          CellType cellType;
          if (cellDim==0) cellType=PointCell;
          if (cellDim==1) cellType=LineCell;
          if (cellDim==2) cellType=TriangleCell;
          if (cellDim==3) cellType=TetCell;
                  
          Array<Point> qPts;
          Array<double> qWts;
          quad.getPoints(cellType, qPts, qWts);

          BasisFamily b1 = new Lagrange(p);
#ifdef HAVE_FIAT
          BasisFamily b2 = new FIATLagrange(p);
#else
          BasisFamily b2 = new Lagrange(p);
#endif
          for (int d=0; d<=maxDiffOrder; d++)
          {
            if (cellDim==0 && d>0) continue;
            Tabs tab4;
            cerr << tab4 << "differentiation order = " << d << endl;
                      
            for (unsigned int dir=0; dir<iPow(cellDim, d); dir++)
            {
              Tabs tab5;
              cerr << tab5 << "direction = " << dir << endl;
              MultiIndex mi;
              mi[dir]=d;
              Array<Array<Array<double> > > values1;
              Array<Array<Array<double> > > values2;
              cerr << tab5 << "computing basis1...";
              b1.ptr()->refEval(cellType, cellType, qPts, mi, values1);
              cerr << "done" << endl;
              cerr << tab5 << "computing basis2...";
              b2.ptr()->refEval(cellType, cellType, qPts, mi, values2);
              cerr << "done" << endl;
              int nNodes1 = b1.ptr()->nReferenceDOFs(cellType, cellType);
              int nNodes2 = b2.ptr()->nReferenceDOFs(cellType, cellType);
              cerr << tab5 << "num nodes: basis1=" << nNodes1
                   << " basis2=" << nNodes2 << endl;
              if (nNodes1 != nNodes2) 
              {
                cerr << "******** ERROR: node counts should be equal" << endl;
                numErrors++;
                continue;
              }
              if (values1.size() != values2.size())
              {
                cerr << "******** ERROR: value array outer sizes should be equal" << endl;
                numErrors++;
                continue;
              }
              if (values1[0].size() != qPts.size())
              {
                cerr << "******** ERROR: value array outer size should be equal to number of quad points" << endl;
                numErrors++;
                continue;
              }
              for (int q=0; q<qPts.length(); q++)
              {
                if (values1[0][q].length() != nNodes1)
                {
                  cerr << "******** ERROR: value array inner size should be equal to number of nodes" << endl;
                  numErrors++;
                  continue;
                }
                Tabs tab6;
                cerr << tab6 << "quad point q=" << q << " pt=" << qPts[q]
                     << endl;
                for (int n=0; n<nNodes1; n++)
                {
                  Tabs tab7;
                  cerr << tab7 << "node n=" << n << " phi1=" 
                       << values1[0][q][n] 
                       << " phi2=" << values2[0][q][n] 
                       << " |phi1-phi2|=" << fabs(values1[0][q][n]-values2[0][q][n]) 
                       << endl;
                  if (fabs(values1[0][q][n]-values2[0][q][n]) > tol) { cout << "ERROR" << endl; numErrors++; }
                }
              }
            }
          }
        }
      }
    }

    cerr << endl << endl << "Summary: detected " << numErrors << " errors " << endl;
      
    if (numErrors == 0)
    {
      cerr << "BasisCheck PASSED" << endl;
    }
    else
    {
      cerr << "BasisCheck FAILED" << endl;
      stat = -1;
    }
    TimeMonitor::summarize();
  }
	catch(exception& e)
  {
    cerr << e.what() << endl;
    cerr << "BasisCheck FAILED" << endl;
    stat = -1;
  }

  return stat;
}
