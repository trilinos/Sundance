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
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceFieldWriter.hpp"
#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceTestFunction.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceLagrange.hpp"
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
#include "TSFVectorType.hpp"
#include "TSFEpetraVectorType.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

bool leftPointTest(const Point& x) {return fabs(x[0]) < 1.0e-10;}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      MeshType meshType = new BasicSimplicialMeshType();

      int np = MPIComm::world().getNProc();

      

      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, 2*np, np,
                                                         0.0, 1.0, 2, 1,
                                                         meshType);

      Mesh mesh = mesher.getMesh(); 
      FieldWriter w = new VerboseFieldWriter();
      w.addMesh(mesh);
      w.write();

      Array<int> funcs = tuple(0);

      //verbosity<CellSet>() = VerbExtreme;

      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellPredicate leftEdgeFunc = new PositionalCellPredicate(leftPointTest);
      CellFilter leftEdge = edges.subset(leftEdgeFunc);
      
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      //      Expr u0 = new DiscreteFunction(new Lagrange(1), "u0");
      Expr u0 = new ZeroExpr();
      Expr v = new TestFunction(new Lagrange(1), "v");
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      

      QuadratureFamily quad = new GaussianQuadrature(2);
      Expr eqn = Integral(interior, (dx*v)*(dx*u) + (1.0+x*x)*v, quad);
      Expr bc = EssentialBC(leftEdge, v*u, quad);

      RefCountPtr<EquationSet> eqnSet 
        = rcp(new EquationSet(eqn, bc, v, u, u0));

      EquationSet::classVerbosity() = VerbHigh;
      HomogeneousDOFMap::classVerbosity() = VerbExtreme;
      Expr::showAllParens() = true;

      VectorType<double> vecType = new EpetraVectorType();

      Assembler assembler(mesh, eqnSet, vecType); 

      RefCountPtr<DOFMapBase> rowMap = assembler.rowMap();
      cerr << "DOF Map" << endl;
      for (int c=0; c<mesh.numCells(0); c++)
        {
          Array<int> dofs;
          rowMap->getDOFsForCell(0, c, 0, dofs);
          cerr << c << " " << dofs << endl;
        }
      cerr << "bc rows " << endl << *(assembler.bcRows()) << endl;
      
      LinearOperator<double> A;
      Vector<double> b;

    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
