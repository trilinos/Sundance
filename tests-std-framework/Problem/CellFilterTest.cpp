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

#include "Sundance.hpp"
#include "SundanceCFMeshPair.hpp"
#include "SundanceDOFMapBuilder.hpp"


/** 
 * Tests logical operations on cell filters
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})

CELL_PREDICATE(ATest, {return x[0] >= 0.1 && x[0] <= 0.2;})

CELL_PREDICATE(BTest, {return x[0] >= 0.4 && x[0] <= 0.6;})

CELL_PREDICATE(CTest, {return x[0] >= 0.5 && x[0] <= 0.8;})

CELL_PREDICATE(DTest, {return x[0] >= 0.6 && x[0] <= 0.75;})

CELL_PREDICATE(ETest, {return x[0] >= 0.1 && x[0] <= 0.3;})

int main(int argc, void** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 20;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();
      for (int i=nx/3; i<nx/2; i++)
        {
          mesh.setLabel(1, i, 1);
        }
      for (int i=nx/2; i<2*nx/3; i++)
        {
          mesh.setLabel(1, i, 2);
        }
      for (int i=nx/7; i<nx/3; i++)
        {
          mesh.setLabel(1, i, 3);
        }

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter A = interior.subset(new ATest());
      CellFilter B = interior.subset(new BTest());
      CellFilter C = interior.labeledSubset(3);
      CellFilter D = interior.labeledSubset(1);
      CellFilter E = interior.labeledSubset(2);

      Set<int> fI = makeSet(0,1,2);
      Set<int> fA = makeSet(1);
      Set<int> fB = makeSet(2);
      Set<int> fC = makeSet(2);
      Set<int> fD = makeSet(1);
      Set<int> fE = makeSet(2);
      Set<int> fP = makeSet(0, 4);

      SundanceUtils::Map<CellFilter, Set<int> > fmap;
      fmap.put(interior, fI);
      fmap.put(B, fB);
      fmap.put(C, fC);
      fmap.put(D, fD);
      fmap.put(E, fE);
      fmap.put(leftPoint, fP);
      SundanceUtils::Map<CellFilter, SundanceUtils::Map<Set<int>, CellSet> > inputChildren;

      Array<SundanceUtils::Map<Set<int>, CellFilter> >djf = DOFMapBuilder::funcDomains(mesh, fmap, inputChildren);

      cout << "disjoint filters: " << endl;
      if (mesh.numCells(1) < 1000)
        {
          for (int d=0; d<=mesh.spatialDim(); d++)
            {
              cout << "dimension = " << d << endl;
              for (SundanceUtils::Map<Set<int>, CellFilter>::const_iterator
                     i=djf[d].begin(); i!=djf[d].end(); i++)
                {
                  cout << "cells = " << i->second.getCells(mesh)
                       << ", funcs=" << i->first << endl << endl;
                }
            }

          for (SundanceUtils::Map<CellFilter,SundanceUtils::Map<Set<int>, CellSet> >::const_iterator
                 i=inputChildren.begin(); i!=inputChildren.end(); i++)
            {
              cout << "input filter = " << i->first
                   << ", children = " << i->second << endl;
            }
        }

      
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
