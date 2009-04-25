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
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceRaviartThomas.hpp"
#include "SundanceFunctionSupportResolver.hpp"


int main(int argc, char** argv)
{
  int stat = 0;
  try
  {
    Sundance::init(&argc, &argv);

    /* Create a mesh. It will be of type BasisSimplicialMesh, and will
     * be built using a PartitionedLineMesher. */
    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedRectangleMesher(0.0, 10.0, 2, 1,
      0.0, 10.0, 2, 1,
      meshType);
    Mesh mesh = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();

    /* */
    Expr v = new TestFunction(new RaviartThomas(2));
    /* */
    Expr u = new UnknownFunction(new RaviartThomas(2));

    Out::os() << "v = " << describeFunction(v) << endl;
    Out::os() << "u = " << describeFunction(u) << endl;



    QuadratureFamily quad = new GaussianQuadrature(2);
    Expr eqn = Integral(interior, v*u, quad);
    Expr dum;

    RCP<FunctionSupportResolver> fsr 
      = rcp(new FunctionSupportResolver(eqn, dum, tuple(v), tuple(u),
          dum, dum, tuple(dum), false, 4));
          
    
    DOFMapBuilder builder(mesh, fsr, false);

    for (unsigned int br = 0; br<builder.rowMap().size(); br++)
    {
      RCP<DOFMapBase> rm = builder.rowMap()[br];
      rm->print(Out::os());
    }
  }
	catch(exception& e)
  {
    stat = -1;
    cerr << "RT dof test FAILED" << endl;
    cerr << e.what() << endl;
  }

  return stat;
  
}
