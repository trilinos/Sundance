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

/** 
 * Exo2Triangle converts a NCDF mesh to a Triangle format mesh
 */


int main(int argc, void** argv)
{
  
  try
		{
      string inputFile;
      string outputFile;
      int indexOffset=0;

      Sundance::setOption("i", inputFile, "Input mesh file");
      Sundance::setOption("o", outputFile, "Output mesh file");
      Sundance::setOption("offset", indexOffset, "Index offset");

      Sundance::init(&argc, &argv);

      MeshType meshType = new BasicSimplicialMeshType();

      /* Read the input mesh */
      MeshSource mesher 
        = new ExodusNetCDFMeshReader(inputFile, meshType);
      Mesh mesh = mesher.getMesh();


      FieldWriter w = new TriangleWriter(outputFile, indexOffset);
      w.addMesh(mesh);
      w.write();


    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  TimeMonitor::summarize();
  Sundance::finalize();
}
