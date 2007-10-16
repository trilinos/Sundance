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

#include "SundanceExodusWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_XMLObject.hpp"

#ifdef HAVE_EXODUS 
#include "exodusII.h"
#endif


using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;



void ExodusWriter::write() const 
{
#ifdef HAVE_EXODUS
  int ws = 8;
  int exoid = ex_create(filename().c_str(), EX_CLOBBER, &ws, &ws);

  TEST_FOR_EXCEPTION(exoid < 0, RuntimeError, "failure to create file "
    << filename());

  writeMesh(exoid);

  if (nProc() > 1 && myRank()==0) 
  {
    // do parallel header stuff
  }

  ex_close(exoid);
#else
  TEST_FOR_EXCEPTION(true, RuntimeError, "Exodus not enabled");
#endif
}



void ExodusWriter::writeMesh(int exoid) const
{
#ifdef HAVE_EXODUS

  int ierr = 0;

  int dim = mesh().spatialDim();
  int nElems = mesh().numCells(dim);

  ierr = ex_put_init(
    exoid, 
    filename().c_str(), 
    dim,
    mesh().numCells(0), 
    nElems,
    mesh().numLabels(dim), 
    mesh().numLabels(0), 
    mesh().numLabels(dim-1)
    );

  char* qa_record[1][4];
  qa_record[0][0] = "Sundance";
  qa_record[0][1] = "sundance";
  qa_record[0][2] = "date";
  qa_record[0][3] = "time";

  ierr = ex_put_qa(exoid, 1, qa_record);

  int nPts = mesh().numCells(0);
  Array<double> x(nPts);
  Array<double> y(nPts);
  Array<double> z(nPts * (dim > 2));

  for (int n=0; n<nPts; n++)
  {
    const Point& P = mesh().nodePosition(n);
    x[n] = P[0];
    y[n] = P[1];
    if (dim==3) z[n] = P[2];
  }

  if (dim==2)
  {
    ierr = ex_put_coord(exoid, &(x[0]), &(y[0]), (void*) 0);
  }
  else
  {
    ierr = ex_put_coord(exoid, &(x[0]), &(y[0]), &(z[0]));
  }

  if (dim==2)
  {
    char* coordNames[2];
    coordNames[0] = "x";
    coordNames[1] = "y";
    ierr = ex_put_coord_names(exoid, coordNames);
  }
  else
  {
    char* coordNames[3];
    coordNames[0] = "x";
    coordNames[1] = "y";
    coordNames[2] = "z";
    ierr = ex_put_coord_names(exoid, coordNames);
  }


  /* write the element blocks */
  Array<int> blockLabels = mesh().getAllLabelsForDimension(dim).elements();
  int nodesPerElem = dim+1;
  std::string eType = elemType(mesh().cellType(dim));
  for (unsigned int b=0; b<blockLabels.size(); b++)
  {
    int numBlockAttr = 0;
    Array<int> blockElemLIDs;
    Array<int> nodeLIDs;
    Array<int> orient;
    mesh().getLIDsForLabel(dim, blockLabels[b], blockElemLIDs);
    int numElemsThisBlock = blockElemLIDs.size();
    mesh().getFacetLIDs(dim, blockElemLIDs, 0, nodeLIDs, orient);
    ierr = ex_put_elem_block(
      exoid, blockLabels[b], eType.c_str(), 
      numElemsThisBlock, nodesPerElem, numBlockAttr
      );

    ierr = ex_put_elem_conn(exoid, blockLabels[b], &(nodeLIDs[0]));
  }

  
  /* write the side sets */
  Array<int> ssLabels = mesh().getAllLabelsForDimension(dim-1).elements();

  for (unsigned int ss=0; ss<ssLabels.size(); ss++)
  {
    Array<int> sideLIDs;
    Array<int> elemLIDs;
    Array<int> facets;

    mesh().getLIDsForLabel(dim-1, ssLabels[ss], sideLIDs);

    int numSides = sideLIDs.size();
    int numDists = 0;

    ierr = ex_put_side_set_param(exoid, ssLabels[ss], numSides, numDists);
    
    mesh().getMaxCofacetLIDs(dim-1, sideLIDs, elemLIDs, facets);

    ierr = ex_put_side_set(exoid, ssLabels[ss], &(elemLIDs[0]), &(facets[0]));
  }


  
  /* write the node sets */
  Array<int> nsLabels = mesh().getAllLabelsForDimension(0).elements();

  for (unsigned int ns=0; ns<nsLabels.size(); ns++)
  {
    Array<int> nodeLIDs;

    mesh().getLIDsForLabel(0, nsLabels[ns], nodeLIDs);

    int numNodes = nodeLIDs.size();
    int numDists = 0;

    ierr = ex_put_node_set_param(exoid, nsLabels[ns], numNodes, numDists);
    
    ierr = ex_put_node_set(exoid, nsLabels[ns], &(nodeLIDs[0]));
  }
#else
  TEST_FOR_EXCEPTION(true, RuntimeError, "Exodus not enabled");
#endif
}


std::string ExodusWriter::elemType(const CellType& type) const
{
  switch(type)
  {
    case TriangleCell:
      return "TRIANGLE";
    case TetCell:
      return "TETRA";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError, "cell type=" << type << " cannot be used as a "
        "maximal-dimension cell in exodus");
  }
  return "NULL"; //-Wall
}

