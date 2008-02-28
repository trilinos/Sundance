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
using namespace std;


void ExodusWriter::write() const 
{
  string exoFile = filename();
  string parFile = filename();
  if (nProc() > 1) 
  {
    exoFile = exoFile + "-" + Teuchos::toString(nProc()) + "-" + Teuchos::toString(myRank());
    parFile = parFile + "-" + Teuchos::toString(nProc()) + "-" + Teuchos::toString(myRank());
  }
  exoFile = exoFile + ".exo";
  parFile = parFile + ".pxo";

  if (nProc() > 1) writeParallelInfo(parFile);
#ifdef HAVE_EXODUS
  int ws = 8;
  int exoid = ex_create(exoFile.c_str(), EX_CLOBBER, &ws, &ws);

  TEST_FOR_EXCEPTION(exoid < 0, RuntimeError, "failure to create file "
    << filename());

  writeMesh(exoid);

  ex_close(exoid);
#else
  TEST_FOR_EXCEPTION(true, RuntimeError, "Exodus not enabled");
#endif
}


void ExodusWriter::offset(Array<int>& x) const
{
  for (unsigned int i=0; i<x.size(); i++) x[i]++;
}


void ExodusWriter::writeMesh(int exoid) const
{
#ifdef HAVE_EXODUS

  int ierr = 0;

  int dim = mesh().spatialDim();
  int nElems = mesh().numCells(dim);
  
  Array<int> ssLabels = mesh().getAllLabelsForDimension(dim-1).elements();
  int numSS = 0;

  for (int ss=0; ss<ssLabels.size(); ss++) 
  {
    if (ssLabels[ss] != 0) numSS++;
  }
  
  Array<int> nsLabels = mesh().getAllLabelsForDimension(0).elements();
  int numNS = 0;

  for (int ns=0; ns<nsLabels.size(); ns++) 
  {
    if (nsLabels[ns] != 0) numNS++;
  }

  
  ierr = ex_put_init(
    exoid, 
    filename().c_str(), 
    dim,
    mesh().numCells(0), 
    nElems,
    mesh().numLabels(dim), 
    numNS,
    numSS
    );

  char* qa_record[1][4];
  qa_record[0][0] = "Sundance";
  qa_record[0][1] = "sundance";
  qa_record[0][2] = "date";
  qa_record[0][3] = "time";

  ierr = ex_put_qa(exoid, 1, qa_record);
  TEST_FOR_EXCEPT(ierr < 0);

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

  TEST_FOR_EXCEPT(ierr < 0);

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

  TEST_FOR_EXCEPT(ierr < 0);

  /* write the element blocks */
  Array<int> blockLabels = mesh().getAllLabelsForDimension(dim).elements();
  cout << "block labels =" << blockLabels << endl;
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
    offset(nodeLIDs);
    ierr = ex_put_elem_block(
      exoid, blockLabels[b], eType.c_str(), 
      numElemsThisBlock, nodesPerElem, numBlockAttr
      );

    ierr = ex_put_elem_conn(exoid, blockLabels[b], &(nodeLIDs[0]));
  }


  TEST_FOR_EXCEPT(ierr < 0);
  
  /* write the side sets */
  cout << "writing side sets" << endl;
  
  for (unsigned int ss=0; ss<ssLabels.size(); ss++)
  {
    if (ssLabels[ss]==0) continue;
    cout << "ss=" << ss << " label=" << ssLabels[ss] << endl;
    Array<int> sideLIDs;
    Array<int> elemLIDs;
    Array<int> facets;

    mesh().getLIDsForLabel(dim-1, ssLabels[ss], sideLIDs);
    mesh().getMaxCofacetLIDs(dim-1, sideLIDs, elemLIDs, facets);

    int numSides = sideLIDs.size();
    int numDists = 0;

    offset(sideLIDs);
    offset(elemLIDs);
    offset(facets);

    ierr = ex_put_side_set_param(exoid, ssLabels[ss], numSides, numDists);
    ierr = ex_put_side_set(exoid, ssLabels[ss], &(elemLIDs[0]), &(facets[0]));
  }
  
  TEST_FOR_EXCEPT(ierr < 0);
  cout << "writing node sets" << endl;

  
  /* write the node sets */

  for (unsigned int ns=0; ns<nsLabels.size(); ns++)
  {
    if (nsLabels[ns]==0) continue;
    cout << "ns=" << ns << " label=" << nsLabels[ns] << endl;
    Array<int> nodeLIDs;

    mesh().getLIDsForLabel(0, nsLabels[ns], nodeLIDs);

    int numNodes = nodeLIDs.size();
    int numDists = 0;

    offset(nodeLIDs);

    ierr = ex_put_node_set_param(exoid, nsLabels[ns]+1, numNodes, numDists);
    
    ierr = ex_put_node_set(exoid, nsLabels[ns]+1, &(nodeLIDs[0]));
  }

  TEST_FOR_EXCEPT(ierr < 0);
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



void ExodusWriter::writeParallelInfo(const string& parfile) const 
{
  std::ofstream os(parfile.c_str());

  int dim = mesh().spatialDim();
  int nCells = mesh().numCells(dim);
  int nPts = mesh().numCells(0);

  os << myRank() << " " << nProc() << std::endl;

  os << nPts << std::endl;
  for (int i=0; i<nPts; i++)
    {
      os << i << " " << mesh().mapLIDToGID(0,i) 
         << " " << mesh().ownerProcID(0,i) << std::endl;
    }

  os << nCells << std::endl;
  for (int c=0; c<nCells; c++)
    {
      os << c << " " << mesh().mapLIDToGID(dim,c) 
         << " " << mesh().ownerProcID(dim,c) << std::endl;
    }

  for (int i=0; i<comments().length(); i++)
    {
      os << "# " << comments()[i] << std::endl;
    }
}
