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

#include "SundanceSpatiallyHomogeneousDOFMapBase.hpp"
#include "SundanceTabs.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;


SpatiallyHomogeneousDOFMapBase::SpatiallyHomogeneousDOFMapBase(const Mesh& mesh)
  : DOFMapBase(mesh)
{}

void SpatiallyHomogeneousDOFMapBase::getHomogeneousCellBatches(const RefCountPtr<Array<int> >& inputCells,
                                                           Array<RefCountPtr<Array<int> > >& cellBatches,
                                                           Array<int>& homogSubregionIndices) const
{
  homogSubregionIndices.resize(1);
  cellBatches.resize(1);
  homogSubregionIndices[0] = 0;
  cellBatches[0] = inputCells;
}


void SpatiallyHomogeneousDOFMapBase::print(ostream& os) const
{
  int myRank = mesh().comm().getRank();

  Tabs tabs;
  int dim = mesh().spatialDim();

  for (int p=0; p<mesh().comm().getNProc(); p++)
    {
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      if (p == myRank)
        {
          os << tabs << 
            "========= DOFMap on proc p=" << p << " =============" << endl;
          for (int d=dim; d>=0; d--)
            {
              Tabs tabs1;
              os << tabs1 << "dimension = " << d << endl;
              for (int c=0; c<mesh().numCells(d); c++)
                {
                  Tabs tabs2;
                  os << tabs2 << "Cell d=" << d << " LID=" << c << " GID=" 
                     << mesh().mapLIDToGID(d, c);
                  if (d==0) 
                    {
                      os << " x=" << mesh().nodePosition(c) << endl;
                    }
                  else 
                    {
                      Array<int> facetLIDs;
                      Array<int> facetDirs;
                      mesh().getFacetArray(d, c, 0, facetLIDs, facetDirs);
                      Array<int> facetGIDs(facetLIDs.size());
                      for (unsigned int v=0; v<facetLIDs.size(); v++)
                        {
                          facetGIDs[v] = mesh().mapLIDToGID(0, facetLIDs[v]);
                        }
                      os << " nodes LIDs=" << facetLIDs << " GIDs=" << facetGIDs
                         << endl;
                    }
                  for (int b=0; b<nBasisChunks(0); b++)
                    {
                      for (unsigned int f=0; f<funcID(0,b).size(); f++)
                        {
                          Tabs tabs3;
                          Array<int> dofs;
                          getDOFsForCell(d, c, funcID(0,b)[f], dofs);
                          os << tabs3 << "f=" << funcID(0,b)[f] << " " 
                             << dofs << endl;
                          if (false)
                            {
                              os << tabs3 << "{";
                              for (unsigned int i=0; i<dofs.size(); i++)
                                {
                                  if (i != 0) os << ", ";
                                  if (isLocalDOF(dofs[i])) os << "L";
                                  else os << "R";
                                }
                              os << "}" << endl;
                            }
                        }
                    }
                }
            }
        }
      mesh().comm().synchronize();
      mesh().comm().synchronize();
    }
}
