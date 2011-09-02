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

#include "SundancePeriodicSingleCellMesh1D.hpp"

#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "PlayaMPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;
using Playa::MPIComm;
using Playa::MPIContainerComm;


PeriodicSingleCellMesh1D::PeriodicSingleCellMesh1D(double xMin, double xMax)
  : MeshBase(1, MPIComm::self(), ExodusMeshOrder),
    xMin_(xMin),
    xMax_(xMax),
    vert_(0),
    labels_(2)
{
  labels_[0] = 0;
  labels_[1] = 0;
}

int PeriodicSingleCellMesh1D::numCells(int cellDim) const
{
  switch(cellDim)
  {
    case 0 :
      return 1;
    case 1:
      return 1;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return -1; // -Wall
}


Point PeriodicSingleCellMesh1D::nodePosition(int i) const 
{
  TEST_FOR_EXCEPT(i != 0);
  return xMin_;
}


const double* PeriodicSingleCellMesh1D::nodePositionView(int i) const 
{
  TEST_FOR_EXCEPT(i != 0);
  return &(xMin_);
}

void PeriodicSingleCellMesh1D::getJacobians(int cellDim, const Array<int>& cellLID,
    CellJacobianBatch& jBatch) const
{
  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), std::logic_error,
    "cellDim=" << cellDim 
    << " is not in expected range [0, " << spatialDim()
    << "]");

  jBatch.resize(cellLID.size(), spatialDim(), cellDim);

  int nCells = cellLID.size();

  TEST_FOR_EXCEPT(nCells != 1);

  if (cellDim==0)
  {
    for (int i=0; i<nCells; i++)
    {
      double* detJ = jBatch.detJ(i);
      *detJ = 1.0;
    }
  }
  else
  {
    for (int i=0; i<nCells; i++)
    {
      double* J = jBatch.jVals(i);
      J[0] = fabs(xMin_-xMax_);
    }
  }
}


void PeriodicSingleCellMesh1D::getCellDiameters(int cellDim, const Array<int>& cellLID,
  Array<double>& cellDiameters) const
{
  cellDiameters.resize(1);
  
  TEST_FOR_EXCEPT(cellDim != 1);

  cellDiameters[0] = ::fabs(xMax_-xMin_);
}

void PeriodicSingleCellMesh1D::pushForward(int cellDim, const Array<int>& cellLID,
    const Array<Point>& refQuadPts,
    Array<Point>& physQuadPts) const
{
  TEST_FOR_EXCEPT(cellDim < 0 || cellDim > 1);

  TEST_FOR_EXCEPT(cellLID.size() > 1);

  if (cellDim==1)
  {
    if (physQuadPts.size() > 0) physQuadPts.resize(0);
    physQuadPts.reserve(refQuadPts.size() * cellLID.size());
    
    for (int i=0; i<cellLID.size(); i++)
    {
      double h = xMax_ - xMin_;
      for (int q=0; q<refQuadPts.size(); q++)
      {
        physQuadPts.append(xMin_ + refQuadPts[q][0] * h);
      }
    }
  }
  else
  {
    for (int i=0; i<cellLID.size(); i++)
    {
      physQuadPts.append(xMin_);
    }
  }
}

int PeriodicSingleCellMesh1D::numFacets(int cellDim, int cellLID,
    int facetDim) const
{
  TEST_FOR_EXCEPT(cellLID != 0);
  if (cellDim == 1 && facetDim==0) return 1;
  return 0;
}


    
int PeriodicSingleCellMesh1D::facetLID(int cellDim, int cellLID,
  int facetDim, int facetIndex,
  int& facetOrientation) const
{
  TEST_FOR_EXCEPT(cellLID < 0 || cellLID >= 1);

  TEST_FOR_EXCEPT(cellDim != 1);
  TEST_FOR_EXCEPT(facetDim != 0);
  TEST_FOR_EXCEPT(facetIndex < 0);
  TEST_FOR_EXCEPT(facetIndex > 1);

  return vert_;
}

void PeriodicSingleCellMesh1D::getFacetLIDs(int cellDim,
    const Array<int>& cellLID,
    int facetDim,
    Array<int>& facetLID,
    Array<int>& facetSign) const
{
  facetLID.resize(2*cellLID.size());
  facetSign.resize(2*cellLID.size());

  for (int i=0; i<cellLID.size(); i++) 
  {
    facetLID[2*i] = this->facetLID(cellDim, cellLID[i], facetDim, 0, facetSign[2*i]);
    facetLID[2*i+1] = this->facetLID(cellDim, cellLID[i], facetDim, 1, facetSign[2*i]);
  }
}


const int* PeriodicSingleCellMesh1D::elemZeroFacetView(int cellLID) const
{
  TEST_FOR_EXCEPT(cellLID != 0);
  return &vert_;
}


int PeriodicSingleCellMesh1D::numMaxCofacets(int cellDim, int cellLID) const
{
  TEST_FOR_EXCEPT(cellDim != 0);
  return 1;
}

int PeriodicSingleCellMesh1D::maxCofacetLID(int cellDim, int cellLID,
    int cofacetIndex,
    int& facetIndex) const
{
  TEST_FOR_EXCEPT(cellDim != 0 || cellLID != 0);
  facetIndex = 0;
  return 0;
}

void PeriodicSingleCellMesh1D::getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const
{
  TEST_FOR_EXCEPT(true);
}

void PeriodicSingleCellMesh1D::getCofacets(int cellDim, int cellLID,
    int cofacetDim, Array<int>& cofacetLIDs) const
{
  TEST_FOR_EXCEPT(cellDim != 0);
  TEST_FOR_EXCEPT(cofacetDim != 1);

  cofacetLIDs.resize(1);
  cofacetLIDs[0] = 0;
}

int PeriodicSingleCellMesh1D::mapGIDToLID(int cellDim, int globalIndex) const
{
  return globalIndex;
}

bool PeriodicSingleCellMesh1D::hasGID(int cellDim, int globalIndex) const
{
  return globalIndex==0;
}

int PeriodicSingleCellMesh1D::mapLIDToGID(int cellDim, int localIndex) const
{
  return localIndex;
}

CellType PeriodicSingleCellMesh1D::cellType(int cellDim) const
{
  if (cellDim==0) return PointCell;
  else if (cellDim==1) return LineCell;
  else return NullCell;
}

int PeriodicSingleCellMesh1D::label(int cellDim, int cellLID) const
{
  return labels_[cellDim];
}

void PeriodicSingleCellMesh1D::getLabels(int cellDim, const Array<int>& cellLID,
  Array<int>& labels) const
{
  labels.resize(cellLID.size());
  for (int i=0; i<cellLID.size(); i++)
  {
    labels[i] = labels_[cellDim];
  }
}

Set<int> PeriodicSingleCellMesh1D::getAllLabelsForDimension(int cellDim) const
{
  Set<int> rtn;

  rtn.put(labels_[cellDim]);
    
  return rtn;
}

void PeriodicSingleCellMesh1D::setLabel(int cellDim, int cellLID, int label)
{
  labels_[cellDim] = label;
}


void PeriodicSingleCellMesh1D::getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const
{
  cellLIDs.resize(0);
  if (label != labels_[cellDim]) cellLIDs.append(0);
}
