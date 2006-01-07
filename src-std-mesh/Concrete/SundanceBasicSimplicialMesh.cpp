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

#include "SundanceBasicSimplicialMesh.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "TSFObjectWithVerbosity.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;

static Time& batchedFacetGrabTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("batched facet grabbing"); 
  return *rtn;
}

static Time& getJacobianTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("cell Jacobian grabbing"); 
  return *rtn;
}


//#define SKIP_FACES

BasicSimplicialMesh::BasicSimplicialMesh(int dim, const MPIComm& comm)
	: IncrementallyCreatableMesh(dim, comm),
    numCells_(dim+1),
    points_(),
    edgeVerts_(2),
    faceVertLIDs_(dim),
    faceVertGIDs_(dim),
    faceEdges_(dim),
    faceEdgeSigns_(dim),
    elemVerts_(dim+1),
    elemEdges_(),
    elemEdgeSigns_(),
    elemFaces_(dim+1),
    elemFaceRotations_(dim+1),
    vertexSetToFaceIndexMap_(),
    edgeCofacets_(),
    faceCofacets_(),
    vertEdges_(),
    vertCofacets_(),
    vertEdgePartners_(),
    LIDToGIDMap_(dim+1),
    GIDToLIDMap_(dim+1),
    labels_(dim+1),
    ownerProcID_(dim+1),
    faceVertGIDBase_(1),
    hasEdgeGIDs_(false),
    hasFaceGIDs_(false)
{
  verbosity() = MeshSource::classVerbosity();
  estimateNumVertices(1000);
  estimateNumElements(1000);

  /* Set up the pointer giving a view of the face vertex array.
     /* Resize to 1 so that phony dereference will work, then resize to zero to make the
     * new array logically empty */
  faceVertGIDs_.resize(1);
  faceVertGIDBase_[0] = &(faceVertGIDs_.value(0,0));
  faceVertGIDs_.resize(0);

  /* size the element edge lists as appropriate to the mesh's dimension */
  if (spatialDim()==2) 
    {
      elemEdges_.setTupleSize(3);
      elemEdgeSigns_.setTupleSize(3);
    }
  if (spatialDim()==3) 
    {
      elemEdges_.setTupleSize(6);
      elemEdgeSigns_.setTupleSize(6);
    }
}


void BasicSimplicialMesh::getJacobians(int cellDim, const Array<int>& cellLID,
                                       CellJacobianBatch& jBatch) const
{
  TimeMonitor timer(getJacobianTimer());

  int flops = 0 ;
  int nCells = cellLID.size();

  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
                     "cellDim=" << cellDim 
                     << " is not in expected range [0, " << spatialDim()
                     << "]");

  jBatch.resize(cellLID.size(), spatialDim(), cellDim);

  if (cellDim < spatialDim())
    {
      switch(cellDim)
        {
        case 1:
          flops += 3*nCells;
          break;
        case 2:
          // 4 flops for two pt subtractions, 10 for a cross product
          flops += (4 + 10)*nCells;
          break;
        default:
          break;
        }

      for (int i=0; i<nCells; i++)
        {
          int lid = cellLID[i];
          double* detJ = jBatch.detJ(i);
          switch(cellDim)
            {
            case 0:
              *detJ = 1.0;
              break;
            case 1:
              {
                int a = edgeVerts_.value(lid, 0);
                int b = edgeVerts_.value(lid, 1);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                Point dx = pb-pa;
                *detJ = sqrt(dx*dx);
              }
              break;
            case 2:
              {
                int a = faceVertLIDs_.value(lid, 0);
                int b = faceVertLIDs_.value(lid, 1);
                int c = faceVertLIDs_.value(lid, 2);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                const Point& pc = points_[c];
                Point directedArea = cross(pc-pa, pb-pa);
                *detJ = sqrt(directedArea*directedArea);
              }
              break;
            default:
              TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
                                 "cellDim=" << cellDim 
                                 << " in BasicSimplicialMesh::getJacobians()");
            }
        }
    }
  else
    {
      Array<double> J(cellDim*cellDim);
  
      flops += cellDim*cellDim*nCells;

      for (unsigned int i=0; i<cellLID.size(); i++)
        {
          int lid = cellLID[i];
          double* J = jBatch.jVals(i);
          switch(cellDim)
            {
            case 0:
              J[0] = 1.0;
              break;
            case 1:
              {
                int a = elemVerts_.value(lid, 0);
                int b = elemVerts_.value(lid, 1);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                J[0] = fabs(pa[0]-pb[0]);
              }
              break;
            case 2:
              {
                int a = elemVerts_.value(lid, 0);
                int b = elemVerts_.value(lid, 1);
                int c = elemVerts_.value(lid, 2);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                const Point& pc = points_[c];
                J[0] = pb[0] - pa[0];
                J[1] = pc[0] - pa[0];
                J[2] = pb[1] - pa[1];
                J[3] = pc[1] - pa[1];
            
              }
              break;
            case 3:
              {
                int a = elemVerts_.value(lid, 0);
                int b = elemVerts_.value(lid, 1);
                int c = elemVerts_.value(lid, 2);
                int d = elemVerts_.value(lid, 3);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                const Point& pc = points_[c];
                const Point& pd = points_[d];
                J[0] = pb[0] - pa[0];
                J[1] = pc[0] - pa[0];
                J[2] = pd[0] - pa[0];
                J[3] = pb[1] - pa[1];
                J[4] = pc[1] - pa[1];
                J[5] = pd[1] - pa[1];
                J[6] = pb[2] - pa[2];
                J[7] = pc[2] - pa[2];
                J[8] = pd[2] - pa[2];
              }
              break;
            default:
              TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
                                 "cellDim=" << cellDim 
                                 << " in BasicSimplicialMesh::getJacobians()");
            }
        }
    }
  CellJacobianBatch::addFlops(flops);
}



void BasicSimplicialMesh::getCellDiameters(int cellDim, const Array<int>& cellLID,
                                           Array<double>& cellDiameters) const
{
  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
                     "cellDim=" << cellDim 
                     << " is not in expected range [0, " << spatialDim()
                     << "]");

  cellDiameters.resize(cellLID.size());

  if (cellDim < spatialDim())
    {
      for (unsigned int i=0; i<cellLID.size(); i++)
        {
          int lid = cellLID[i];
          switch(cellDim)
            {
            case 0:
              cellDiameters[i] = 1.0;
              break;
            case 1:
              {
                int a = edgeVerts_.value(lid, 0);
                int b = edgeVerts_.value(lid, 1);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                Point dx = pb-pa;
                cellDiameters[i] = sqrt(dx*dx);
              }
              break;
            case 2:
              {
                int a = faceVertLIDs_.value(lid, 0);
                int b = faceVertLIDs_.value(lid, 1);
                int c = faceVertLIDs_.value(lid, 2);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                const Point& pc = points_[c];
                Point directedArea = cross(pc-pa, pb-pa);
                cellDiameters[i] = sqrt(directedArea*directedArea);
              }
              break;
            default:
              TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
                                 "cellDim=" << cellDim 
                                 << " in BasicSimplicialMesh::getCellDiameters()");
            }
        }
    }
  else
    {
      for (unsigned int i=0; i<cellLID.size(); i++)
        {
          int lid = cellLID[i];
          switch(cellDim)
            {
            case 0:
              cellDiameters[i] = 1.0;
              break;
            case 1:
              {
                int a = elemVerts_.value(lid, 0);
                int b = elemVerts_.value(lid, 1);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                cellDiameters[i] = fabs(pa[0]-pb[0]);
              }
              break;
            case 2:
              {
                int a = elemVerts_.value(lid, 0);
                int b = elemVerts_.value(lid, 1);
                int c = elemVerts_.value(lid, 2);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                const Point& pc = points_[c];
                cellDiameters[i] 
                  = (pa.distance(pb) + pb.distance(pc) + pa.distance(pc))/3.0;
              }
              break;
            case 3:
              {
                int a = elemVerts_.value(lid, 0);
                int b = elemVerts_.value(lid, 1);
                int c = elemVerts_.value(lid, 2);
                int d = elemVerts_.value(lid, 3);
                const Point& pa = points_[a];
                const Point& pb = points_[b];
                const Point& pc = points_[c];
                const Point& pd = points_[d];
                
                cellDiameters[i] 
                  = (pa.distance(pb) + pa.distance(pc) + pa.distance(pd)
                     + pb.distance(pc) + pb.distance(pd)
                     + pc.distance(pd))/6.0;
              }
              break;
            default:
              TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
                                 "cellDim=" << cellDim 
                                 << " in BasicSimplicialMesh::getCellDiameters()");
            }
        }
    }
}


void BasicSimplicialMesh::pushForward(int cellDim, const Array<int>& cellLID,
                                      const Array<Point>& refQuadPts,
                                      Array<Point>& physQuadPts) const
{
  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
                     "cellDim=" << cellDim 
                     << " is not in expected range [0, " << spatialDim()
                     << "]");
  int flops = 0;
  int nQuad = refQuadPts.size();
  int nCells = cellLID.size();
  Array<double> J(cellDim*cellDim);

  if (physQuadPts.size() > 0) physQuadPts.resize(0);
  physQuadPts.reserve(cellLID.size() * refQuadPts.size());

  switch(cellDim)
    {
    case 1:
      flops += nCells * (1 + 2*nQuad);
      break;
    case 2:
      if (spatialDim()==2)
        {
          flops += nCells*(4 + 8*nQuad);
        }
      else
        {
          flops += 18*nCells*nQuad;
        }
      break;
    case 3:
      flops += 27*nCells*nQuad;
      break;
    default:
      break;
    }


  for (unsigned int i=0; i<cellLID.size(); i++)
    {
      int lid = cellLID[i];
      switch(cellDim)
        {
        case 0:
          physQuadPts.append(points_[lid]);
          break;
        case 1:
          {
            int a, b;
            if (spatialDim()==1)
              {
                a = elemVerts_.value(lid, 0);
                b = elemVerts_.value(lid, 1);
              }
            else
              {
                a = edgeVerts_.value(lid, 0);
                b = edgeVerts_.value(lid, 1);
              }
            const Point& pa = points_[a];
            const Point& pb = points_[b];
            Point dx = pb-pa;
            for (int q=0; q<nQuad; q++)
              {
                physQuadPts.append(pa + refQuadPts[q][0]*dx);
              }
          }
          break;
        case 2:
          {
            int a,b,c;
            if (spatialDim()==2)
              {
                a = elemVerts_.value(lid, 0);
                b = elemVerts_.value(lid, 1);
                c = elemVerts_.value(lid, 2);
              }
            else
              {
                a = faceVertLIDs_.value(lid, 0);
                b = faceVertLIDs_.value(lid, 1);
                c = faceVertLIDs_.value(lid, 2);
              }
            const Point& pa = points_[a];
            const Point& pb = points_[b];
            const Point& pc = points_[c];

            if (spatialDim()==2)
              {
                J[0] = pb[0] - pa[0];
                J[1] = pc[0] - pa[0];
                J[2] = pb[1] - pa[1];
                J[3] = pc[1] - pa[1];
                for (int q=0; q<nQuad; q++)
                  {
                    physQuadPts.append(pa 
                                       + Point(J[0]*refQuadPts[q][0] +J[1]*refQuadPts[q][1],
                                               J[2]*refQuadPts[q][0] +J[3]*refQuadPts[q][1]) );
                  }
              }
            else
              {
                for (int q=0; q<nQuad; q++)
                  {
                    physQuadPts.append(pa 
                                       + (pb-pa)*refQuadPts[q][0] 
                                       + (pc-pa)*refQuadPts[q][1]);
                  }
              }
            
          }
          break;
        case 3:
          {
            int a = elemVerts_.value(lid, 0);
            int b = elemVerts_.value(lid, 1);
            int c = elemVerts_.value(lid, 2);
            int d = elemVerts_.value(lid, 3);
            const Point& pa = points_[a];
            const Point& pb = points_[b];
            const Point& pc = points_[c];
            const Point& pd = points_[d];
            J[0] = pb[0] - pa[0];
            J[1] = pc[0] - pa[0];
            J[2] = pd[0] - pa[0];
            J[3] = pb[1] - pa[1];
            J[4] = pc[1] - pa[1];
            J[5] = pd[1] - pa[1];
            J[6] = pb[2] - pa[2];
            J[7] = pc[2] - pa[2];
            J[8] = pd[2] - pa[2];
            
            for (int q=0; q<nQuad; q++)
              {
                physQuadPts.append(pa 
                                   + Point(J[0]*refQuadPts[q][0] 
                                           + J[1]*refQuadPts[q][1]
                                           + J[2]*refQuadPts[q][2],
                                           J[3]*refQuadPts[q][0] 
                                           + J[4]*refQuadPts[q][1]
                                           + J[5]*refQuadPts[q][2],
                                           J[6]*refQuadPts[q][0] 
                                           + J[7]*refQuadPts[q][1]
                                           + J[8]*refQuadPts[q][2]));

              }
            
          }
          break;
        default:
          TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
                             "in BasicSimplicialMesh::getJacobians()");
        }
    }

  CellJacobianBatch::addFlops(flops);
}

void BasicSimplicialMesh::estimateNumVertices(int nPts)
{
  points_.reserve(nPts);
  vertCofacets_.reserve(nPts);
  vertEdges_.reserve(nPts);
  vertEdgePartners_.reserve(nPts);

  ownerProcID_[0].reserve(nPts);
  LIDToGIDMap_[0].reserve(nPts);
  GIDToLIDMap_[0] = Hashtable<int,int>(nPts, 0.6);
  labels_[0].reserve(nPts);
}

void BasicSimplicialMesh::estimateNumElements(int nElems)
{
  int nEdges = 0;
  int nFaces = 0;

  if (spatialDim()==3) 
    {
      nFaces = 5*nElems;
      nEdges = 5*nElems;
      labels_[2].reserve(nFaces);
    }
  else if (spatialDim()==2)
    {
      nEdges = 3*nElems;
    }
  
  labels_[1].reserve(nEdges);
  vertexSetToFaceIndexMap_ = Hashtable<VertexView, int>(nFaces);

  edgeVerts_.reserve(nEdges);
  faceVertLIDs_.reserve(nFaces);
  faceVertGIDs_.reserve(nFaces);
  faceEdges_.reserve(nFaces);
  elemVerts_.reserve(nElems);
  elemEdges_.reserve(nElems);
  elemEdgeSigns_.reserve(nElems);
  elemFaces_.reserve(nElems);
  edgeCofacets_.reserve(nEdges);
  faceCofacets_.reserve(nFaces);

  ownerProcID_[spatialDim()].reserve(nElems);
  LIDToGIDMap_[spatialDim()].reserve(nElems);
  GIDToLIDMap_[spatialDim()] = Hashtable<int,int>(nElems, 0.6);
  labels_[spatialDim()].reserve(nElems);

  /* resize to 1 so that phony dereference will work, then resize to zero to make the
   * new array logically empty */
  faceVertGIDs_.resize(1);
  faceVertGIDBase_[0] = &(faceVertGIDs_.value(0,0));
  faceVertGIDs_.resize(0);
}



int BasicSimplicialMesh::numCells(int cellDim) const
{
  return numCells_[cellDim];
}

int BasicSimplicialMesh::ownerProcID(int cellDim, int cellLID) const
{
  return ownerProcID_[cellDim][cellLID];
}

int BasicSimplicialMesh::numFacets(int cellDim, int cellLID, 
                                   int facetDim) const
{
  if (cellDim==1)
    {
      return 2;
    }
  else if (cellDim==2)
    {
      return 3;
    }
  else 
    {
      if (facetDim==0) return 4;
      if (facetDim==1) return 6;
      return 4;
    }
}

void BasicSimplicialMesh::getFacetLIDs(int cellDim, 
                                       const Array<int>& cellLID,
                                       int facetDim,
                                       Array<int>& facetLID,
                                       Array<int>& facetSign) const 
{
  TimeMonitor timer(batchedFacetGrabTimer());

  int nf = numFacets(cellDim, cellLID[0], facetDim);
  facetLID.resize(cellLID.size() * nf);
  if (facetDim > 0) facetSign.resize(cellLID.size() * nf);

  
  if (facetDim==0)
    {
      if (cellDim == spatialDim())
        {
          int ptr=0;
          for (unsigned int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(elemVerts_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                }
            }
        }
      else if (cellDim==1) 
        {
          int ptr=0;
          for (unsigned int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(edgeVerts_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                }
            }
        }
      else if (cellDim==2) 
        {
          int ptr=0;
          for (unsigned int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(faceVertLIDs_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                }
            }
        }
    }
  else if (facetDim==1)
    {
      int ptr=0;
      if (cellDim == spatialDim())
        {
          for (unsigned int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(elemEdges_.value(cellLID[c], 0));
              const int* fsPtr = &(elemEdgeSigns_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                  facetSign[ptr] = fsPtr[f];
                }
            }
        }
      else
        {
          for (unsigned int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(faceEdges_.value(cellLID[c], 0));
              //    const int* fsPtr = &(faceEdgeSigns_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                  //  facetSign[ptr] = fsPtr[f];
                }
            }
        }
    }
  else
    {
      int ptr=0;
      for (unsigned int c=0; c<cellLID.size(); c++)
        {
          const int* fPtr = &(elemFaces_.value(cellLID[c], 0));
          const int* fsPtr = &(elemFaceRotations_.value(cellLID[c], 0));
          for (int f=0; f<nf; f++, ptr++)
            {
              facetLID[ptr] = fPtr[f];
              facetSign[ptr] = fsPtr[f];
            }
        }
    }
}

int BasicSimplicialMesh::facetLID(int cellDim, int cellLID,
                                  int facetDim, int facetIndex, int& facetSign) const 
{

  if (facetDim==0)
    {
      if (cellDim == spatialDim())
        {
          return elemVerts_.value(cellLID, facetIndex);
        }
      else if (cellDim==1) return edgeVerts_.value(cellLID, facetIndex);
      else if (cellDim==2) return faceVertLIDs_.value(cellLID, facetIndex);
    }
  if (facetDim==1)
    {
      if (cellDim==spatialDim())
        {
          facetSign = elemEdgeSigns_.value(cellLID, facetIndex);
          return elemEdges_.value(cellLID, facetIndex);
        }
      else
        {
          //facetSign = faceEdgeSigns_.value(cellLID, facetIndex);
          return faceEdges_.value(cellLID, facetIndex);
        }
    }
  else
    {
      facetSign = elemFaceRotations_.value(cellLID, facetIndex);
      return elemFaces_.value(cellLID, facetIndex);
    }
}


int BasicSimplicialMesh::numCofacets(int cellDim, int cellLID) const 
{

  if (cellDim==0)
    {
      return vertCofacets_[cellLID].length();
    }
  if (cellDim==1)
    {
      return edgeCofacets_[cellLID].length();
    }
  if (cellDim==2)
    {
      return faceCofacets_[cellLID].length();
    }
  return -1; // -Wall
}

int BasicSimplicialMesh::cofacetLID(int cellDim, int cellLID,
                                    int cofacetIndex) const
{

  if (cellDim==0)
    {
      return vertCofacets_[cellLID][cofacetIndex];
    }
  if (cellDim==1)
    {
      return edgeCofacets_[cellLID][cofacetIndex];
    }
  if (cellDim==2)
    {
      return faceCofacets_[cellLID][cofacetIndex];
    }
  return -1; // -Wall
}

int BasicSimplicialMesh::mapGIDToLID(int cellDim, int globalIndex) const
{
  return GIDToLIDMap_[cellDim].get(globalIndex);
}

int BasicSimplicialMesh::mapLIDToGID(int cellDim, int localIndex) const
{
  return LIDToGIDMap_[cellDim][localIndex];
}

CellType BasicSimplicialMesh::cellType(int cellDim) const
{

  switch(cellDim)
    {
    case 0:
      return PointCell;
    case 1:
      return LineCell;
    case 2:
      return TriangleCell;
    case 3:
      return TetCell;
    default:
      return NullCell; // -Wall
    }
}

int BasicSimplicialMesh::label(int cellDim, 
                               int cellLID) const
{
  return labels_[cellDim][cellLID];
}

int BasicSimplicialMesh::addVertex(int globalIndex, const Point& x,
                                   int ownerProcID, int label)
{

  SUNDANCE_OUT(this->verbosity() > VerbLow,
               "BSM points are " << points_.toString());

  int lid = points_.length();
  points_.append(x);

  SUNDANCE_OUT(this->verbosity() > VerbLow,
               "BSM added point " << x << " lid = " << lid);

  numCells_[0]++;

  LIDToGIDMap_[0].append(globalIndex);
  GIDToLIDMap_[0].put(globalIndex, lid);

  ownerProcID_[0].append(ownerProcID);
  labels_[0].append(label);

  vertCofacets_.resize(points_.length());
  vertEdges_.resize(points_.length());
  vertEdgePartners_.resize(points_.length());
  
  return lid;
}

int BasicSimplicialMesh::addElement(int globalIndex, 
                                    const Array<int>& vertGID,
                                    int ownerProcID, int label)
{
  SUNDANCE_VERB_HIGH("adding element " << globalIndex 
                     << " with vertices:" << vertGID);

  /* 
   * do basic administrative steps for the new element: 
   * set LID, label, and procID; update element count. 
   */
  int lid = elemVerts_.length();
  elemEdgeSigns_.resize(lid+1);

  numCells_[spatialDim()]++;

  LIDToGIDMap_[spatialDim()].append(globalIndex);
  GIDToLIDMap_[spatialDim()].put(globalIndex, lid);
  labels_[spatialDim()].append(label);
  ownerProcID_[spatialDim()].append(ownerProcID);


  /* these little arrays will get used repeatedly, so make them static
   * to save a few cycles. */
  static Array<int> edges;
  static Array<int> faces;
  static Array<int> faceRotations;
  static Array<int> vertLID;

  /* find the vertex LIDs given the input GIDs */
  vertLID.resize(vertGID.size());
  for (unsigned int i=0; i<vertGID.size(); i++) 
    {
      vertLID[i] = GIDToLIDMap_[0].get(vertGID[i]);
    }
  
  /* 
   * Now comes the fun part: creating edges and faces for the 
   * new element, and registering it as a cofacet of its
   * lower-dimensional facets. 
   */
  

  if (spatialDim()==1)  
    {
      /* In 1D, there are neither edges nor faces. */
      edges.resize(0);
      faces.resize(0);
      /* register the new element as a cofacet of its vertices. */
      vertCofacets_[vertLID[0]].append(lid);
      vertCofacets_[vertLID[1]].append(lid);
    }
  if (spatialDim()==2)
    {
      int edgeSign;

      /* In 2D, we need to define edges but not faces for the new element. */
      edges.resize(3);
      faces.resize(0);

      /* add the edges and define the relative orientations of the edges */
      edges[0] = addEdge(vertLID[0], vertLID[1], lid, 0);
      edges[1] = addEdge(vertLID[1], vertLID[2], lid, 1);
      edges[2] = addEdge(vertLID[2], vertLID[0], lid, 2);

      /* register the new element as a cofacet of its vertices. */
      vertCofacets_[vertLID[0]].append(lid);
      vertCofacets_[vertLID[1]].append(lid);
      vertCofacets_[vertLID[2]].append(lid);

    }
  else if (spatialDim()==3)
    {
      int edgeSign;

      /* In 3D, we need to define edges and faces for the new element. */
      edges.resize(6);
      faces.resize(4);
      faceRotations.resize(4);

      /* add the edges and define the relative orientations of the edges */
      edges[0] = addEdge(vertLID[0], vertLID[1], lid, 0);
      edges[1] = addEdge(vertLID[1], vertLID[2], lid, 1);
      edges[2] = addEdge(vertLID[2], vertLID[0], lid, 2);
      edges[3] = addEdge(vertLID[0], vertLID[3], lid, 3);
      edges[4] = addEdge(vertLID[1], vertLID[3], lid, 4);
      edges[5] = addEdge(vertLID[2], vertLID[3], lid, 5);

      /* register the new element as a cofacet of its vertices. */
      vertCofacets_[vertLID[0]].append(lid);
      vertCofacets_[vertLID[1]].append(lid);
      vertCofacets_[vertLID[2]].append(lid);
      vertCofacets_[vertLID[3]].append(lid);

      /* add the faces and define the relative orientations of the faces */
      faces[0] = addFace(vertLID[0], vertLID[1], vertLID[3], 
                         edges[0], edges[4], edges[3],
                         faceRotations[0]);


      faces[1] = addFace(vertLID[1], vertLID[2], vertLID[3], 
                         edges[1], edges[5], edges[4],
                         faceRotations[1]);

      faces[2] = addFace(vertLID[0], vertLID[3], vertLID[2], 
                         edges[3], edges[5], edges[2],
                         faceRotations[2]);

      faces[3] = addFace(vertLID[2], vertLID[1], vertLID[0], 
                         edges[1], edges[0], edges[2],
                         faceRotations[3]);
    }

  elemVerts_.append(vertLID);
  if (edges.length() > 0) elemEdges_.append(edges);

  if (faces.length() > 0) 
    {
      elemFaces_.append(faces);
      elemFaceRotations_.append(faceRotations);
    }
  for (int i=0; i<edges.length(); i++) edgeCofacets_[edges[i]].append(lid);


  for (int i=0; i<faces.length(); i++) faceCofacets_[faces[i]].append(lid);

  return lid;
}

int BasicSimplicialMesh::addFace(int v1, int v2, int v3, 
                                 int e1, int e2, int e3,
                                 int& rotation)
{
  /* create static data for workspace arrays and pointers that will be
   * used repeatedly */
  static Array<int> sortedVertGIDs(3);
  static Array<int> reorderedVertLIDs(3);
  static Array<int> reorderedEdgeLIDs(3);
  //  static Array<int> edgeSigns(3);
  static int* sortedGIDs = &(sortedVertGIDs[0]);
  static int* reorderedLIDs = &(reorderedVertLIDs[0]);
  static int* reorderedEdges = &(reorderedEdgeLIDs[0]);
  //static int* edgeSigns = &(edgeSigns[0]);

  /* First we check whether the face already exists, and
   * along the way determine the orientation of the new element's 
   * face relative to the absolute orientation of the face. */
  int lid = checkForExistingFace(v1, v2, v3, e1, e2, e3,
                                 sortedGIDs, reorderedLIDs, reorderedEdges,
                                  rotation);

  if (lid >= 0) /* if the face already exists, we're done */
    {
      /* return the LID of the pre-existing face */
      return lid;
    }
  else /* we need to register the face */
    {
      /* get a LID for the new face */
      lid = faceVertLIDs_.length();
      
      /* record the new face's vertex sets (both GID and LID, ordered by GID)
       * and its edges (reordered to the the sorted-GID orientation) */
      faceVertGIDs_.append(sortedGIDs, 3);
      faceVertLIDs_.append(reorderedLIDs, 3);
      faceEdges_.append(reorderedEdges, 3);
      ownerProcID_[2].append(ownerProcID_[0][reorderedLIDs[0]]);


      /* the face doesn't yet have a label, so set it to zero */
      labels_[2].append(0);
      
      /* update the view pointer to stay in synch with the resized
       * face vertex GID array */
      faceVertGIDBase_[0] = &(faceVertGIDs_.value(0,0));
      
      /* create a vertex set view that points to the new face's 
       * vertex GIDs. The first argument gives a pointer to the base of the
       * vertex GID array. The view is offset from the base by the face's
       * LID, and is of length 3. */
      VertexView face(&(faceVertGIDBase_[0]), lid, 3);
      
      /* record this vertex set in the hashtable of 
       * existing face vertex sets */
      vertexSetToFaceIndexMap_.put(face, lid);
      
      /* create space for the new face's cofacets */
      faceCofacets_.resize(lid+1);
      
      /* update the cell count */
      numCells_[spatialDim()-1]++;

      /* return the LID of the new face */
      return lid;
    }

}


void BasicSimplicialMesh::getSortedFaceVertices(int a, int b, int c,
                                                int la, int lb, int lc,
                                                int ea, int eb, int ec, 
                                                int* sortedVertGIDs,
                                                int* reorderedVertLIDs,
                                                int* reorderedEdgeLIDs,
                                                int& rotation) const
{
  /* 
   * Do a hand-coded sort of a 3-tuple of vertex GIDs, reordering tuples of vertex LIDs
   * and edge LIDs accordingly.
   */
  if (a < b)
    {
      if (c < a) 
        {
          fillSortedArray(c,a,b,sortedVertGIDs); 
          fillSortedArray(lc,la,lb,reorderedVertLIDs); 
          fillSortedArray(ec,ea,eb,reorderedEdgeLIDs); 
          rotation=3;
        }
      else if (c > b) 
        {
          fillSortedArray(a,b,c,sortedVertGIDs); 
          fillSortedArray(la,lb,lc,reorderedVertLIDs); 
          fillSortedArray(ea,eb,ec,reorderedEdgeLIDs); 
          rotation=1;
        }
      else 
        {
          fillSortedArray(a,c,b,sortedVertGIDs); 
          fillSortedArray(la,lc,lb,reorderedVertLIDs); 
          fillSortedArray(ea,ec,eb,reorderedEdgeLIDs); 
          rotation=-1;
        }
    }
  else /* b < a */
    {
      if (c < b) 
        {
          fillSortedArray(c,b,a,sortedVertGIDs); 
          fillSortedArray(lc,lb,la,reorderedVertLIDs); 
          fillSortedArray(ec,eb,ea,reorderedEdgeLIDs); 
          rotation=-3;
        }
      else if (c > a) 
        {
          fillSortedArray(b,a,c,sortedVertGIDs); 
          fillSortedArray(lb,la,lc,reorderedVertLIDs); 
          fillSortedArray(eb,ea,ec,reorderedEdgeLIDs); 
          rotation=-2;
        }
      else
        {
          fillSortedArray(b,c,a,sortedVertGIDs); 
          fillSortedArray(lb,lc,la,reorderedVertLIDs); 
          fillSortedArray(eb,ec,ea,reorderedEdgeLIDs); 
          rotation=2;
        }
    }
}


void BasicSimplicialMesh::getSortedFaceVertices(int a, int b, int c,
                                                int* sortedVertGIDs) const
{
  /* 
   * Do a hand-coded sort of a 3-tuple of vertex GIDs
   */
  if (a < b)
    {
      if (c < a) 
        {
          fillSortedArray(c,a,b,sortedVertGIDs); 
        }
      else if (c > b) 
        {
          fillSortedArray(a,b,c,sortedVertGIDs); 
        }
      else 
        {
          fillSortedArray(a,c,b,sortedVertGIDs); 
        }
    }
  else /* b < a */
    {
      if (c < b) 
        {
          fillSortedArray(c,b,a,sortedVertGIDs); 
        }
      else if (c > a) 
        {
          fillSortedArray(b,a,c,sortedVertGIDs); 
        }
      else
        {
          fillSortedArray(b,c,a,sortedVertGIDs); 
        }
    }
}


int BasicSimplicialMesh::checkForExistingEdge(int vertLID1, int vertLID2)
{
  SUNDANCE_OUT(this->verbosity() > VerbHigh, 
               "p=" << comm().getRank() << " finding edge for verts="
               << vertLID1 << ", " << vertLID2);
  const Array<int>& edgePartners = vertEdgePartners_[vertLID1];
  for (int i=0; i<edgePartners.length(); i++)
    {
      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   "p=" << comm().getRank() << " checking partner="
                   << edgePartners[i]);
      if (edgePartners[i] == vertLID2)
        {
          SUNDANCE_OUT(this->verbosity() > VerbHigh, "p=" << comm().getRank() 
                       << " found match!");
          return vertEdges_[vertLID1][i];
        }
    }
  SUNDANCE_OUT(this->verbosity() > VerbHigh, "p=" << comm().getRank() 
               << " no match!");
  return -1;
}

int BasicSimplicialMesh::checkForExistingFace(int v1, int v2, int v3,
                                              int e1, int e2, int e3,
                                              int* sortedVertGIDs,
                                              int* reorderedVertLIDs,
                                              int* reorderedEdgeLIDs,
                                              int& rotation) 
{
  /* sort the face's vertices by global ID, reordering the vertex LIDs and
   * edge LIDs to follow. */
  int g1 = LIDToGIDMap_[0][v1];
  int g2 = LIDToGIDMap_[0][v2];
  int g3 = LIDToGIDMap_[0][v3];

  getSortedFaceVertices(g1, g2, g3, v1, v2, v3, e1, e2, e3, 
                        sortedVertGIDs, reorderedVertLIDs, reorderedEdgeLIDs,
                        rotation);

  /* see if the face is already in existence by checking for the
   * presence of the same ordered set of vertices. */
  VertexView inputFaceView(&(sortedVertGIDs), 0, 3);

  if (vertexSetToFaceIndexMap_.containsKey(inputFaceView)) 
    {
      /* return the existing face's LID */
      return vertexSetToFaceIndexMap_.get(inputFaceView);
    }
  else 
    {
      /* return -1 as an indication that the face does not yet exist */
      return -1;
    }
}

int BasicSimplicialMesh::lookupFace(int g1, int g2, int g3) 
{
  static Array<int> sortedVertGIDs(3);
  static int* sortedVertGIDPtr = &(sortedVertGIDs[0]);

  /* sort the face's vertices by global ID */
  getSortedFaceVertices(g1, g2, g3, sortedVertGIDPtr);

  /* see if the face is already in existence by checking for the
   * presence of the same ordered set of vertices. */
  VertexView inputFaceView(&(sortedVertGIDPtr), 0, 3);

  if (vertexSetToFaceIndexMap_.containsKey(inputFaceView)) 
    {
      /* return the existing face's LID */
      return vertexSetToFaceIndexMap_.get(inputFaceView);
    }
  else 
    {
      /* return -1 as an indication that the face does not yet exist */
      return -1;
    }
}
                                                
int BasicSimplicialMesh::addEdge(int v1, int v2, 
                                 int elemLID, int myFacetNumber)
{
  /* 
   * First we ask if this edge already exists in the mesh. If so,
   * we're done. 
   */
  int lid = checkForExistingEdge(v1, v2);

  if (lid >= 0) 
    {
      /* determine the sign of the existing edge */
      int g1 = LIDToGIDMap_[0][v1];
      int g2 = LIDToGIDMap_[0][v2];
      int edgeSign;
      if (g2 > g1) 
        {
          edgeSign = 1;
        }
      else 
        {
          edgeSign = -1;
        }
      elemEdgeSigns_.value(elemLID, myFacetNumber) = edgeSign;
      /* return the LID of the pre-existing edge */
      return lid;
    }
  else
    {
      /* get the LID of the new edge */
      lid = edgeVerts_.length();

      /* create the new edge, oriented from lower to higher vertex GID */
      int g1 = LIDToGIDMap_[0][v1];
      int g2 = LIDToGIDMap_[0][v2];
      int edgeSign;
      if (g2 > g1) 
        {
          edgeVerts_.append(tuple(v1,v2));
          ownerProcID_[1].append(ownerProcID_[0][v2]);
          edgeSign = 1;
        }
      else 
        {
          edgeVerts_.append(tuple(v2,v1));
          ownerProcID_[1].append(ownerProcID_[0][v1]);
          edgeSign = -1;
        }
      elemEdgeSigns_.value(elemLID, myFacetNumber) = edgeSign;

      /* register the new edge with its vertices */
      vertEdges_[v1].append(lid);
      vertEdgePartners_[v1].append(v2);
      vertEdges_[v2].append(lid);
      vertEdgePartners_[v2].append(v1);
      /* create storage for the cofacets of the new edge */
      edgeCofacets_.resize(lid+1);
      
      /* the new edge is so far unlabeled, so set its label to zero */
      labels_[1].append(0);
      
      /* update the edge count */
      numCells_[1]++;
    }
  /* return the LID of the edge */
  return lid;
}


void BasicSimplicialMesh::synchronizeGIDNumbering(int dim, int localCount) 
{
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               "sharing offsets for GID numbering for dim=" << dim);
  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               "I have " << localCount << " cells");
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  Array<int> gidOffsets;
  int total;
  MPIContainerComm<int>::accumulate(localCount, gidOffsets, total, comm());
  int myOffset = gidOffsets[comm().getRank()];

  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               "back from MPI accumulate");

  for (int i=0; i<numCells(dim); i++)
    {
      if (LIDToGIDMap_[dim][i] == -1) continue;
      LIDToGIDMap_[dim][i] += myOffset;
      GIDToLIDMap_[dim].put(LIDToGIDMap_[dim][i], i);
    }
  SUNDANCE_OUT(this->verbosity() > VerbMedium, 
               "done sharing offsets for GID numbering for dim=" << dim);
}


void BasicSimplicialMesh::assignIntermediateCellOwners(int cellDim)
{
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  SUNDANCE_OUT(this->verbosity() > VerbMedium,  
               "synchronizing cells of dimension " << cellDim);

  int myRank = comm().getRank();

  /* list of vertex GIDs which serve to identify the cells whose
   * GIDs are needed */
  Array<Array<int> > outgoingRequests(comm().getNProc());
  Array<Array<int> > incomingRequests(comm().getNProc());

  /* lists of cells for which GIDs are needed from each processor */
  Array<Array<int> > outgoingRequestLIDs(comm().getNProc());

  /* lists of GIDs to be communicated. */
  Array<Array<int> > outgoingGIDs(comm().getNProc());
  Array<Array<int> > incomingGIDs(comm().getNProc());

  /* Define these so we're not constantly dereferencing the owners arrays */
  const Array<int>& ptOwners = ownerProcID_[0];
  Array<int>& cellOwners = ownerProcID_[cellDim];

  /* run through the cells, assigning GID numbers to the locally owned ones
   * and compiling a list of the remotely owned ones. */
  int localCount = 0;
  ArrayOfTuples<int>* cellVerts;
  if (cellDim==2) cellVerts = &(faceVertLIDs_);
  else cellVerts = &(edgeVerts_);

  LIDToGIDMap_[cellDim].resize(numCells(cellDim));
  
  SUNDANCE_OUT(this->verbosity() > VerbMedium,  
               "starting loop over cells");

  for (int i=0; i<numCells(cellDim); i++)
    {  
      SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                   "p=" << myRank 
                   <<" cell " << i);

      /* If the cell is locally owned, give it a temporrary GID. This GID
       * is temporary because it's one of a sequence numbered from zero
       * rather than the lowest GID on this processor. We'll offset
       * all of these GIDs once we know how many cells this processor owns. */
      int owner = cellOwners[i];
      if (owner==myRank)
        {
          LIDToGIDMap_[cellDim][i] = localCount++;
        }
      else /* If the cell is remotely owned, we can't give it a GID now. 
            * Give it a GID of -1 to mark it as unresolved, and then
            * add it to a list of cells for which we'll need remote GIDs. 
            * To get a GID, we have to look up this cell on another processor.
            * We can't use the LID to identify the cell, since other
            * procs won't know our LID numbering. Thus we send the
            * GIDs of the cell's vertices, which are sufficient to 
            * identify the cell on the other processor. */
        {
          LIDToGIDMap_[cellDim][i] = -1;
          SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                       "p=" << myRank << " adding LID=" << i 
                       << " to req list");
          outgoingRequestLIDs[owner].append(i);
          for (int v=0; v<=cellDim; v++)
            {
              int ptLID = cellVerts->value(i, v);
              int ptGID = LIDToGIDMap_[0][ptLID];
              SUNDANCE_OUT(this->verbosity() > VerbHigh, 
                           "p=" << myRank << " adding pt GID=" << ptGID
                           << " to req list");
              outgoingRequests[owner].append(ptGID);
            }
        }
    }

  /* Now that we know how many cells are locally owned, we can 
   * set up a global GID numbering */
  synchronizeGIDNumbering(cellDim, localCount);

  
  /* We now share requests for cells for which remote GIDs are needed */
  SUNDANCE_OUT(this->verbosity() > VerbHigh,
               "p=" << myRank << "sending requests: " << outgoingRequests);
  MPIContainerComm<int>::allToAll(outgoingRequests,
                                  incomingRequests,
                                  comm());
  SUNDANCE_OUT(this->verbosity() > VerbHigh,
               "p=" << myRank << "recv'd requests: " << incomingRequests);

  /* Answer the requests incoming to this processor */
  for (int p=0; p<comm().getNProc(); p++)
    {
      const Array<int>& requestsFromProc = incomingRequests[p];
      
      for (unsigned int c=0; c<requestsFromProc.size(); c+=(cellDim+1))
        {
          SUNDANCE_OUT(this->verbosity() > VerbHigh,
                       "p=" << myRank << "processing request c=" 
                       << c/(cellDim+1));
          /* The request for each cell is a tuple of vertex GIDs. 
           * Find the LID of the local cell that contains them */
          int cellLID;
          if (cellDim==1) 
            {
              int vertLID1 = mapGIDToLID(0, requestsFromProc[c]);
              int vertLID2 = mapGIDToLID(0, requestsFromProc[c+1]);
              SUNDANCE_VERB_HIGH("p=" << myRank << " vertex GIDs are "
                                 <<  requestsFromProc[c+1]
                                 << ", " << requestsFromProc[c]);
              cellLID = checkForExistingEdge(vertLID1, vertLID2);
            }
          else
            {
              
              SUNDANCE_VERB_HIGH("p=" << myRank << " vertex GIDs are "
                                 <<  requestsFromProc[c]
                                 << ", " << requestsFromProc[c+1]
                                 << ", " << requestsFromProc[c+2]);
              cellLID = lookupFace(requestsFromProc[c],
                                   requestsFromProc[c+1],
                                   requestsFromProc[c+2]);
            }
          
          SUNDANCE_VERB_HIGH("p=" << myRank << "cell LID is " << cellLID);
          /* Finally, we get the cell's GID and append to the list 
           * of GIDs to send back as answers. */
          outgoingGIDs[p].append(mapLIDToGID(cellDim, cellLID));
        }
    }


  /* We now share the GIDs for the requested cells */
  SUNDANCE_OUT(this->verbosity() > VerbHigh,
               "p=" << myRank << "sending GIDs: " << outgoingGIDs);
  MPIContainerComm<int>::allToAll(outgoingGIDs,
                                  incomingGIDs,
                                  comm());
  SUNDANCE_OUT(this->verbosity() > VerbHigh,
               "p=" << myRank << "recv'ing GIDs: " << incomingGIDs);

  /* Now assign the cell GIDs we've recv'd from from the other procs */
  for (int p=0; p<comm().getNProc(); p++)
    {
      const Array<int>& gidsFromProc = incomingGIDs[p];
      for (unsigned int c=0; c<gidsFromProc.size(); c++)
        {
          int cellLID = outgoingRequestLIDs[p][c];
          int cellGID = incomingGIDs[p][c];
          SUNDANCE_OUT(this->verbosity() > VerbHigh,
                       "p=" << myRank << 
                       " assigning GID=" << cellGID << " to LID=" 
                       << cellLID);
          LIDToGIDMap_[cellDim][cellLID] = cellGID;
          GIDToLIDMap_[cellDim].put(cellGID, cellLID);
        }
    }

  SUNDANCE_OUT(this->verbosity() > VerbMedium,  
               "p=" << myRank << "done synchronizing cells of dimension " << cellDim);

  if (cellDim==1) hasEdgeGIDs_ = true;
  else hasFaceGIDs_ = true;
}
