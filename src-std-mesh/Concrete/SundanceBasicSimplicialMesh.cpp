/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBasicSimplicialMesh.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceDebug.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

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


//#define SKIP_FACES

BasicSimplicialMesh::BasicSimplicialMesh(int dim, const MPIComm& comm)
	: CreatableMesh(dim, comm),
    numCells_(dim+1),
    points_(),
    edgeVerts_(2),
    faceVerts_(dim),
    faceEdges_(dim),
    elemVerts_(dim+1),
    elemEdges_(),
    elemFaces_(dim+1),
    elemFaceRotations_(dim+1),
    vertexSetToFaceIndexMap_(),
    edgeCofacets_(),
    faceCofacets_(),
    vertEdges_(),
    vertCofacets_(),
    vertEdgePartners_(),
    vertEdgeSigns_(),
    tmpFaceVerts_(),
    tmpFaceEdges_(),
    LIDToGIDMap_(dim+1),
    GIDToLIDMap_(dim+1),
    ownerProcID_(dim+1),
    labels_(dim+1),
    base_(1),
    tmpBase_(1),
    hasEdgeGIDs_(false),
    hasFaceGIDs_(false)
{
  verbosity() = MeshSource::classVerbosity();
  estimateNumVertices(1000);
  estimateNumElements(1000);

  base_[0] = &(faceVerts_.value(0,0));

  tmpFaceVerts_.resize(3);
  tmpBase_[0] = &(tmpFaceVerts_[0]);

  tmpFaceEdges_.resize(3);
  tmpBase_[0] = &(tmpFaceVerts_[0]);

  if (spatialDim()==2) elemEdges_.setTupleSize(3);
  if (spatialDim()==3) elemEdges_.setTupleSize(6);
}


void BasicSimplicialMesh::getJacobians(int cellDim, const Array<int>& cellLID,
                                       CellJacobianBatch& jBatch) const
{

  
  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
                     "cellDim=" << cellDim 
                     << " is not in expected range [0, " << spatialDim()
                     << "]");

  jBatch.resize(cellLID.size(), spatialDim(), cellDim);

  if (cellDim < spatialDim())
    {
      for (int i=0; i<cellLID.size(); i++)
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
                int a = faceVerts_.value(lid, 0);
                int b = faceVerts_.value(lid, 1);
                int c = faceVerts_.value(lid, 2);
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
  

      for (int i=0; i<cellLID.size(); i++)
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
}

void BasicSimplicialMesh::pushForward(int cellDim, const Array<int>& cellLID,
                                      const Array<Point>& refQuadPts,
                                      Array<Point>& physQuadPts) const
{
  TEST_FOR_EXCEPTION(cellDim < 0 || cellDim > spatialDim(), InternalError,
                     "cellDim=" << cellDim 
                     << " is not in expected range [0, " << spatialDim()
                     << "]");

  int nQuad = refQuadPts.size();
  Array<double> J(cellDim*cellDim);

  if (physQuadPts.size() > 0) physQuadPts.resize(0);
  physQuadPts.reserve(cellLID.size() * refQuadPts.size());


  for (int i=0; i<cellLID.size(); i++)
    {
      int lid = cellLID[i];
      switch(cellDim)
        {
        case 0:
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
                a = faceVerts_.value(lid, 0);
                b = faceVerts_.value(lid, 1);
                c = faceVerts_.value(lid, 2);
              }
            const Point& pa = points_[a];
            const Point& pb = points_[b];
            const Point& pc = points_[c];
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
          }
          break;
        default:
          TEST_FOR_EXCEPTION(true, InternalError, "impossible switch value "
                             "in BasicSimplicialMesh::getJacobians()");
        }
    }
}

void BasicSimplicialMesh::estimateNumVertices(int nPts)
{
  points_.reserve(nPts);
  vertCofacets_.reserve(nPts);
  vertEdges_.reserve(nPts);
  vertEdgePartners_.reserve(nPts);
  vertEdgeSigns_.reserve(nPts);

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
    }
  else if (spatialDim()==2)
    {
      nEdges = 3*nElems;
    }
  

  vertexSetToFaceIndexMap_ = Hashtable<VertexSet, int>(nFaces);

  edgeVerts_.reserve(nEdges);
  faceVerts_.reserve(nFaces);
  faceEdges_.reserve(nFaces);
  elemVerts_.reserve(nElems);
  elemEdges_.reserve(nElems);
  elemFaces_.reserve(nElems);
  edgeCofacets_.reserve(nEdges);
  faceCofacets_.reserve(nFaces);

  ownerProcID_[spatialDim()].reserve(nElems);
  LIDToGIDMap_[spatialDim()].reserve(nElems);
  GIDToLIDMap_[spatialDim()] = Hashtable<int,int>(nElems, 0.6);
  labels_[spatialDim()].reserve(nElems);

  base_[0] = &(faceVerts_.value(0,0));
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
                                       Array<int>& facetLID) const 
{
  TimeMonitor timer(batchedFacetGrabTimer());

  int nf = numFacets(cellDim, cellLID[0], facetDim);
  facetLID.resize(cellLID.size() * nf);

  
  if (facetDim==0)
    {
      if (cellDim == spatialDim())
        {
          int ptr=0;
          for (int c=0; c<cellLID.size(); c++)
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
          for (int c=0; c<cellLID.size(); c++)
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
          for (int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(faceVerts_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                }
            }
        }
    }
  else if (facetDim==1)
    {
      if (cellDim == spatialDim())
        {
          int ptr=0;
          for (int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(elemEdges_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                }
            }
        }
      else
        {
          int ptr=0;
          for (int c=0; c<cellLID.size(); c++)
            {
              const int* fPtr = &(faceEdges_.value(cellLID[c], 0));
              for (int f=0; f<nf; f++, ptr++)
                {
                  facetLID[ptr] = fPtr[f];
                }
            }
        }
    }
  else
    {
      int ptr=0;
      for (int c=0; c<cellLID.size(); c++)
        {
          const int* fPtr = &(elemFaces_.value(cellLID[c], 0));
          for (int f=0; f<nf; f++, ptr++)
            {
              facetLID[ptr] = fPtr[f];
            }
        }
    }
}

int BasicSimplicialMesh::facetLID(int cellDim, int cellLID,
                                  int facetDim, int facetIndex) const 
{

  if (facetDim==0)
    {
      if (cellDim == spatialDim())
        {
          return elemVerts_.value(cellLID, facetIndex);
        }
      else if (cellDim==1) return edgeVerts_.value(cellLID, facetIndex);
      else if (cellDim==2) return faceVerts_.value(cellLID, facetIndex);
    }
  if (facetDim==1)
    {
      if (cellDim==spatialDim())
        {
          return elemEdges_.value(cellLID, facetIndex);
        }
      else
        {
          return faceEdges_.value(cellLID, facetIndex);
        }
    }
  else
    {
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

  SUNDANCE_OUT(verbosity() > VerbLow,
               "BSM points are " << points_.toString());

  int lid = points_.length();
  points_.append(x);

  SUNDANCE_OUT(verbosity() > VerbLow,
               "BSM added point " << x << " lid = " << lid);

  numCells_[0]++;

  LIDToGIDMap_[0].append(globalIndex);
  GIDToLIDMap_[0].put(globalIndex, lid);

  ownerProcID_[0].append(ownerProcID);
  labels_[0].append(label);

  vertCofacets_.resize(points_.length());
  vertEdges_.resize(points_.length());
  vertEdgePartners_.resize(points_.length());
  vertEdgeSigns_.resize(points_.length());
  
  return lid;
}

int BasicSimplicialMesh::addElement(int globalIndex, 
                                    const Array<int>& vertLID,
                                    int ownerProcID, int label)
{
  int lid = elemVerts_.length();

  numCells_[spatialDim()]++;

  LIDToGIDMap_[spatialDim()].append(globalIndex);
  GIDToLIDMap_[spatialDim()].put(globalIndex, lid);
  labels_[spatialDim()].append(label);
  ownerProcID_[spatialDim()].append(ownerProcID);

  

  static Array<int> edges;
  static Array<int> faces;
  static Array<int> faceRotations;
  
  if (spatialDim()==1)
    {
      edges.resize(0);
      faces.resize(0);
      vertCofacets_[vertLID[0]].append(lid);
      vertCofacets_[vertLID[1]].append(lid);
    }
  if (spatialDim()==2)
    {
      edges.resize(3);
      faces.resize(0);
      edges[0] = addEdge(vertLID[0], vertLID[1]);
      edges[1] = addEdge(vertLID[1], vertLID[2]);
      edges[2] = addEdge(vertLID[2], vertLID[0]);
      vertCofacets_[vertLID[0]].append(lid);
      vertCofacets_[vertLID[1]].append(lid);
      vertCofacets_[vertLID[2]].append(lid);

    }
  else if (spatialDim()==3)
    {
      edges.resize(6);
      faces.resize(4);
      faceRotations.resize(4);
      edges[0] = addEdge(vertLID[0], vertLID[1]);
      edges[1] = addEdge(vertLID[1], vertLID[2]);
      edges[2] = addEdge(vertLID[2], vertLID[0]);
      edges[3] = addEdge(vertLID[0], vertLID[3]);
      edges[4] = addEdge(vertLID[1], vertLID[3]);
      edges[5] = addEdge(vertLID[2], vertLID[3]);
      vertCofacets_[vertLID[0]].append(lid);
      vertCofacets_[vertLID[1]].append(lid);
      vertCofacets_[vertLID[2]].append(lid);
      vertCofacets_[vertLID[3]].append(lid);
#ifndef SKIP_FACES
      //setWorkFace(vertLID[0], vertLID[1], vertLID[3]);
      faces[0] = addFace(vertLID[0], vertLID[1], vertLID[3], 
                         edges[0], edges[4], edges[3],
                         faceRotations[0]);

      //setWorkFace(vertLID[1], vertLID[2], vertLID[3]);
      faces[1] = addFace(vertLID[1], vertLID[2], vertLID[3], 
                         edges[1], edges[5], edges[4],
                         faceRotations[1]);
      // setWorkFace(vertLID[0], vertLID[3], vertLID[2]);
      faces[2] = addFace(vertLID[0], vertLID[3], vertLID[2], 
                         edges[3], edges[5], edges[2],
                         faceRotations[2]);
      //  setWorkFace(vertLID[2], vertLID[1], vertLID[0]);
      faces[3] = addFace(vertLID[2], vertLID[1], vertLID[0], 
                         edges[1], edges[0], edges[2],
                         faceRotations[3]);
#endif
    }

  elemVerts_.append(vertLID);
  if (edges.length() > 0) elemEdges_.append(edges);

#ifndef SKIP_FACES
  if (faces.length() > 0) 
    {
      elemFaces_.append(faces);
      elemFaceRotations_.append(faceRotations);
    }
#endif
  for (int i=0; i<edges.length(); i++) edgeCofacets_[edges[i]].append(lid);

#ifndef SKIP_FACES
  for (int i=0; i<faces.length(); i++) faceCofacets_[faces[i]].append(lid);
#endif  
  return lid;
}

int BasicSimplicialMesh::addFace(int v1, int v2, int v3, 
                                 int e1, int e2, int e3,
                                 int& rotation)
{
  int rtn = getFaceLIDFromVertLIDs(v1, v2, v3, rotation);
  if (rtn >= 0) 
    {
      return rtn;
    }

  rotation = 0 ;
  rtn = faceVerts_.length();

  faceVerts_.append(&(tmpFaceVerts_[0]), tmpFaceVerts_.length());

  tmpFaceEdges_[0] = e1;
  tmpFaceEdges_[1] = e2;
  tmpFaceEdges_[2] = e3;
  faceEdges_.append(&(tmpFaceEdges_[0]), tmpFaceEdges_.length());
  
  base_[0] = &(faceVerts_.value(0,0));

  VertexSet face(&(base_[0]), rtn, tmpFaceVerts_.length());

  vertexSetToFaceIndexMap_.put(face, rtn);
  
  faceCofacets_.resize(rtn+1);

  numCells_[spatialDim()-1]++;

  return rtn;
}

int BasicSimplicialMesh::getEdgeLIDFromVertLIDs(int v1, int v2) 
{
  SUNDANCE_OUT(verbosity() > VerbHigh, 
               "p=" << comm().getRank() << " finding edge for verts="
               << v1 << ", " << v2);
  const Array<int>& edgePartners = vertEdgePartners_[v1];
  for (int i=0; i<edgePartners.length(); i++)
    {
      SUNDANCE_OUT(verbosity() > VerbHigh, 
                   "p=" << comm().getRank() << " checking partner="
                   << edgePartners[i]);
      if (edgePartners[i] == v2)
        {
          SUNDANCE_OUT(verbosity() > VerbHigh, "p=" << comm().getRank() 
                       << " found match!");
          return vertEdges_[v1][i];
        }
    }
  SUNDANCE_OUT(verbosity() > VerbHigh, "p=" << comm().getRank() 
               << " no match!");
  return -1;
}

int BasicSimplicialMesh::getFaceLIDFromVertLIDs(int v1, int v2, int v3,
                                                int& rotation) 
{
  tmpFaceVerts_[0] = v1;
  tmpFaceVerts_[1] = v2;
  tmpFaceVerts_[2] = v3;

  VertexSet tmpFace(&(tmpBase_[0]), 0, tmpFaceVerts_.length());

  if (vertexSetToFaceIndexMap_.containsKey(tmpFace))
    {
      rotation = tmpFace.rotation();
      return vertexSetToFaceIndexMap_.get(tmpFace);
    }
  rotation = 0;
  return -1;
}


int BasicSimplicialMesh::addEdge(int v1, int v2)
{
  int rtn = getEdgeLIDFromVertLIDs(v1, v2);
  if (rtn >= 0) return rtn;
  
  rtn = edgeVerts_.length();
  edgeVerts_.append(tuple(v1,v2));

  vertEdges_[v1].append(rtn);
  vertEdgePartners_[v1].append(v2);
  vertEdgeSigns_[v1].append(-1);

  vertEdges_[v2].append(rtn);
  vertEdgePartners_[v2].append(v1);
  vertEdgeSigns_[v2].append(-1);

  edgeCofacets_.resize(rtn+1);

  numCells_[1]++;

  return rtn;
}


void BasicSimplicialMesh::synchronizeGIDNumbering(int dim, int localCount) 
{
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  SUNDANCE_OUT(verbosity() > VerbMedium, 
               "sharing offsets for GID numbering for dim=" << dim);
  SUNDANCE_OUT(verbosity() > VerbMedium, 
               "I have " << localCount << " cells");
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  Array<int> gidOffsets;
  int total;
  MPIContainerComm<int>::accumulate(localCount, gidOffsets, total, comm());
  int myOffset = gidOffsets[comm().getRank()];

  SUNDANCE_OUT(verbosity() > VerbMedium, 
               "back from MPI accumulate");

  for (int i=0; i<numCells(dim); i++)
    {
      if (LIDToGIDMap_[dim][i] == -1) continue;
      LIDToGIDMap_[dim][i] += myOffset;
      GIDToLIDMap_[dim].put(LIDToGIDMap_[dim][i], i);
    }
  SUNDANCE_OUT(verbosity() > VerbMedium, 
               "done sharing offsets for GID numbering for dim=" << dim);
}


void BasicSimplicialMesh::assignIntermediateCellOwners(int cellDim)
{
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  comm().synchronize();
  SUNDANCE_OUT(verbosity() > VerbMedium,  
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

  /* allocate a temporary array of vertices we'll use lots of times */
  Array<int> verts(cellDim+1);

  /* Define these so we're not constantly dereferencing the owners arrays */
  const Array<int>& ptOwners = ownerProcID_[0];
  Array<int>& cellOwners = ownerProcID_[cellDim];
  cellOwners.resize(numCells(cellDim));

  /* run through the cells, assigning GID numbers to the locally owned ones
   * and compiling a list of the remotely owned ones. */
  int localCount = 0;
  ArrayOfTuples<int>* cellVerts;
  if (cellDim==2) cellVerts = &(faceVerts_);
  else cellVerts = &(edgeVerts_);

  LIDToGIDMap_[cellDim].resize(numCells(cellDim));
  
  SUNDANCE_OUT(verbosity() > VerbMedium,  
               "starting loop over cells");

  for (int i=0; i<numCells(cellDim); i++)
    {  
      SUNDANCE_OUT(verbosity() > VerbHigh, 
                   "p=" << myRank 
                   <<" cell " << i);
      /* We'll define the owner of an intermediate cell as
       * the max rank of the cell's point facets. With this
       * definition, all processors will agree on who owns each cell. */
      int owner = 0;
      for (int v=0; v<=cellDim; v++)
        {
          verts[v] = cellVerts->value(i, v);
          if (ptOwners[verts[v]] > owner) owner = ptOwners[verts[v]];
        }
      SUNDANCE_OUT(verbosity() > VerbHigh, 
                   "p=" << myRank 
                   << " verts of cell " << i << " are " << verts 
                   << " and its owner is " << owner);
      /* record the owner we just determined */
      cellOwners[i] = owner;

      /* If the cell is locally owned, give it a temporrary GID. This GID
       * is temporary because it's one of a sequence numbered from zero
       * rather than the lowest GID on this processor. We'll offset
       * all of these GIDs once we know how many cells this processor owns. */
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
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       "p=" << myRank << " adding LID=" << i 
                       << " to req list");
          outgoingRequestLIDs[owner].append(i);
          for (int v=0; v<=cellDim; v++)
            {
              int ptGID = LIDToGIDMap_[0][verts[v]];
              SUNDANCE_OUT(verbosity() > VerbHigh, 
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
  SUNDANCE_OUT(verbosity() > VerbHigh,
               "p=" << myRank << "sending requests: " << outgoingRequests);
  MPIContainerComm<int>::allToAll(outgoingRequests,
                                  incomingRequests,
                                  comm());
  SUNDANCE_OUT(verbosity() > VerbHigh,
               "p=" << myRank << "recv'd requests: " << incomingRequests);

  /* Answer the requests incoming to this processor */
  int dummyRotation;
  for (int p=0; p<comm().getNProc(); p++)
    {
      const Array<int>& requestsFromProc = incomingRequests[p];
      
      for (int c=0; c<requestsFromProc.size(); c+=(cellDim+1))
        {
          SUNDANCE_OUT(verbosity() > VerbHigh,
                       "p=" << myRank << "processing request c=" << c/(cellDim+1));
          /* The request for each cell is a tuple of vertex GIDs. 
           * We convert these to LIDs and pack in our verts[] helper array */
          for (int v=0; v<=cellDim; v++)
            {
              verts[v] = mapGIDToLID(0, requestsFromProc[c+v]);
            }
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       "p=" << myRank << "point LIDs are " << verts);
          /* Now we can identify the LID of the cell as known on 
           * this processor */
          int cellLID;
          if (cellDim==1) 
            {
              cellLID = getEdgeLIDFromVertLIDs(verts[0], verts[1]);
            }
          else
            {
              cellLID = getFaceLIDFromVertLIDs(verts[0], verts[1], verts[2],
                                               dummyRotation);
            }
          SUNDANCE_OUT(verbosity() > VerbHigh, 
                       "p=" << myRank << "cell LID is " << cellLID);
          /* Finally, we get the cell's GID and append to the list 
           * of GIDs to send back as answers. */
          outgoingGIDs[p].append(mapLIDToGID(cellDim, cellLID));
        }
    }


  /* We now share the GIDs for the requested cells */
  SUNDANCE_OUT(verbosity() > VerbHigh,
               "p=" << myRank << "sending GIDs: " << outgoingGIDs);
  MPIContainerComm<int>::allToAll(outgoingGIDs,
                                  incomingGIDs,
                                  comm());
  SUNDANCE_OUT(verbosity() > VerbHigh,
               "p=" << myRank << "recv'ing GIDs: " << incomingGIDs);

  /* Now assign the cell GIDs we've recv'd from from the other procs */
  for (int p=0; p<comm().getNProc(); p++)
    {
      const Array<int>& gidsFromProc = incomingGIDs[p];
      for (int c=0; c<gidsFromProc.size(); c++)
        {
          int cellLID = outgoingRequestLIDs[p][c];
          int cellGID = incomingGIDs[p][c];
          SUNDANCE_OUT(verbosity() > VerbHigh,
                       "p=" << myRank << 
                       " assigning GID=" << cellGID << " to LID=" 
                       << cellLID);
          LIDToGIDMap_[cellDim][cellLID] = cellGID;
          GIDToLIDMap_[cellDim].put(cellGID, cellLID);
        }
    }

  SUNDANCE_OUT(verbosity() > VerbMedium,  
               "p=" << myRank << "done synchronizing cells of dimension " << cellDim);

  if (cellDim==1) hasEdgeGIDs_ = true;
  else hasFaceGIDs_ = true;
}
