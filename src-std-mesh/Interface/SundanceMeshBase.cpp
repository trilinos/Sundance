/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMeshBase.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExceptions.hpp"

using namespace Sundance::StdMesh::Internal;
using namespace Sundance;
using namespace Teuchos;


MeshBase::MeshBase(int dim, const MPIComm& comm) 
  : dim_(dim), 
    comm_(comm),
    reorderer_(Mesh::defaultReorderer().createInstance(this)) 
{;}



void MeshBase::getFacetArray(int cellDim, int cellLID, int facetDim, 
                             Array<int>& facetLIDs) const
{
  int nf = numFacets(cellDim, cellLID, facetDim);
  facetLIDs.resize(nf);
  for (int f=0; f<nf; f++) 
    {
      facetLIDs[f] = facetLID(cellDim, cellLID, facetDim, f);
    }
}



