/* @HEADER@ */
/* @HEADER@ */

#include "SundanceMeshBase.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExceptions.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace SundanceStdMesh::Internal;
using namespace SundanceStdMesh;
using namespace Teuchos;
using namespace SundanceUtils;


MeshBase::MeshBase(int dim, const MPIComm& comm) 
  : dim_(dim), 
    comm_(comm),
    reorderer_(Mesh::defaultReorderer().createInstance(this)) 
{;}


static Time& facetGrabTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("facet grabbing"); 
  return *rtn;
}

void MeshBase::getFacetArray(int cellDim, int cellLID, int facetDim, 
                             Array<int>& facetLIDs) const
{
  TimeMonitor timer(facetGrabTimer());

  int nf = numFacets(cellDim, cellLID, facetDim);
  facetLIDs.resize(nf);
  for (int f=0; f<nf; f++) 
    {
      facetLIDs[f] = facetLID(cellDim, cellLID, facetDim, f);
    }
}



