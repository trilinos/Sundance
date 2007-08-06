#include "SundanceExodusMeshReader.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"

#ifdef HAVE_EXODUS
#include "exodusII.h"
#endif 

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;


ExodusMeshReader::ExodusMeshReader(const string& fname,
                                               const MeshType& meshType,
                                               const MPIComm& comm)
  : MeshReaderBase(fname, meshType, comm)
{
  TEST_FOR_EXCEPTION(nProc() > 1, RuntimeError,
                     "ExodusMeshReader not useable with parallel meshes");
}

ExodusMeshReader::ExodusMeshReader(const ParameterList& params)
  : MeshReaderBase(params)
{
  TEST_FOR_EXCEPTION(nProc() > 1, RuntimeError,
                     "ExodusMeshReader not useable with parallel meshes");
}




Mesh ExodusMeshReader::fillMesh() const 
{
  Mesh mesh;
#ifndef HAVE_EXODUS
  TEST_FOR_EXCEPTION(true, RuntimeError, 
    "ExodusMeshReader called for a build without ExodusII");
#else

  int CPU_word_size = 8;
  int IO_word_size = 0;
  float version;

  if (verbosity() > VerbMedium) ex_opts(EX_DEBUG | EX_VERBOSE);

  int exoID = ex_open(filename().c_str(), EX_READ, 
    &CPU_word_size, &IO_word_size, &version);

  TEST_FOR_EXCEPTION(exoID < 0, RuntimeError, "ExodusMeshReader unable to "
    "open file: " << filename());

  TEST_FOR_EXCEPT(IO_word_size != 8 || CPU_word_size != 8);

  char title[MAX_LINE_LENGTH+1];

  int dim = 0;
  int numNodes = 0;
  int numElems = 0;
  int numElemBlocks = 0;
  int numNodeSets = 0;
  int numSideSets = 0;

  int ierr = ex_get_init(exoID, title, &dim, &numNodes, &numElems,
    &numElemBlocks, &numNodeSets, &numSideSets);

  TEST_FOR_EXCEPTION(numNodes <= 0, RuntimeError, "invalid numNodes=" 
    << numNodes);
  TEST_FOR_EXCEPTION(numElems <= 0, RuntimeError, "invalid numElems=" 
    << numElems);

  /* now we can build the mesh */
  mesh = createMesh(dim);


  /* Read the points */
  Array<double> x(numNodes);
  Array<double> y(numNodes);
  Array<double> z(numNodes * (dim > 2));

  if (dim == 2)
  {
    ierr = ex_get_coord(exoID, &(x[0]), &(y[0]), (void*) 0);
    TEST_FOR_EXCEPT(ierr < 0);
    TEST_FOR_EXCEPT(ierr > 0);
  }
  else if (dim==3)
  {
    ierr = ex_get_coord(exoID, (void*) &(x[0]), (void*) &(y[0]), (void*) &(z[0]));
    TEST_FOR_EXCEPT(ierr < 0);
  }
  else 
  {
    TEST_FOR_EXCEPTION(dim < 2 || dim > 3, RuntimeError, 
      "invalid dimension=" << dim << " in ExodusMeshReader");
  }

  /* add the points to the mesh */
  for (int n=0; n<numNodes; n++)
    {
      Point p;
      if (dim==2)
        {
          p = Point(x[n], y[n]);
        }
      else
        {
          p = Point(x[n], y[n], z[n]);
        }
      mesh.addVertex(n, p, 0, 0);
    }

  /* Read the elements for each block */
  int lid=0;
  Array<int> blockIDs(numElemBlocks);
  if (numElemBlocks > 0)
  {
    ierr = ex_get_elem_blk_ids(exoID, &(blockIDs[0]));
  }
  TEST_FOR_EXCEPT(ierr < 0);
  for (int b=0; b<numElemBlocks; b++)
  {
    char elemType[MAX_LINE_LENGTH+1];
    int elsInBlock;
    int nodesPerEl;
    int numAttrs;
    int bid = blockIDs[b];

    ierr = ex_get_elem_block(exoID, bid, elemType, &elsInBlock,
      &nodesPerEl, &numAttrs);
    TEST_FOR_EXCEPT(ierr < 0);

    Array<int> connect(elsInBlock * nodesPerEl);

    ierr = ex_get_elem_conn(exoID, bid, &(connect[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    int n=0;

    for (int e=0; e<elsInBlock; e++, n+=nodesPerEl)
    {
      if (dim==2)
      {
        mesh.addElement(lid, tuple(connect[n]-1, connect[n+1]-1, connect[n+2]-1), 0, bid);
        SUNDANCE_VERB_HIGH("adding element=("
          << connect[n]-1 << ", " << connect[n+1]-1
          << ", " << connect[n+2]-1 << ")");
      }
      else
      {
        mesh.addElement(lid, 
          tuple(connect[n]-1, connect[n+1]-1, connect[n+2]-1, connect[n+3]-1),
          0, bid);
      }
    }
  }


  /* Read the node sets */
  Array<int> nsIDs(numNodeSets);
  if (numNodeSets > 0)
  {
    ierr = ex_get_node_set_ids(exoID, &(nsIDs[0]));
  }
  TEST_FOR_EXCEPT(ierr < 0);
  for (int ns=0; ns<numNodeSets; ns++)
  {
    int nNodes;
    int nDist;
    int nsID = nsIDs[ns];
    ierr = ex_get_node_set_param(exoID, nsID, &nNodes, &nDist);
    TEST_FOR_EXCEPT(ierr < 0);
    Array<int> nodes(nNodes);
    ierr = ex_get_node_set(exoID, nsID, &(nodes[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    for (int n=0; n<nNodes; n++)
    {
      mesh.setLabel(0, nodes[n]-1, nsID);
    }
  }

    
  /* Read the side sets */
  Array<int> ssIDs(numSideSets);
  if (numSideSets > 0)
  {
    ierr = ex_get_side_set_ids(exoID, &(ssIDs[0]));
  }
  TEST_FOR_EXCEPT(ierr < 0);
  for (int ss=0; ss<numSideSets; ss++)
  {
    int nSides;
    int nDist;
    int ssID = ssIDs[ss];
    ierr = ex_get_side_set_param(exoID, ssID, &nSides, &nDist);
    TEST_FOR_EXCEPT(ierr < 0);
    Array<int> sides(nSides);
    Array<int> elems(nSides);
    ierr = ex_get_side_set(exoID, ssID, &(elems[0]), &(sides[0]));
    TEST_FOR_EXCEPT(ierr < 0);
    for (int n=0; n<nSides; n++)
    {
      int elemID = elems[n];
      int facetNum = sides[n];
      int fsign;
      int sideLID = mesh.facetLID(dim, elemID-1, dim-1, facetNum-1,fsign);
      mesh.setLabel(dim-1, sideLID, ssID);
    }
  }

  ierr = ex_close(exoID);
  TEST_FOR_EXCEPT(ierr < 0);
#endif
	return mesh;
}

