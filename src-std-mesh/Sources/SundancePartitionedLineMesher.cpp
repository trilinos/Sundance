#include "SundancePartitionedLineMesher.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;


PartitionedLineMesher::PartitionedLineMesher(const ParameterList& params)
  : MeshSourceBase(params),
    ax_(params.get<double>("ax")),
    bx_(params.get<double>("bx")),
    nx_(params.get<int>("nx"))
{;}

Mesh PartitionedLineMesher::fillMesh() const
{
  SUNDANCE_OUT(this->verbosity() > VerbSilent,
               "PartitionedLineMesher::fillLocalMesh() is meshing "
               "interval [" << ax_ << ", " << bx_ << "]");

  SUNDANCE_OUT(this->verbosity() == VerbHigh,
               "PartitionedLineMesher::fillLocalMesh() starting creation "
               "of empty mesh");

  Mesh mesh = createMesh(1);

  SUNDANCE_OUT(this->verbosity() == VerbHigh,
               "PartitionedLineMesher::fillLocalMesh() done creation of "
               "empty mesh");
  
  /* compute number of points per proc */

  int np = nProc();
	int nppx = nx_/np;

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
               "PartitionedLineMesher::fillLocalMesh() has " << nppx
               << " points per proc");

	int px = myRank();

	int lowestVisiblePtX = px*nppx-1;
	if (lowestVisiblePtX < 0) lowestVisiblePtX = 0;
	
	int highestVisiblePtX = lowestVisiblePtX + nppx + 1;
	if (highestVisiblePtX > nx_) highestVisiblePtX = nx_;

  SUNDANCE_OUT(this->verbosity() > VerbSilent,
               "index range is [" << lowestVisiblePtX << ", " << 
               highestVisiblePtX << "]");

	Array<int> pts(highestVisiblePtX-lowestVisiblePtX+1); 
	int globalIndex = 0;

	/* add the visible points into the mesh */
  for (int i=0; i<=nx_; i++, globalIndex++)
    {
			if (i < lowestVisiblePtX || i > highestVisiblePtX) continue;
			int pointOwner = i/nppx;
			if (i==nx_) pointOwner--;
			Point x( ax_ + ((double) i)*(bx_-ax_)/((double) nx_)); 

      SUNDANCE_OUT(this->verbosity() > VerbLow, "adding point GID=" 
                   << globalIndex << " x=" << x << " owner=" << pointOwner); 
			int lid = mesh.addVertex(globalIndex, x, pointOwner, 0);
			pts[i-lowestVisiblePtX] = globalIndex;
      SUNDANCE_OUT(this->verbosity() ==  VerbHigh,
                   "point " << x << " registered with LID=" << lid);
    }

	/* add the visible cells to the mesh */
	globalIndex = 0 ;

  for (int i=0; i<nx_; i++, globalIndex++)
    {
			if (i < lowestVisiblePtX || i >= highestVisiblePtX) continue;
			int a = pts[i-lowestVisiblePtX];
			int b = pts[i-lowestVisiblePtX+1];
			int cellOwner = i/nppx;
      SUNDANCE_OUT(this->verbosity() > VerbLow, "adding elem GID=" 
                   << globalIndex << " nodes=" << tuple(a,b) 
                   << " owner=" << cellOwner); 

      int lid = mesh.addElement(globalIndex, tuple(a,b), cellOwner, 0);
      SUNDANCE_OUT(this->verbosity() ==  VerbHigh,
                   "elem " << tuple(a,b) << " registered with LID=" << lid);
    }
  
  return mesh;
}
