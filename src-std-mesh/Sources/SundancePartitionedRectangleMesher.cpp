#include "SundancePartitionedRectangleMesher.hpp"
#include "SundanceOut.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;

Mesh PartitionedRectangleMesher::fillMesh() const
{
  SUNDANCE_OUT(verbosity() > VerbSilent,
               "PartitionedRectangleMesher::fillLocalMesh() is meshing "
               "rectangle [" << ax_ << ", " << bx_ << "] by ["
               << ay_ << ", " << by_ << "]");

  SUNDANCE_OUT(verbosity() == VerbHigh,
               "PartitionedRectangleMesher::fillLocalMesh() starting creation "
               "of empty mesh");

  Mesh mesh = createMesh(2);

  SUNDANCE_OUT(verbosity() == VerbHigh,
               "PartitionedRectangleMesher::fillLocalMesh() done creation of "
               "empty mesh");
  
  /* compute number of points per proc */

  int np = nProc();
	int rank = myRank();

	TEST_FOR_EXCEPTION(npx_ * npy_ != np, RuntimeError,
                     "PartitionedRectangleMesher::fillLocalMesh(): product "
                     "of npx=" << npx_ << " and npy=" << npy_
                     << " is not equal to np=" << np);

	/* compute number of points per proc */

	int nppx = nx_;
	int nppy = ny_;
  int nxTot = nx_*npx_;
  int nyTot = ny_*npy_;

	int px = rank/npy_;
	int py = rank % npy_;

	int lowestVisiblePtX = px*nppx-1;
	int lowestVisiblePtY = py*nppy-1;
	if (lowestVisiblePtX < 0) lowestVisiblePtX = 0;
	if (lowestVisiblePtY < 0) lowestVisiblePtY = 0;

	
	int highestVisiblePtX = lowestVisiblePtX + nppx + 1;
	int highestVisiblePtY = lowestVisiblePtY + nppy + 1;
	if (highestVisiblePtX > nxTot) highestVisiblePtX = nxTot;
	if (highestVisiblePtY > nyTot) highestVisiblePtY = nyTot;


	Array<Array<int> > pts(highestVisiblePtX-lowestVisiblePtX+1, 
												 highestVisiblePtY-lowestVisiblePtY+1);
	int globalIndex = 0;

	/* add the visible points into the mesh */
  for (int i=0; i<=nxTot; i++)
    {
      for (int j=0; j<=nyTot; j++, globalIndex++)
        {
					if (i < lowestVisiblePtX || i > highestVisiblePtX) continue;
					if (j < lowestVisiblePtY || j > highestVisiblePtY) continue;

					int ip = i/nppx;
					if (i==nxTot) ip--;
					int jp = j/nppy;
					if (j==nyTot) jp--;
					int pointOwner = ip*npy_ + jp;

          Point x( ax_ + ((double) i)*(bx_-ax_)/((double) nxTot) ,
                   ay_ + ((double) j)*(by_-ay_)/((double) nyTot));
          SUNDANCE_OUT(verbosity() > VerbLow, "adding point GID=" 
                       << globalIndex << " x=" << x << " owner=" 
                       << pointOwner);
          int lid = mesh.addVertex(globalIndex, x, pointOwner, 0);
          pts[i-lowestVisiblePtX][j-lowestVisiblePtY] = lid;
          
          SUNDANCE_OUT(verbosity() ==  VerbHigh,
                       "point " << x << " registered with LID=" << lid);
        }
    }

	/* add the visible cells to the mesh */
	globalIndex = 0 ;

  for (int i=0; i<nxTot; i++)
    {
      for (int j=0; j<nyTot; j++, globalIndex+=2)
				{
					if (i < lowestVisiblePtX || i >= highestVisiblePtX) continue;
					if (j < lowestVisiblePtY || j >= highestVisiblePtY) continue;

					int a = pts[i-lowestVisiblePtX][j-lowestVisiblePtY];
					int b = pts[i+1-lowestVisiblePtX][j-lowestVisiblePtY];
					int c = pts[i+1-lowestVisiblePtX][j+1-lowestVisiblePtY];
					int d = pts[i-lowestVisiblePtX][j+1-lowestVisiblePtY];

					int ip = i/nppx;
					int jp = j/nppy;
					int cellOwner = ip*npy_ + jp;
          Array<int> tri1;
          Array<int> tri2;
					if (i%2 == j%2)
						{
              tri1 = tuple(a,b,c);
              tri2 = tuple(a,c,d);
						}
					else
						{
              tri1 = tuple(a,b,d);
              tri2 = tuple(b,c,d);
						}
          int lid1 = mesh.addElement(globalIndex, tri1, cellOwner, 0);
          SUNDANCE_OUT(verbosity() ==  VerbHigh,
                       "elem " << tri1 
                       << " registered with LID=" << lid1);
          int lid2 = mesh.addElement(globalIndex+1, tri2, cellOwner, 0);
          SUNDANCE_OUT(verbosity() ==  VerbHigh,
                       "elem " << tri2 
                       << " registered with LID=" << lid2);
          
				}
    }
  return mesh;
}