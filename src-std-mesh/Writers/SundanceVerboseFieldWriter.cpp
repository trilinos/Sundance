/* @HEADER@ */
/* @HEADER@ */

#include "SundanceVerboseFieldWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"


using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;

void VerboseFieldWriter::write() const 
{
  int nProc = mesh().comm().getNProc();
  int myRank = mesh().comm().getRank();

  RefCountPtr<ostream> osp;
  if (filename().length()==0)
    {
      osp = rcp(&cout, false);
    }
  else 
    {
      string f = filename() + ".txt";
      if (nProc > 1) f = f + "." + Teuchos::toString(myRank);
      osp = rcp(new ofstream(f.c_str()));
    }
  ostream& os = *osp;

  if (myRank==0) os << "VerboseFieldWriter output" << endl;
  for (int p=0; p<nProc; p++)
    {
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      mesh().comm().synchronize();
      if (p != myRank) continue;
      os << "======== processor " << p << " ============================ "
         << endl;
      Tabs tab0;
      int dim = mesh().spatialDim();
      int nPts = mesh().numCells(0);
      int nElems = mesh().numCells(dim);
      os << tab0 << "spatial dimension = " << dim << endl;
      os << tab0 << "num points = " << nPts << endl;
      os << tab0 << "num elements = " << nElems << endl;
      os << tab0 << "Point list: " << endl;

      for (int i=0; i<nPts; i++)
        {
          Tabs tab1;
          os << tab1 << "L=" << i 
             << " G=" << mesh().mapLIDToGID(0, i) 
             << " x=" << mesh().nodePosition(i) 
             << " owner=" << mesh().ownerProcID(0,i) 
             << " label=" << mesh().label(0,i) << endl;
          int nc = mesh().numCofacets(0,i);
          Tabs tab2;
          os << tab2 << "num cofacets=" << nc << " cofs = {";
          for (int c=0; c<nc; c++)
            {
              if (c==0) os << " " ;
              else os << ", ";
              os << mesh().mapLIDToGID(dim, mesh().cofacetLID(0,i,c));
            }
          os << "}" << endl;
        }

      
      os << tab0 << "Element list: " << endl;

      for (int i=0; i<nElems; i++)
        {
          Tabs tab1;
          os << tab1 << "L=" << i 
             << " G=" << mesh().mapLIDToGID(dim, i) 
             << ", nodes L={";
          int numNodes = mesh().numFacets(dim, i, 0);
          for (int j=0; j<numNodes; j++)
            {
              if (j != 0) os << ", ";
              os << mesh().facetLID(dim, i, 0, j);
            }
          os << "}, G={";
          for (int j=0; j<numNodes; j++)
            {
              if (j != 0) os << ", ";
              os << mesh().mapLIDToGID(0, mesh().facetLID(dim, i, 0, j));
            }
          os << "}, owner=" << mesh().ownerProcID(dim,i)
             << ", label=" << mesh().label(dim,i) << endl;
          for (int fd=1; fd<dim; fd++)
            {
              Tabs tab2;
              os << tab2 << "facets of dimension " << fd << endl;
              int nf = mesh().numFacets(dim, i, fd);
              for (int f=0; f<nf; f++)
                {
                  Tabs tab3;
                  int flid = mesh().facetLID(dim, i, fd, f);
                  int fgid = -1;
                  int fowner = -1;
                  if (mesh().hasIntermediateGIDs(fd))
                    { 
                      fgid = mesh().mapLIDToGID(fd, flid);
                      fowner = mesh().ownerProcID(fd, flid);
                    }
                  os << tab3 << "f#=" << f << " L=" << flid
                     << " G=" << fgid << " owner=" << fowner
                     << " nodes={";
                  int nfn = mesh().numFacets(fd, flid, 0);
                  for (int fn=0; fn<nfn; fn++)
                    {
                      if (fn != 0) os << ", ";
                      os << mesh().facetLID(fd, flid, 0, fn);
                    }
                  os << "}" << endl;
                }
              
            }
        }
      
    }
}


