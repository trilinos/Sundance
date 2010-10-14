#include "SundanceVertexSort.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <algorithm>

using namespace Sundance;
using namespace Teuchos;
using std::cout;
using std::cerr;
using std::endl;


int main(int argc, char** argv)
{
  int stat = 0;
  try
  {
    bool bad = false;
    GlobalMPISession session(&argc, &argv);

    Out::os() << endl;
    Array<int> tetVerts = tuple(234, 102, 371, 259);
    Array<int> triVerts = tuple(71, 43, 183);

    Array<Array<int> > testVerts = tuple(triVerts, tetVerts);
      
    for (int test=0; test<testVerts.size(); test++)
    {
      bool more = true;

      Array<int> exVerts = testVerts[test];
      int dim = exVerts.size()-1;

      while (more)
      {
        Out::os() << "exVerts = " << exVerts << endl;
        Array<int> ufcVerts = exVerts;
        int key = -1;
        vertexSort(ufcVerts, &key);
        Out::os() << "ufcVerts = " << ufcVerts << endl;

        for (int i=0; i<exVerts.size(); i++)
        {
          int u = exVertPosToUFCVertPos(dim, key, i);
          Out::os() << "ex vert=" << exVerts[i] << " is at pos=" 
                    << u << " in UFC array and has vertID=" << ufcVerts[u] << endl;
        }

        for (int f=0; f<dim+1; f++)
        {
          Out::os() << endl << "------------------------ " << endl;
          int missingEx = mapExSideToMissingVertex(dim, f);
          int missingUfc = exVertPosToUFCVertPos(dim, key, missingEx);
          Array<int> exFc = exSide(f, exVerts);
          Array<int> ufcFc = ufcSide(missingUfc, ufcVerts);
          Out::os() << "exodus face=" << f << ", missing vert pos="
                    << missingEx << ", id=" << exVerts[missingEx]<< endl;
          Out::os() << "verts = " << exFc << endl;
          Out::os() << "missing ufc vertex pos = "
                    << missingUfc << ", id=" << ufcVerts[missingUfc]
                    << endl;

          Out::os() << "ufc verts=" << ufcFc << endl;
          bool ok = true;
          insertionSort(exFc);
          for (int j=0; j<ufcFc.size(); j++)
          {
            if (exFc[j] != ufcFc[j]) ok = false;
          }
          if (ok) Out::os() << "check OK" << endl;
          else
          {
            Out::os() << "check FAILED: sort ex=" << exFc << ", ufc=" << ufcFc
                      << endl;
            bad = true;
          }
        }
        Out::os() << endl << endl;
        more = std::next_permutation(exVerts.begin(), exVerts.end());
      }
    }
    if (!bad) 
    {
      std::cerr << "all tests PASSED" << std::endl;
    }
    else
    {
      stat = -1;
      std::cerr << "a test has FAILED" << std::endl;
    }
  }
	catch(std::exception& e)
  {
    stat = -1;
    std::cerr << "detected exception " << e.what() << std::endl;
  }

  return stat;
}
