/* @HEADER@ */
/* @HEADER@ */

#include "SundanceElementIntegral.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;

ElementIntegral::ElementIntegral(int dim, 
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 int testDerivOrder)
  : dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(nNodesTest_),
    isTwoForm_(false)
{;}



ElementIntegral::ElementIntegral(int dim,
                         const CellType& cellType,
                         const BasisFamily& testBasis,
                         int testDerivOrder,
                         const BasisFamily& unkBasis,
                         int unkDerivOrder)
  : dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)), 
    unkDerivOrder_(unkDerivOrder), 
    nRefDerivUnk_(ipow(dim, unkDerivOrder)),
    nNodesUnk_(unkBasis.nNodes(cellType)), 
    nNodes_(nNodesTest_*nNodesUnk_),
    isTwoForm_(true)
{;}

int ElementIntegral::ipow(int base, int power) 
{
  int rtn = 1;
  for (int i=0; i<power; i++) rtn *= base;
  return rtn;
}

