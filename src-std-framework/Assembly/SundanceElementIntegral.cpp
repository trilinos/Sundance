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
                                 const CellType& cellType)
  : dim_(dim),
    testDerivOrder_(-1), 
    nRefDerivTest_(-1),
    nNodesTest_(-1),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(-1),
    order_(1),
    alpha_(),
    beta_()
{;}

ElementIntegral::ElementIntegral(int dim, 
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 const Array<int>& alpha,
                                 int testDerivOrder)
  : dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(nNodesTest_),
    order_(1),
    alpha_(alpha),
    beta_()
{;}



ElementIntegral::ElementIntegral(int dim,
                                 const CellType& cellType,
                                 const BasisFamily& testBasis,
                                 const Array<int>& alpha,
                                 int testDerivOrder,
                                 const BasisFamily& unkBasis,
                                 const Array<int>& beta,
                                 int unkDerivOrder)
  : dim_(dim),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(dim, testDerivOrder)),
    nNodesTest_(testBasis.nNodes(cellType)), 
    unkDerivOrder_(unkDerivOrder), 
    nRefDerivUnk_(ipow(dim, unkDerivOrder)),
    nNodesUnk_(unkBasis.nNodes(cellType)), 
    nNodes_(nNodesTest_*nNodesUnk_),
    order_(2),
    alpha_(alpha),
    beta_(beta)
{;}

int ElementIntegral::ipow(int base, int power) 
{
  int rtn = 1;
  for (int i=0; i<power; i++) rtn *= base;
  return rtn;
}

