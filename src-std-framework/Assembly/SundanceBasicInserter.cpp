/* @HEADER@ */
/* @HEADER@ */

#include "SundanceBasicInserter.hpp"
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

void BasicInserter::insert(int cellDim, const RefCountPtr<Array<int> >& cells,
                           bool isBC,
                           const LocalMatrixContainer* localMat,
                           LinearOperator<double>& A,
                           Vector<double>& b) const
{}


void BasicInserter::configureMat(LinearOperator<double>& A) const
{
}

