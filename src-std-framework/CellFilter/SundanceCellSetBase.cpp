/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCellSetBase.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

CellSetBase::CellSetBase(const Mesh& mesh, int cellDim,
                         const CellType& cellType)
  : mesh_(mesh), cellType_(cellType), dim_(cellDim)
{}

