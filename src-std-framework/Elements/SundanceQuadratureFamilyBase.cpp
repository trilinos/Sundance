/* @HEADER@ */
/* @HEADER@ */

#include "SundanceQuadratureFamilyBase.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;

void QuadratureFamilyBase::getPoints(const CellType& cellType, 
                                     Array<Point>& quadPoints,
                                     Array<double>& quadWeights) const
{
  switch(cellType)
    {
    case PointCell:
      quadPoints = tuple(Point());
      quadWeights = tuple(1.0);
      break;
    case LineCell:
      getLineRule(quadPoints, quadWeights);
      break;
    case TriangleCell:
      getTriangleRule(quadPoints, quadWeights);
      break;
    case QuadCell:
      getQuadRule(quadPoints, quadWeights);
      break;
    case TetCell:
      getTetRule(quadPoints, quadWeights);
      break;
    case BrickCell:
      getBrickRule(quadPoints, quadWeights);
      break;
    default:
      SUNDANCE_ERROR("cell type " << cellType << " not handled in "
                     "QuadratureFamilyBase::getPoints()");
    }
}

void QuadratureFamilyBase::getLineRule(Array<Point>& /* quadPoints */,
                                       Array<double>& /* quadWeights */) const 
{
  SUNDANCE_ERROR("Line rule not available for " << toXML());
}

void QuadratureFamilyBase::getTriangleRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
  SUNDANCE_ERROR("Triangle rule not available for " << toXML());
}

void QuadratureFamilyBase::getQuadRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
  SUNDANCE_ERROR("Quad cell rule not available for " << toXML());
}

void QuadratureFamilyBase::getTetRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
  SUNDANCE_ERROR("Tet cell rule not available for " << toXML());
}

void QuadratureFamilyBase::getBrickRule(Array<Point>& /* quadPoints */,
                                           Array<double>& /* quadWeights */) const 
{
  SUNDANCE_ERROR("Brick rule not available for " << toXML());
}
