/* @HEADER@ */
/* @HEADER@ */

#include "SundanceAssembler.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceHomogeneousDOFMap.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"


using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;


Assembler::Assembler(const Mesh& mesh, 
                     const RefCountPtr<EquationSet>& eqn)
  : mesh_(mesh),
    eqn_(eqn),
    rowMap_(),
    colMap_(),
    bcRows_(),
    rqc_(),
    isBCRqc_(),
    rqcExprs_(),
    rqcDerivSet_()
{
  DOFMapBuilder mapBuilder(mesh, eqn);
  rowMap_ = mapBuilder.rowMap();
  colMap_ = mapBuilder.colMap();
  bcRows_ = mapBuilder.bcRows();

  for (int r=0; r<eqn_->regionQuadCombos().size(); r++)
    {
      rqc_.append(eqn_->regionQuadCombos()[r]);
      isBCRqc_.append(false);
      rqcExprs_.append(eqn_->expr(eqn_->regionQuadCombos()[r]));
      rqcDerivSet_.append(eqn_->nonzeroFunctionalDerivs(eqn_->regionQuadCombos()[r]));
    }
  for (int r=0; r<eqn_->bcRegionQuadCombos().size(); r++)
    {
      rqc_.append(eqn_->bcRegionQuadCombos()[r]);
      isBCRqc_.append(true);
      rqcExprs_.append(eqn_->bcExpr(eqn_->bcRegionQuadCombos()[r]));
      rqcDerivSet_.append(eqn_->nonzeroBCFunctionalDerivs(eqn_->bcRegionQuadCombos()[r]));
    }
}

void Assembler::print(ostream& os) const 
{
  for (int r=0; r<rqc_.size(); r++)
    {
      os << "Region/Quad combination: " << rqc_[r] << endl;
      os << "isBC = " << Teuchos::toString(isBCRqc_[r]) << endl;
      os << "nonzero derivs = " << rqcDerivSet_[r] << endl;
    }
}

