/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLDIAMETEREXPR_H
#define SUNDANCE_CELLDIAMETEREXPR_H

#include "SundanceLeafExpr.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCellDiameterEvaluator.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace SundanceCore::Internal;

  /** */
  class CellDiameterExpr 
    : public GenericEvaluatorFactory<CellDiameterExpr, CellDiameterExprEvaluator>,
      virtual public LeafExpr
  {
  public:
    /** */
    CellDiameterExpr(const string& name="h");
    
    /** */
    virtual ~CellDiameterExpr() {;}

    /** */
    virtual XMLObject toXML() const ;

    const string& name() const {return name_;}

    /** Write a simple text description suitable 
     * for output to a terminal */
    virtual ostream& toText(ostream& os, bool paren) const ;
    
    /** Write in a form suitable for LaTeX formatting */
    virtual ostream& toLatex(ostream& os, bool paren) const ;

    /*
     * Determine which functional and spatial derivatives are nonzero in the
     * given context. We also keep track of which functional derivatives
     * are known to be constant, which can simplify evaluation. 
     */
    virtual void findNonzeros(const EvalContext& context,
                              const Set<MultiIndex>& multiIndices,
                              bool regardFuncsAsConstant) const ;

    /** */
    virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}
  private:
    string name_;
  };
}

#endif
