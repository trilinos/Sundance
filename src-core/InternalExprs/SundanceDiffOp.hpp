/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_DIFFOP_H
#define SUNDANCE_DIFFOP_H

#include "SundanceDefs.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  using namespace Teuchos;

  using std::string;
  using std::ostream;

  namespace Internal
  {
    /**
     *
     */
    class DiffOp : public UnaryExpr
    {
    public:
      /** ctor */
      DiffOp(const MultiIndex& op, const RefCountPtr<ScalarExpr>& arg);

      /** virtual destructor */
      virtual ~DiffOp() {;}

      /** Write a simple text description suitable
       * for output to a terminal */
      virtual ostream& toText(ostream& os, bool paren) const ;

      /** Write in a form suitable for LaTeX formatting */
      virtual ostream& toLatex(ostream& os, bool paren) const ;

      /** Write in XML */
      virtual XMLObject toXML() const ;

      /**
       * Indicate whether the given functional derivative is nonzero
       */
      virtual bool hasNonzeroDeriv(const MultipleDeriv& f) const ;

      /**
       * Find all functions and their derivatives beneath my level
       * in the tree. For FuncElementBase, this appends a derivative
       * wrt me.
       */
      virtual void getRoughDependencies(Set<Deriv>& funcs) const ;

      /** */
      const Deriv& myCoordDeriv() const {return myCoordDeriv_;}

      /** */
      const MultiIndex& mi() const {return mi_;}

      /** Get the functions that are required in the evaluation
       * of the multiple deriv d */
      const SundanceUtils::Set<Deriv>& requiredFunctions(const MultipleDeriv& d) const 
      {return requiredFunctions_[d];}

      /** */
      bool requiresFunctionsToEval(const MultipleDeriv& d) const 
      {return requiredFunctions_.containsKey(d);}

    

      /** Return the set of derivatives required by the operands
       * of this expression given that this expression
       * requires the set d. For all expressions other than DiffOp,
       * the operand derivative set is identical to the input derivative
       * set. DiffOp will require a different set of derivatives from
       * its operand. */
      virtual Array<DerivSet>
      derivsRequiredFromOperands(const DerivSet& d) const ;

      /** */
      virtual RefCountPtr<ExprBase> getRcp() {return rcp(this);}



    private:

      MultiIndex mi_;


      Deriv myCoordDeriv_;

      mutable Map<MultipleDeriv, SundanceUtils::Set<Deriv>, 
                  increasingOrder<MultipleDeriv> > requiredFunctions_;

    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
