/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FUNCTIONALEVALUATOR_H
#define SUNDANCE_FUNCTIONALEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceAssembler.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
  using namespace SundanceUtils;
  using namespace SundanceStdMesh;
  using namespace SundanceStdMesh::Internal;
  using namespace SundanceCore;
  using namespace SundanceCore::Internal;

  /** */
  double evaluateIntegral(const Mesh& mesh, const Expr& expr);

  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * 
     */
    class FunctionalEvaluator 
      : public TSFExtended::ObjectWithVerbosity<FunctionalEvaluator>
    {
    public:
      /** */
      FunctionalEvaluator(const Mesh& mesh, 
                          const Expr& integral);
      /** */
      FunctionalEvaluator(const Mesh& mesh, 
                          const Expr& integral,
                          const Expr& bcs,
                          const Expr& var,
                          const Expr& varEvalPts,
                          const VectorType<double>& vectorType);
      /** */
      FunctionalEvaluator(const Mesh& mesh, 
                          const Expr& integral,
                          const Expr& bcs,
                          const Expr& vars,
                          const Expr& varEvalPts,
                          const Expr& fields,
                          const Expr& fieldValues,
                          const VectorType<double>& vectorType);


      /** */
      double evaluate() const ;

      /** */
      Expr evalGradient(double& value) const ;

      /** */
      double fdGradientCheck(double h) const ;
      
    private:

      /** */
      Vector<double> evalGradientVector(double& value) const ;
      
      /** */
      RefCountPtr<Assembler> assembler_;
      
      /** */
      mutable Expr varValues_;

      /** */
      VectorType<double> vecType_;
      
      /** */
      mutable Vector<double> gradient_;
      
    };
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
