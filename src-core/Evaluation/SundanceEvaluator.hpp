/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALUATOR_H
#define SUNDANCE_EVALUATOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceEvalVectorArray.hpp"
#include "TSFObjectWithVerbosity.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore 
{
  class CoordExpr;
  class UserDefOp;

  namespace Internal
  {
    class EvalContext;
  }
  
  using namespace Internal;

  namespace Internal 
  {
    class EvalManager;
    class SumExpr;
    class ProductExpr;
    class DiffOp;
    class UnaryMinus;
    class SpatiallyConstantExpr;
    class SymbolicFuncElement;
    class DiscreteFuncElement;
    class NonlinearUnaryOp;


    /**
     * Base class for evaluator objects. Sundance uses pluggable objects
     * to evaluate expressions, so that higher-performance evaluators
     * can be developed and introduced with minimal disruption to the
     * core expression code. 
     */
    class Evaluator : public TSFExtended::ObjectWithVerbosity<Evaluator>
    {
    public:
      /** */
      Evaluator();

      /** */
      virtual ~Evaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const = 0 ;


    protected:
    private:
    };

    /**
     * 
     */
    template <class ExprType> class SubtypeEvaluator : public Evaluator
    {
    public:
      /** */
      SubtypeEvaluator(const ExprType* expr)
        : Evaluator(), expr_(expr) {;}

      /** */
      virtual ~SubtypeEvaluator(){;}


    protected:
      /** */
      const ExprType* expr() const {return expr_;}

    private:
      const ExprType* expr_;
    };

    typedef SubtypeEvaluator<SumExpr> SumEvaluator;

    typedef SubtypeEvaluator<ProductExpr> ProductEvaluator;

    typedef SubtypeEvaluator<DiffOp> DiffOpEvaluator;

    typedef SubtypeEvaluator<UnaryMinus> UnaryMinusEvaluator;

    typedef SubtypeEvaluator<NonlinearUnaryOp> NonlinearUnaryOpEvaluator;

    typedef SubtypeEvaluator<UserDefOp> UserDefOpEvaluator;

    /** 
     *
     */
    class CoordExprEvaluator : public SubtypeEvaluator<CoordExpr>
    {
    public:
      /** */
      CoordExprEvaluator(const CoordExpr* expr);

      /** */
      virtual ~CoordExprEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;
      
    };

    /** 
     *
     */
    class ConstantEvaluator : public SubtypeEvaluator<SpatiallyConstantExpr>
    {
    public:
      /** */
      ConstantEvaluator(const SpatiallyConstantExpr* expr);

      /** */
      virtual ~ConstantEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;
      
    }; 

    /** 
     *
     */
    class SymbolicFuncElementEvaluator 
      : public SubtypeEvaluator<SymbolicFuncElement>
    {
    public:
      /** */
      SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr);

      /** */
      virtual ~SymbolicFuncElementEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;
      
    };

    /** 
     *
     */
    class DiscreteFuncElementEvaluator 
      : public SubtypeEvaluator<DiscreteFuncElement>
    {
    public:
      /** */
      DiscreteFuncElementEvaluator(const DiscreteFuncElement* expr);

      /** */
      virtual ~DiscreteFuncElementEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;
      
    };

  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
