#ifndef SUNDANCE_BRUTEFORCEEVALUATOR_H
#define SUNDANCE_BRUTEFORCEEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "Teuchos_TimeMonitor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

Teuchos::Time& sumEvalTimer() ;
Teuchos::Time& productEvalTimer() ;
Teuchos::Time& diffOpEvalTimer() ;
Teuchos::Time& unaryMinusEvalTimer() ;
Teuchos::Time& nonlinearUnaryExprEvalTimer() ;

namespace SundanceCore 
{
  namespace Internal 
  {



    /**
     * Factory to produce brute-force evaluators.
     */
    class BruteForceEvaluatorFactory : public EvaluatorFactory
    {
    public:
      /** */
      BruteForceEvaluatorFactory();
      
      /** */
      virtual ~BruteForceEvaluatorFactory(){;}
      
      /** */
      virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                         int derivSetIndex) const ;
      
    private:
    };

    
    /**
     * BruteForceSumEvaluator evaluates a sum expression with little to
     * no attention to use of sparsity, minimization of temporaries
     * and vector copies, or other such niceties of optimal performance. 
     * It is intended as a bulletproof fallback for testing, and should
     * not be used in production code.
     */
    class BruteForceSumEvaluator : public SumEvaluator
    {
    public:
      /** */
      BruteForceSumEvaluator(const SumExpr* expr)
        : SumEvaluator(expr) {;}

      /** */
      virtual ~BruteForceSumEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
    }; 


    /**
     * BruteForceProductEvaluator evaluates a product expression with
     * little to no attention to use of sparsity, minimization of
     * temporaries and vector copies, or other such niceties of
     * optimal performance. It is intended as a bulletproof fallback
     * for testing, and should not be used in production code.
     */
    class BruteForceProductEvaluator : public ProductEvaluator
    {
    public:
      /** */
      BruteForceProductEvaluator(const ProductExpr* expr)
        : ProductEvaluator(expr) {;}

      /** */
      virtual ~BruteForceProductEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
    };

    /**
     * BruteForceDiffOpEvaluator evaluates a diff op expression with
     * little to no attention to use of sparsity, minimization of
     * temporaries and vector copies, or other such niceties of
     * optimal performance. It is intended as a bulletproof fallback
     * for testing, and should not be used in production code.
     */
    class BruteForceDiffOpEvaluator : public DiffOpEvaluator
    {
    public:
      /** */
      BruteForceDiffOpEvaluator(const DiffOp* expr)
        : DiffOpEvaluator(expr) {;}

      /** */
      virtual ~BruteForceDiffOpEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
    };

    /**
     * BruteForceUnaryMinusEvaluator evaluates a unary minus
     * expression with little to no attention to use of sparsity,
     * minimization of temporaries and vector copies, or other such
     * niceties of optimal performance. It is intended as a
     * bulletproof fallback for testing, and should not be used in
     * production code.
     */
    class BruteForceUnaryMinusEvaluator : public UnaryMinusEvaluator
    {
    public:
      /** */
      BruteForceUnaryMinusEvaluator(const UnaryMinus* expr)
        : UnaryMinusEvaluator(expr) {;}

      /** */
      virtual ~BruteForceUnaryMinusEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
    };


    /**
     *
     */
    class BruteForceNonlinearUnaryOpEvaluator 
      : public NonlinearUnaryOpEvaluator
    {
    public:
      /** */
      BruteForceNonlinearUnaryOpEvaluator(const NonlinearUnaryOp* expr)
        : NonlinearUnaryOpEvaluator(expr) {;}

      /** */
      virtual ~BruteForceNonlinearUnaryOpEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
    }; 


    /**
     *
     */
    class BruteForceUserDefOpEvaluator 
      : public UserDefOpEvaluator
    {
    public:
      /** */
      BruteForceUserDefOpEvaluator(const UserDefOp* expr)
        : UserDefOpEvaluator(expr) {;}

      /** */
      virtual ~BruteForceUserDefOpEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
