#ifndef SUNDANCE_INSTRUCTIONCACHINGEVALUATOR_H
#define SUNDANCE_INSTRUCTIONCACHINGEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceEvaluatorFactory.hpp"
#include "SundanceSparsityPattern.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore {
  namespace Internal {

    /**
     * Factory to produce instruction caching evaluators.
     */
    class InstructionCachingEvaluatorFactory : public EvaluatorFactory
    {
    public:
      /** */
      InstructionCachingEvaluatorFactory();
      
      /** */
      virtual ~InstructionCachingEvaluatorFactory(){;}
      
      /** */
      virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                         int derivSetIndex) const ;
      
    private:
    };

    
    /**
     */
    class InstructionCachingSumEvaluator : public SumEvaluator
    {
    public:
      /** */
      InstructionCachingSumEvaluator(const SumExpr* expr,
                                     int derivSetIndex);

      /** */
      virtual ~InstructionCachingSumEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
      void init();
      Array<int> leftIndex_;
      Array<int> rightIndex_;
      int derivSetIndex_;
      int leftDerivSetIndex_;
      int rightDerivSetIndex_;
      bool isNegative_;
      const EvaluatableExpr* leftExpr_;
      const EvaluatableExpr* rightExpr_;
      SparsityPattern* sparsity_;
      mutable bool needsInit_;
    }; 


    /**
     */
    class InstructionCachingProductEvaluator : public ProductEvaluator
    {
    public:
      /** */
      InstructionCachingProductEvaluator(const ProductExpr* expr,
                                         int derivSetIndex);


      /** */
      virtual ~InstructionCachingProductEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:
      void init();
      Array<Array<int> > leftIndex_;
      Array<Array<int> > rightIndex_;
      int derivSetIndex_;
      int leftDerivSetIndex_;
      int rightDerivSetIndex_;
      const EvaluatableExpr* leftExpr_;
      const EvaluatableExpr* rightExpr_;
      SparsityPattern* sparsity_;
      mutable bool needsInit_;
    };

    

    /**
     */
    class InstructionCachingUnaryMinusEvaluator : public UnaryMinusEvaluator
    {
    public:
      /** */
      InstructionCachingUnaryMinusEvaluator(const UnaryMinus* expr,
                                            int derivSetIndex);

      /** */
      virtual ~InstructionCachingUnaryMinusEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:

      int derivSetIndex_;
      SparsityPattern* sparsity_;
      
    };


    /**
     *
     */
    class InstructionCachingNonlinearUnaryOpEvaluator 
      : public NonlinearUnaryOpEvaluator
    {
    public:
      /** */
      InstructionCachingNonlinearUnaryOpEvaluator(const NonlinearUnaryOp* expr,
                                                  int derivSetIndex);

      /** */
      virtual ~InstructionCachingNonlinearUnaryOpEvaluator(){;}

      /** */
      virtual void eval(const EvalManager& mgr,
                        RefCountPtr<EvalVectorArray>& results) const ;

    private:


      int derivSetIndex_;
      SparsityPattern* sparsity_;
      int maxOrder_;
      int zeroDerivIndex_;
    };
  }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
