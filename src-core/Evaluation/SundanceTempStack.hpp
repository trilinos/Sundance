/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TEMPSTACK_H
#define SUNDANCE_TEMPSTACK_H

#include "SundanceDefs.hpp"
#include "SundanceEvalVectorArray.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceSparsityPattern.hpp"
#include <stack>

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
    {



      /**
       * TempStack provides a stack of temporary variables for use during
       * evaluation. 
       *
       * During the course of evaluating an expression, it is often necessary
       * to create temporary variables. For example, in evaluating
       * \code
       * a += b*(c+d)
       * \endcode
       * it is required to create three temporaries (ignoring for
       * explanatory purposes any copies made in cases where one of
       * the vectors will be used elsewhere). We can see this
       * by breaking the operation
       * down into the following steps:
       * \code
       * 1. Create a temporary variable t1
       * 2. Evaluate expression b into t1
       * 3. Create a temporary variable t2
       * 4. Evaluate expression c into t2
       * 3. Create a temporary variable t3
       * 4. Evaluate expression d into t3
       * 5. Carry out t2 += t3
       * 6. Carry out t1 *= t2
       * 7. Carry out a += t1
       * \endcode
       * The number of temporaries required for a given expression
       * will depend on the graph of the expression. In general, we want to
       * create exactly as many temporaries as are needed, and reuse any
       * temporaries that are no longer needed. This is a well-known problem
       * in compiler design, and can be accomplished by maintaining a
       * stack of temporaries. When a new temporary is needed, it is popped
       * from the stack; if the stack is empty, a new temporary is allocated.
       * When a step of a calculation is done, any temporaries used are
       * put back on the stack for further use.
       *
       * In our case, however, there are several additional wrinkles:
       * <ul>
       * <li> Our vectors can be "trivial", i.e., zero or constant with
       * no need to allocate a vector representation,
       * or they can be "full", i.e., have a full vector allocated. We
       * don't want to allocate a vector unless necessary. To deal with this
       * issue, we maintain two stacks, one for trivial vectors and one
       * for full vectors, with two methods, popTrivialVector() and
       * popFullVector(), to pop a vector from the appropriate stack.
       * Since we will always know from our sparsity data which kind of vector
       * we need, we will know which method to call. There is a single
       * method, pushVector(), for pushing vectors back onto the stack; 
       * it takes responsibility for determining which stack the 
       * vector should be pushed onto.
       * <li> We often allocate arrays of vectors at a time. Thus we have
       * popVectorArray() and pushVectorArray() to pop and push arrays
       * of vectors. The popVectorArray() method takes a SparsityPattern
       * as input because it needs to know whether each of its elements
       * should be trivial or full.
       * </ul>
       */
      class TempStack : public Noncopyable
        {
        public:
          /** Empty ctor. This stack will supply vectors for string
           * evaluation */
          TempStack();

          /** Construct with an initial vector size */
          TempStack(int vecSize);

          /** Pop a trivial (constant-valued) vector from the stack. */
          RefCountPtr<EvalVector> popTrivialVector() ;

          /** Pop a full (non-constant-valued) vector from the stack. */
          RefCountPtr<EvalVector> popFullVector() ;

          /** Push a vector onto the stack */
          void pushVector(const RefCountPtr<EvalVector>& vec) ;

          /** Pop a vector array from the stack, with each vector
           * allocated as trivial or full as
           * indicated by the sparsity pattern */
          RefCountPtr<EvalVectorArray> 
          popVectorArray(const SparsityPattern* sparsity) ;

          /** Push a vector array onto the stack */
          void pushVectorArray(const RefCountPtr<EvalVectorArray>& vecs) ;

          /** */
          void setVecSize(int vecSize) {vecSize_ = vecSize;}

          /** */
          void resetCounters() ;

          /** */
          int numTrivialVecsAccessed() const {return trivialVecsAccessed_;}

          /** */
          int numFullVecsAccessed() const {return fullVecsAccessed_;}

          /** */
          int numVecArraysAccessed() const {return vecArraysAccessed_;}

          /** */
          int numTrivialVecsAllocated() const {return trivialVecsAllocated_;}

          /** */
          int numFullVecsAllocated() const {return fullVecsAllocated_;}

          /** */
          int numVecArraysAllocated() const {return vecArraysAllocated_;}

        private:
          
          bool numericalEval_;

          int vecSize_;

          std::stack<RefCountPtr<EvalVector> > fullVecStack_;

          std::stack<RefCountPtr<EvalVector> > trivialVecStack_;

          std::stack<RefCountPtr<EvalVectorArray> > vecArrayStack_;

          int trivialVecsAllocated_;

          int fullVecsAllocated_;

          int vecArraysAllocated_;

          int trivialVecsAccessed_;

          int fullVecsAccessed_;

          int vecArraysAccessed_;
        };
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
