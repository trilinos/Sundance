/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALVECTOR_H
#define SUNDANCE_EVALVECTOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceNoncopyable.hpp"

namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    using namespace Teuchos;

    /**
     *
     */
    class EvalVector : public Noncopyable,
      public TSFExtended::ObjectWithVerbosity<EvalVector>
    {
      friend class EvalManager;
      friend class TempStack;

    private:

      /** */
      EvalVector(TempStack* s);

      /** */
      EvalVector(TempStack* s, const RefCountPtr<Array<double> >& data,
                 const string& str);


    public:
      /** 
       * EvalVector has a nontrivial destructor. Upon destruction, 
       * the vector's underlying data object is not destroyed, but rather
       * is put back on the stack of temporary vectors. 
       */
      ~EvalVector();

      /** \name Mathematical operations */
      //@{

      /** */
      void add_SV(const double& alpha, 
                  const EvalVector* B) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this + alpha*B*C
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void add_SVV(const double& alpha,
                   const EvalVector* B,
                   const EvalVector* C) ;

      /** */
      void add_V(const EvalVector* A) ;

      /** */
      void add_S(const double& alpha);

      /**
       * Perform the operation 
       * \f[ 
       * this = this + A*B
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void add_VV(const EvalVector* A,
                  const EvalVector* B) ;


      /**  
       * Perform a scaled addition with another vector,
       * \f[ 
       * this = \alpha this + \beta C
       * \f]
       * The operation is done in-place, overwriting the old values of the
       * vector. 
       */
      void multiply_S_add_SV(const double& alpha, 
                             const double& beta,
                             const EvalVector* C) ;

      /** Scale and add a constant to this vector. 
       * The operation is done in-place, overwriting the old values of
       * the vector. Each element x[i] is updated as:
       * \f[
       * this = alpha * this + beta
       * \f]
       */
      void multiply_S_add_S(const double& alpha,
                            const double& beta) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A + B*C*D
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V_add_VVV(const EvalVector* A,
                              const EvalVector* B,
                              const EvalVector* C,
                              const EvalVector* D) ;


      /**
       * Perform the operation 
       * \f[ 
       * this = this*A + beta*C*D
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V_add_SVV(const EvalVector* A,
                              const double& beta,
                              const EvalVector* C,
                              const EvalVector* D) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A + beta*C
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V_add_SV(const EvalVector* A,
                             const double& beta,
                             const EvalVector* C) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A*B
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_VV(const EvalVector* A,
                       const EvalVector* B) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*alpha*B
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_SV(const double& alpha,
                       const EvalVector* B) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*A
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_V(const EvalVector* A) ;

      /**
       * Perform the operation 
       * \f[ 
       * this = this*alpha
       * \f]
       * which shows up in the chain rule expansion of a second derivative.
       * 
       */
      void multiply_S(const double& alpha) ;

      /**
       *
       */
      void setTo_S_add_SVV(const double& alpha,
                           const double& beta,
                           const EvalVector* C,
                           const EvalVector* D);

      /**
       *
       */
      void setTo_S_add_VV(const double& alpha,
                          const EvalVector* B,
                          const EvalVector* C);

      /**
       *
       */
      void setTo_S_add_SV(const double& alpha,
                           const double& beta,
                           const EvalVector* C);

      /** 
       *
       */
      void setTo_S_add_V(const double& alpha,
                         const EvalVector* B);


      /**
       *
       */
      void setTo_V(const EvalVector* A);

      /**
       *
       */
      void setTo_VV(const EvalVector* A,
                    const EvalVector* B);

      /**
       *
       */
      void setTo_SV(const double& alpha,
                    const EvalVector* B);

      /**
       *
       */
      void setTo_SVV(const double& alpha,
                     const EvalVector* B,
                     const EvalVector* C);

      


      /**
       * Set every element to a constant value
       */
      void setToConstant(const double& alpha) ;

      /** 
       * Apply a unary function
       */
      void applyUnaryOperator(const UnaryFunctor* func, 
                              Array<RefCountPtr<EvalVector> >& opDerivs);
      
      
      /** */
      RefCountPtr<EvalVector> clone() const ;

      /** */
      void resize(int n);

      /** */
      int length() const {return data_->size();}
      
      /** */
      void print(ostream& os) const ;

      /** */
      const double * const start() const {return &((*data_)[0]);}

      /** */
      double * const start() {return &((*data_)[0]);}

      const string& str() const {return str_;}

      void setString(const string& str) {str_ = str;}

      static bool& shadowOps() {static bool rtn = false; return rtn;}

      bool isValid() const {return data_.get() != 0 && s_ != 0;}
      //@}

      

      static double& totalFlops() {static double rtn = 0; return rtn;}

    private:

      static void addFlops(const double& flops) {totalFlops() += flops;}

      mutable TempStack* s_;

      RefCountPtr<Array<double> > data_;

      string str_;

    };
  }
}

namespace std
{
  inline ostream& operator<<(ostream& os, 
                             const SundanceCore::Internal::EvalVector& vec)
  {
    vec.print(os);
    return os;
  }
}

#endif
