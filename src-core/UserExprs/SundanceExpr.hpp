/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXPR_H
#define SUNDANCE_EXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExprBase.hpp"
#include "TSFHandle.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace SundanceCore
{
  using namespace SundanceUtils;
  /**
   * User-level expression class. Expr is a handle to a
   * reference-counted pointer to a ExprBase subtype. As such,
   * expression copies and assignments are shallow.
   */
  class Expr
    {
    public:
      /** empty ctor */
      Expr() : ptr_() {;}

      /** Construct with a pointer to an ExprBase subtype. */
      Expr(Internal::ExprBase* ptr);

      /** Construct with a constant. Creates a ConstantExpr. */
      Expr(const double& c);

      /** Add two expressions. The operands must have identical list
       * structures.
       */
      Expr operator+(const Expr& other) const ;
      /** Subtract two expressions. The operands must have identical
       *  list structures. */
      Expr operator-(const Expr& other) const ;
      /** Multiply two expressions. The operands must have list
       * structures such that the product can be interpreted as a
       * scalar-vector product or as an inner product between vectors
       * or tensors. The multiplication operator is also used to
       * represent the application of a differential operator.
       */
      Expr operator*(const Expr& other) const ;
      /** Divide one expression by another. The right operand must be
          a scalar. */
      Expr operator/(const Expr& other) const ;

      /** Unary minus operator */
      Expr operator-() const ;

      /** List element accessor */
      const Expr& operator[](int i) const ;
      
      /** Number of elements in top level of list */
      int size() const ;

      /** Total number of elements in list */
      int totalSize() const ;

      /** Append a new element to this list */
      void append(const Expr& expr);

      /** Flatten this list */
      Expr flatten() const ;

      /** */
      string toLatex() const ;

      /** */
      string toString() const ;

      /** */
      XMLObject toXML() const ;


      /**
       *
       */
      void showDerivs() const;


#ifndef DOXYGEN_DEVELOPER_ONLY

      /**
       * Expose the pointer underlying this object.
       */
      const RefCountPtr<Internal::ExprBase>& ptr() const {return ptr_;}

      /**
       * Turn evaluation caching on
       */
      static bool& evaluationCachingOn() {static bool rtn = true; return rtn;}

      /**
       * Show parentheses around every pair of operands
       */
      static bool& showAllParens() {static bool rtn = false; return rtn;}

      /** Create a new handle for an existing ptr.  */
      static Expr handle(const RefCountPtr<Internal::ExprBase>& ptr)
      {return Expr(ptr);}

      /** */
      static Time& opTimer() 
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Expr symbolic ops"); 
        return *rtn;
      }
      /** */
      static Time& outputTimer() 
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Expr output"); 
        return *rtn;
      }


    private:
      RefCountPtr<Internal::ExprBase> ptr_;

      /** Construct an expression handle given an existing RCP */
      Expr(const RefCountPtr<Internal::ExprBase>& ptr);

      /** Add two scalar expressions */
      Expr sum(const Expr& other, int sign) const ;

      /** Multiply two scalar expressions */
      Expr multiply(const Expr& other) const ;

      /** Divide two scalar expressions */
      Expr divide(const Expr& other) const ;

#endif /* DOXYGEN_DEVELOPER_ONLY */
    };

  /** \relates Expr */
  inline ostream& operator<<(ostream& os, const Expr& e)
    {
      return e.ptr()->toText(os, false);
    }

  /** \relates Expr */
  inline Expr operator+(const double& a, const Expr& x)
    {return Expr(a) + x;}

  /** \relates Expr */
  inline Expr operator-(const double& a, const Expr& x)
    {return Expr(a) - x;}

  /** \relates Expr */
  inline Expr operator*(const double& a, const Expr& x)
    {return Expr(a) * x;}

  /** \relates Expr */
  inline Expr operator/(const double& a, const Expr& x)
    {return Expr(a) / x;}

  /** \relates Expr */
  Expr List(const Expr& a);

  /** \relates Expr */
  Expr List(const Expr& a, const Expr& b);

  /** \relates Expr */
  Expr List(const Expr& a, const Expr& b, const Expr& c);

  /** \relates Expr */
  Expr List(const Expr& a, const Expr& b, const Expr& c,
            const Expr& d);

  
}

#endif
