/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EVALVECTOR_H
#define SUNDANCE_EVALVECTOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "TSFObjectWithVerbosity.hpp"
#include "SundanceUnaryFunctor.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceCore
{
  using namespace SundanceUtils;
  namespace Internal
  {
    using namespace Teuchos;

    /**
     * EvalVector is used internally by Sundance for numerical evaluation
     * of expressions and their functional derivatives.
     *
     * Why yet another vector class? EvalVector has two significant
     * features not supported by existing vector classes: (1) since we
     * will often be able to identify symbolically terms that are
     * zero, one, or constant in certain regions, EvalVector has the
     * ability to represent those values exactly and save vector
     * calculations; (2) EvalVector has the ability to "shadow" a numerical
     * calculation with a string calculation which produces a string
     * showing all the steps carried out during the calculation. 
     *
     * The methods isZero(), isOne(), and isConstant() indicate whether
     * the vector has a constant or trivial value that can be used to simplify
     * calculations and avoid full vector operations. 
     *
     * Whatever framework Sundance is connected to 
     * will need to insert numerical values into an
     * EvalVector. This is done using the LoadableVector interface,
     * implemented by EvalVector.
     */
    class EvalVector : public TSFExtended::ObjectWithVerbosity<EvalVector>
    {
      friend class TempStack;

    private:
      /** Create a vector, initialized to zero. Ctors are private; the
       * friend class TempStack can create EvalVectors */
      EvalVector();
      
      /** Create a vector of length n. Ctors are private; the
       * friend class TempStack can create EvalVectors*/
      EvalVector(int n);

    public:

      /** \name Element loading functions */
      //@{
      /** Change the size of the vector to newSize */
      void resize(int newSize) {vectorVal_.resize(newSize);}
          
      /** Return the length of the vector */
      int length() const {return vectorVal_.length();}

      /** Set the i-th element to x */
      void setElement(int i, const double& x) {vectorVal_[i] = x;}

      /** Return a pointer to the physical start of the vector. */
      double* const start() {return &(vectorVal_[0]);}

      /** Return a pointer to the physical start of the vector. */
      const double* const start() const {return &(vectorVal_[0]);}
      //@}

      /** Return true if we are shadowing numerical operations with
       * string operations */
      static bool& shadowOps() {static bool rtn = false; return rtn;}

      /** \name Methods for handling special values */
      //@{
      /** Return true if this vector is a constant */
      bool isConstant() const {return isConstant_;}
      
      /** Return true if this vector is zero */
      bool isZero() const {return isZero_;}
      
      /** Return true if this vector is one */
      bool isOne() const {return isOne_;}

      /** Set this vector to zero */
      void setToZero();

      /** Set this vector to one */
      void setToOne();

      /** Set this vector to a constant value */
      void setToConstantValue(const double& constantVal);

      /** Return the constant value of this vector */
      const double& getConstantValue() const {return constantVal_;}

      /** Make this vector use its vector values rather than
       * any constant values. The elements of the vector are assumed
       * to be set elsewhere; this method simply sets isConstant,
       * isZero, and isOne to false. */
      void setToVectorValue();
      //@}

      /** \name Methods for working with string values */
      //@{
      /** Return the string value of this vector */
      string getStringValue() const ;

      /** Set the string value of this vector, concurrently setting
       * isConstant, isZero, and isOne to false. */
      void setStringValue(const string& stringVal);
      //@}

      /** \name Mathematical operations */
      //@{
      /** Add another vector to this one. Special cases in which
       * either operand is zero or constant are checked. If numerical()
       * is true, the calculation is done numerically; otherwise,
       * a string computation is performed. */
      void addScaled(const RefCountPtr<EvalVector>& other,
                     const double& scalar) ;

      /** Elementwise multiply this vector with another vector.
       * Special cases in which either operand is zero, one, or
       * constant are checked. If numerical() is true, the
       * calculation is done numerically; otherwise, a string
       * computation is performed. */
      void multiply(const RefCountPtr<EvalVector>& other) ;

      /** Add a product a*b to this vector. Special cases in which
       * either operand is zero, one, or constant are checked. If
       * numerical() is true, the calculation is done numerically;
       * otherwise, a string computation is performed.*/
      void addProduct(const RefCountPtr<EvalVector>& a,
                      const RefCountPtr<EvalVector>& b) ;

      /** Take the square root of this vector. If numerical() is
       * true, the calculation is done numerically; otherwise, a
       * string computation is performed. */
      void sqrt() ;

      /** Apply the given unary function to this vector, returning
       * the value and the requested number of derivatives */
      void applyUnaryFunction(const UnaryFunctor* func,
                              Array<RefCountPtr<EvalVector> >& funcDerivs) const ;

      /** Copy a vector into this vector. If numerical() is true,
       * the calculation is done numerically; otherwise, a string
       * computation is performed.*/
      void copy(const RefCountPtr<EvalVector>& other) ;

      /** negate this vector.  */
      void unaryMinus() ;
      //@}

      /** */
      void print(ostream& os) const ;

    private:
      /** The elements of this vector, used if isConstant() is false and
       * we are doing numerical calculations. */
      Array<double> vectorVal_;

      /** the string value of this vector, used if isConstant()
       * is false and we are doing string, rather than numerical,
       * calculations */
      string stringVal_;

      /** the constant value of this vector, used if isConstant() is true */
      double constantVal_;

      /** flag indicating whether all elements of this vector are
       * equal to a real constant */
      bool isConstant_;

      /** flag indicating whether all elements of this vector are
       * equal to zero (0.0) */
      bool isZero_;

      /** flag indicating whether all elements of this vector are
       * equal to one (1.0) */
      bool isOne_;

    };
  }
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
