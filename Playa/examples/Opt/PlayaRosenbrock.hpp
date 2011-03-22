/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_ROSENBROCK_HPP
#define PLAYA_ROSENBROCK_HPP

#include "PlayaObjectiveBase.hpp"
#include "PlayaVectorType.hpp"

namespace Playa
{
using std::ostringstream;

/**
 * Class Rosenbrock is an example of a user implementation of an objective
 * function by derivation from the ObjectiveBase abstract interface. 
 *
 * The extended Rosenbrock function is 
 * \f[
 * f({\bf x}) = \sum_{m=0}^{M-1} \alpha \left(x_{2m+1}-x_{2m}^2\right)^2
 * + \left(x_{2m}-1\right)^2 
 * \f]
 * This function has a curved valley with a unique 
 * minimum at \f$[1,1,\cdots,1,1]^T\f$. 
 * The parameter \f$\alpha\f$ tunes the shape of the valley, with larger 
 * values of \f$\alpha\f$ presenting greater difficulty for numerical
 * optimization. 
 *
 * This class must implement the pure virtual functions from the
 * ObjectiveBase class. These are
 * <ul>
 * <li> eval()
 * <li> evalGrad()
 * <li> getInit()
 * </ul> 
 * Additionally, it is necessary to write a constructor that 
 * creates a Rosenbrock function with the desired values of \f$\alpha\f$
 * and \f$M\f$.
 */
class Rosenbrock : public ObjectiveBase
{
public:
  /** 
   * \brief Constructor for the Rosenbrock function object. A user can use this
   * constructor to make a Rosenbrock function with any choice of parameters
   * \f$\alpha\f$ and \f$M\f$ and with a user-specified low-level
   * vector representation, indicated by the VectorType input argument.
   * 
   */
  Rosenbrock(int N, double alpha, const VectorType<double>& vecType)
    : N_(N), vs_(vecType.createEvenlyPartitionedSpace(MPIComm::self(), 2*N)),
      alpha_(alpha)
    {}

  /** \name Implementations of the ObjectiveBase pure virtual functions */
  //@{

  /**
   * \brief Evaluate the function and its gradient
   *
   * \param x [in] the point at which the function is to be evaluated.
   * \param f [out] the computed value of \f$f(x)\f$
   * \param grad [out] the computed value of \f$\nabla f(x)\f$
   */
  void  evalGrad(const Vector<double>& x, double& f, 
    Vector<double>& grad) const ;

  /** 
   * \brief Evaluate the function 
   * 
   * \param x [in] the point at which the function is to be evalulated
   * \param f [out] the computed value of \f$f(x)\f$
   */
  void eval(const Vector<double>& x, double& f) const ;

  /** 
   * \brief Return a vector to be used as the initial guess for the 
   * optimization loop.
   */
  Vector<double> getInit() const ;
  //@}

  /** \brief Return a short description of this object */
  string description() const 
    {
      ostringstream oss;
      oss << "Rosenbrock[n=" << N_ << ", alpha=" << alpha_ << "]";
      return oss.str();
    }

private:
  int N_;
  VectorSpace<double> vs_;
  double alpha_;
};


}

#endif
