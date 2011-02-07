/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTOROPSIMPL_HPP
#define PLAYA_VECTOROPSIMPL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaVectorOpsDecl.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_RCP.hpp"

namespace PlayaFunctors
{
template <class Scalar> class BoundedMinLocFunctor;

template <class Scalar> class BoundedMaxLocFunctor;
}

namespace Playa
{

using namespace PlayaFunctors;

/* */
template <class Scalar> inline
Scalar minloc(const Vector<Scalar>& x, int& gni)
{
  return minlocWithBound(-HUGE_VAL, x, gni);
}

/* */
template <class Scalar> inline
Scalar maxloc(const Vector<Scalar>& x, int& gni)
{
  return maxlocWithBound(-HUGE_VAL, x, gni);
}

/* */
template <class Scalar> inline
Scalar minlocWithBound(const Scalar& lowerBound, 
  const Vector<Scalar>& x, int& gni)
{
  IndexedValue<Scalar> y = 
    x.applyUnaryReductionFunctor(BoundedMinLocFunctor<Scalar>(x.comm(), lowerBound, x.space().baseGlobalNaturalIndex()));
  gni = y.where;
  return y.what;
}

/* */
template <class Scalar> inline
Scalar maxlocWithBound(const Scalar& upperBound, 
  const Vector<Scalar>& x, int& gni)
{
  IndexedValue<Scalar> y = 
    x.applyUnaryReductionFunctor(BoundedMaxLocFunctor<Scalar>(x.comm(), upperBound, x.space().baseGlobalNaturalIndex()));
  gni = y.where;
  return y.what;
}

} // end namespace Playa

namespace PlayaFunctors
{


using namespace Playa;

/** \brief Find minimum element above a lower bound, returning value and
 * location */
template <class Scalar>
class BoundedMinLocFunctor : public ReductionFunctorBase<Scalar>
{
public:
  /** */
  BoundedMinLocFunctor(const MPIComm& comm, const Scalar& bound,
    int baseGNI)
    : ReductionFunctorBase<Scalar>(comm), min_(), 
      bound_(bound), baseGNI_(baseGNI)
    {
      min_.what = HUGE_VAL;
      min_.where = -1;
    }

  /** */
  void step(int i, const Scalar& x) const 
    {
      if (x < min_.what && x > bound_) 
      {
        min_.what = x;
        min_.where = i;
      }
    }

  /** */
  void postProc() const 
    {
      min_.where += baseGNI_;

      IndexedValue<Scalar> out = min_;

      this->comm().allReduce(&min_, &out, 1, MPIComm::DOUBLE_INT, 
        MPIComm::MINLOC);
      min_ = out;
    }

  /** */
  IndexedValue<Scalar> result() const 
    {
      return min_;
    }

private:
  MPIComm comm_;
  mutable IndexedValue<Scalar> min_;
  Scalar bound_;
  int baseGNI_;
};


/** \brief Find maximum element below an upper bound, returning value and
 * location */
template <class Scalar>
class BoundedMaxLocFunctor : public ReductionFunctorBase<Scalar>
{
public:
  /** */
  BoundedMaxLocFunctor(const MPIComm& comm, const Scalar& bound,
    int baseGNI)
    : ReductionFunctorBase<Scalar>(comm), max_(), 
      bound_(bound), baseGNI_(baseGNI)
    {
      max_.what = -HUGE_VAL;
      max_.where = -1;
    }

  /** */
  void step(int i, const Scalar& x) const 
    {
      if (x > max_.what && x < bound_) 
      {
        max_.what = x;
        max_.where = i;
      }
    }

  /** */
  void postProc() const 
    {
      max_.where += baseGNI_;

      IndexedValue<Scalar> out = max_;

      this->comm().allReduce(&max_, &out, 1, MPIComm::DOUBLE_INT, 
        MPIComm::MAXLOC);
      max_ = out;
    }

  /** */
  IndexedValue<Scalar> result() const 
    {
      return max_;
    }

private:
  MPIComm comm_;
  mutable IndexedValue<Scalar> max_;
  Scalar bound_;
  int baseGNI_;
};



/** \brief Specify return type of BoundedMinLocFunctor */
template <class Scalar>
class VectorFunctorTraits<Scalar, BoundedMinLocFunctor<Scalar> >
{
public:
  typedef IndexedValue<Scalar> ReturnType ;
};


/** \brief Specify return type of BoundedMaxLocFunctor */
template <class Scalar>
class VectorFunctorTraits<Scalar, BoundedMaxLocFunctor<Scalar> >
{
public:
  typedef IndexedValue<Scalar> ReturnType ;
};

  
}

#endif
