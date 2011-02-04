/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_VECTORFUNCTORSIMPL_HPP
#define PLAYA_VECTORFUNCTORSIMPL_HPP


#include "PlayaDefs.hpp"
#include "PlayaVectorFunctorsDecl.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaRand.hpp"

namespace Playa
{

/**
 * IndexedValue is the return type for reduction operations such
 * as MinLoc that return a location and a value. 
 */
template <class Scalar>
struct IndexedValue
{
  /** Value */
  Scalar what;
  /** Index */
  GeneralizedIndex where;
};

/** Functor for elementwise absolute value */
template <class Scalar>
class AbsFunctor
{
public:
  /** */
  AbsFunctor() {}

  /** */
  Scalar operator()(const Scalar& x) const 
    {return abs(x);}
};

/** Functor for elementwise reciprocal */
template <class Scalar>
class ReciprocalFunctor
{
public:
  /** */
  ReciprocalFunctor() {}

  /** */
  Scalar operator()(const Scalar& x) const 
    {return 1.0/x;}
};

/** Functor for elementwise reciprocal */
template <class Scalar>
class RandomFunctor
{
public:
  /** */
  RandomFunctor() {}

  /** */
  Scalar operator()(const Scalar& x) const 
    {return Rand::val();}
};

/** Functor for multiplication by a scalar */
template <class Scalar>
class ScalarMultFunctor
{
public:
  /** */
  ScalarMultFunctor(const Scalar& alpha) : alpha_(alpha) {}

  /** */
  Scalar operator()(const Scalar& x) const 
    {return alpha_*x;}
private:
  Scalar alpha_;
};

/** Identity functor, used for copying */
template <class Scalar>
class IdentityFunctor
{
public:
  /** */
  IdentityFunctor() {}

  /** */
  Scalar operator()(const Scalar& x) const 
    {return x;}
};

/** Functor for setting all elements to a constant */
template <class Scalar>
class ConstantFunctor
{
public:
  /** */
  ConstantFunctor(const Scalar& alpha) : alpha_(alpha) {}

  /** */
  Scalar operator()(const Scalar& x) const 
    {return alpha_;}
private:
  Scalar alpha_;
};


/** Functor for elementwise product */
template <class Scalar>
class DotStarFunctor
{
public:
  /** */
  DotStarFunctor() {}

  /** */
  Scalar operator()(const Scalar& x, const Scalar& y) const 
    {return x*y;}
};


/** Functor for elementwise quotient */
template <class Scalar>
class DotSlashFunctor
{
public:
  /** */
  DotSlashFunctor() {}

  /** */
  Scalar operator()(const Scalar& x, const Scalar& y) const 
    {return x*y;}
};

/** Functor for linear combination of two vectors */
template <class Scalar>
class LCFunctor2
{
public:
  /** */
  LCFunctor2(const Scalar& a, const Scalar& b) : a_(a), b_(b) {}

  /** */
  Scalar operator()(const Scalar& x, const Scalar& y) const 
    {return a_*x + b_*y;}

private:
  Scalar a_;
  Scalar b_;
};


/** Functor for linear combination of two vectors */
template <class Scalar>
class LCFunctor3
{
public:
  /** */
  LCFunctor3(const Scalar& a, const Scalar& b, const Scalar& c)
    : a_(a), b_(b), c_(c) {}

  /** */
  Scalar operator()(const Scalar& x, const Scalar& y, const Scalar& z) const 
    {return a_*x + b_*y + c_*z;}

private:
  Scalar a_;
  Scalar b_;
  Scalar c_;
};


/** Functor for Euclidean norm of a vector */
template <class Scalar>
class Norm2Functor : public ReductionFunctorBase<Scalar>
{
public:
  Norm2Functor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x) const 
    {
      val_ += x*x;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::SUM);
      val_ = final;
    }

  Scalar result() const 
    {
      return ::sqrt(val_);
    }

private:
  mutable Scalar val_;
};

/** Functor for weighted 2-norm of a vector */
template <class Scalar>
class WeightedNorm2Functor : public ReductionFunctorBase<Scalar>
{
public:
  WeightedNorm2Functor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x, const Scalar& y) const 
    {
      val_ += y*x*x;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::SUM);
      val_ = final;
    }

  Scalar result() const 
    {
      return ::sqrt(val_);
    }

private:
  MPIComm comm_;
  mutable Scalar val_;
};

/** Functor for 1-norm */
template <class Scalar>
class Norm1Functor : public ReductionFunctorBase<Scalar>
{
public:
  Norm1Functor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x) const 
    {
      val_ += ::fabs(x);
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::SUM);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

private:
  mutable Scalar val_;
};

/** Functor for infinity norm of a vector */
template <class Scalar>
class NormInfFunctor : public ReductionFunctorBase<Scalar>
{
public:
  NormInfFunctor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(-1.0) {}

  void step(int i, const Scalar& x) const 
    {
      Scalar z = ::fabs(x);
      if (z > val_) val_ = z;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::MAX);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

private:
  mutable Scalar val_;
};

/** Functor for dot product of two vectors */
template <class Scalar>
class DotProductFunctor : public ReductionFunctorBase<Scalar>
{
public:
  DotProductFunctor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(0.0) {}

  void step(int i, const Scalar& x, const Scalar& y) const 
    {
      val_ += x*y;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::SUM);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

private:
  mutable Scalar val_;
};


/** Functor to compute minimum element of a vector */
template <class Scalar>
class MinFunctor : public ReductionFunctorBase<Scalar>
{
public:
  MinFunctor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(HUGE_VAL) {}

  void step(int i, const Scalar& x) const 
    {
      if (x < val_) val_ = x;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::MIN);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

private:
  mutable Scalar val_;
};

/** Functor to compute minimum element of a vector */
template <class Scalar>
class MaxFunctor : public ReductionFunctorBase<Scalar>
{
public:
  MaxFunctor(const MPIComm& comm)
    : ReductionFunctorBase<Scalar>(comm), val_(-HUGE_VAL) {}

  void step(int i, const Scalar& x) const 
    {
      if (x > val_) val_ = x;
    }

  void postProc() const 
    {
      Scalar final = val_;
      this->comm().allReduce(&val_, &final, 1, MPIComm::DOUBLE, MPIComm::MAX);
      val_ = final;
    }

  Scalar result() const 
    {
      return val_;
    }

private:
  mutable Scalar val_;
};




/** Functor to find value and location of minimum element greater than
 * a specified bound */
template <class Scalar>
class BoundedMinLocFunctor : public ReductionFunctorBase<Scalar>
{
public:
  BoundedMinLocFunctor(const MPIComm& comm, const Scalar& bound)
    : ReductionFunctorBase<Scalar>(comm), min_(), bound_(bound)
    {
      min_.what = HUGE_VAL;
    }
      

  void step(int i, const Scalar& x) const 
    {
      if (x < min_.what && x > bound_) 
      {
        min_.what = x;
        min_.where.setLocalIndex(i);
      }
    }

  void postProc() const 
    {
      IndexedValue<Scalar> out = min_;
      this->comm().allReduce(&min_, &out, 1, MPIComm::DOUBLE_INT, MPIComm::MINLOC);
      min_ = out;
    }

  IndexedValue<Scalar> result() const 
    {
      return min_;
    }

private:
  MPIComm comm_;
  mutable IndexedValue<Scalar> min_;
  Scalar bound_;
};




template <class Scalar>
class VectorFunctorTraits<Scalar, BoundedMinLocFunctor<Scalar> >
{
public:
  typedef IndexedValue<Scalar> ReturnType ;
};


}


#endif



