#ifndef PLAYA_OPT_STATE_H
#define PLAYA_OPT_STATE_H


#include "PlayaVectorDecl.hpp"

namespace Playa
{

/** 
 * OptStatus provides diagnostic information on the current state
 * of an optimization run.
 */
enum OptStatus
{
  Opt_Continue,
  Opt_Converged,
  Opt_DirectionFailure,
  Opt_ExceededMaxiters,
  Opt_LineSearchFailed,
  Opt_Crashed
};

/** \relates OptStatus */
inline std::ostream& operator<<(std::ostream& os, const OptStatus& s)
{
  switch (s)
  {
    case Opt_Continue:
      os << "Opt_Continue"; break;
    case Opt_Converged:
      os << "Opt_Converged"; break;
    case Opt_DirectionFailure:
      os << "Opt_DirectionFailure"; break;
    case Opt_ExceededMaxiters:
      os << "Opt_ExceededMaxiters"; break;
    case Opt_LineSearchFailed:
      os << "Opt_LineSearchFailed"; break;
    default:
      os << "Opt_Crashed";
  }
  return os;
}


/** 
 * OptState encapsulates the current state of an optimization run, for
 * use in convergence testing.
 */
class OptState
{
public:
  /** */
  OptState(const Vector<double>& xCur,
    const double& fCur,
    const Vector<double>& gradCur);

  /** */
  OptStatus status() const {return status_;}

  /** */
  void setStatus(const OptStatus status) {status_ = status;} 

  /** Return the current iteration count */
  int iter() const {return iter_;}

  /** Return the current objective function value */
  double fCur() const {return fCur_;}

  /** Return the previous objective function value */
  double fPrev() const {return fPrev_;}

  /** Return the current evaluation point */
  Vector<double> xCur() const {return xCur_;}

  /** Return the previous evaluation point */
  Vector<double> xPrev() const {return xPrev_;}

  /** Return the current gradient */
  Vector<double> gradCur() const {return gradCur_;}

  /** Return the previous gradientx */
  Vector<double> gradPrev() const {return gradPrev_;}

  /** */
  void update(const Vector<double>& xNew, const Vector<double>& gradNew, 
    const double& fNew);

private:
  OptStatus status_;
  int iter_;

  Vector<double> xCur_;
  Vector<double> xPrev_;

  Vector<double> gradCur_;
  Vector<double> gradPrev_;

  double fCur_;
  double fPrev_;

};

}

#endif

