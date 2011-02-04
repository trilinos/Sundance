/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_SOLVERSTATE_HPP
#define PLAYA_SOLVERSTATE_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
using namespace Teuchos;

/** 
 *
 */
enum SolverStatusCode {SolveCrashed, SolveFailedToConverge, SolveConverged};


/**
 *
 */
template <class Scalar>
class SolverState
{
public:
  /** */
  SolverState(SolverStatusCode finalState, const std::string& msg, 
    int finalIters, const Scalar& finalResid)
    : finalState_(finalState),
      finalResid_(finalResid),
      finalIters_(finalIters),
      msg_(msg)
    {;}

  /** */
  SolverState() {;}

  /** */
  const Scalar& finalResid() const {return finalResid_;}

  /** */
  int finalIters() const {return finalIters_;}

  /** */
  const SolverStatusCode& finalState() const {return finalState_;}

  /** */
  const std::string& finalMsg() const {return msg_;}

  /** */
  std::string stateDescription() const 
    {
      switch (finalState_)
      {
        case SolveCrashed:
          return "Crashed";
        case SolveFailedToConverge:
          return "Failed to converge";
        case SolveConverged:
          return "Converged";
      }
      return "Crashed";
    }

private:
    
  SolverStatusCode finalState_;

  Scalar finalResid_;

  int finalIters_;

  std::string msg_;
};


template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, 
  const Playa::SolverState<Scalar>& state)
{
  os << "Solver final state: " << state.stateDescription() << std::endl;
  os << "message: " << state.finalMsg() << std::endl;
  os << "iters taken: " << state.finalIters() << std::endl;
  os << "final residual: " << state.finalResid() << std::endl;
  return os;
}

}


#endif
