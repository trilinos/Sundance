#ifndef PLAYA_STEEPESTDESCENT_H
#define PLAYA_STEEPESTDESCENT_H


#include "PlayaLineSearchBasedOptBase.hpp"


namespace Playa
{
/** 
 * Implements the steepest descent method, used here mainly for debugging
 * objective functions and gradients. 
 *
 * @author Kevin Long
 */
class SteepestDescent : public LineSearchBasedOptBase
{
public:
  /** */
  SteepestDescent(const ParameterList& params);
        
  /** */
  std::string description() const {return "SteepestDescent";}

  /** */
  RCP<DirectionGeneratorBase> makeDirectionGenerator() const ; 
};

/** */
class SteepestDescentDirection : public DirectionGeneratorBase
{
public:
  /** */
  SteepestDescentDirection() {}

  /** */
  ~SteepestDescentDirection() {}

  /** */
  bool generateDirection(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xCur,
    const Vector<double>& gradCur,
    const double& fCur,
    Vector<double>& p) ;
};

}

#endif
