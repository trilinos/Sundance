#ifndef PLAYA_LINE_SEARCH_BASED_OPT_BASE_H
#define PLAYA_LINE_SEARCH_BASED_OPT_BASE_H

#include "PlayaUnconstrainedOptimizerBase.hpp"
#include "PlayaLineSearchBase.hpp"
#include "PlayaOptConvergenceTestBase.hpp"
#include "PlayaObjectiveBase.hpp"

namespace Playa
{

/** 
 * 
 */
class DirectionGeneratorBase : public ObjectWithVerbosity
{
public:
  /** */
  DirectionGeneratorBase() {}

  /** */
  ~DirectionGeneratorBase() {}

  /** */
  virtual bool generateDirection(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xCur,
    const Vector<double>& gradCur,
    const double& fCur,
    Vector<double>& p) = 0 ;
};


/** 
 * Base class for optimizers based on line search methods.
 *
 * @author Kevin Long
 */
class LineSearchBasedOptBase : public UnconstrainedOptimizerBase 
{
public:
  /** */
  LineSearchBasedOptBase(const ParameterList& params);
  /** */
  virtual ~LineSearchBasedOptBase(){}
        
  /** Main method to apply the algorithm starting with x and
      returning the result in x */
  OptState run(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xInit,
    const RCP<ConvergenceMonitor>& convMonitor = null) const ;

  /** */
  virtual RCP<DirectionGeneratorBase> makeDirectionGenerator() const = 0 ; 

  /** */
  virtual void print(std::ostream& os) const ;

protected:
  /** */
  const RCP<LineSearchBase>& lineSearch() const {return lineSearch_;}

  /** */
  const RCP<OptConvergenceTestBase>& convTest() const {return convTest_;}

private:
  RCP<LineSearchBase> lineSearch_;
  RCP<OptConvergenceTestBase> convTest_;
  
};


}

#endif
