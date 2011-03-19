#ifndef PLAYA_UNCONSTRAINED_OPTIMIZER_BASE_H
#define PLAYA_UNCONSTRAINED_OPTIMIZER_BASE_H

#include "PlayaObjectiveBase.hpp"
#include "PlayaOptState.hpp"
#include "PlayaConvergenceMonitor.hpp"
#include "Teuchos_ENull.hpp"

namespace Playa
{

/** 
 * Base class for unconstrained optimizers
 *
 * @author Kevin Long
 */
class UnconstrainedOptimizerBase : public ObjectWithVerbosity,
                                   public Describable,
                                   public Printable
{
public:
  /** */
  UnconstrainedOptimizerBase(){}
  /** */
  virtual ~UnconstrainedOptimizerBase(){}
        

  /** Main method to apply the algorithm starting with x and
      returning the result in x */
  virtual OptState run(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xInit,
    const RCP<ConvergenceMonitor>& convMonitor = null) const = 0 ;



};
}

#endif
