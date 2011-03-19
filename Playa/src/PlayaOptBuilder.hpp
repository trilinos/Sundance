#ifndef PLAYA_OPT_BUILDER_H
#define PLAYA_OPT_BUILDER_H

#include "PlayaUnconstrainedOptimizerBase.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos
{
class ParameterList;
}

namespace Playa
{
using Teuchos::ParameterList;
using Teuchos::RCP;

/** 
 * Builder class for optimizer objects
 *
 * @author Kevin Long
 */
class OptBuilder
{
public:
  /** Construct an optimizer object from a parameter specification */
  static RCP<UnconstrainedOptimizerBase> 
  createOptimizer(const ParameterList& params,
    int verb=0);


  /** Construct an optimizer object from an XML file */
  static RCP<UnconstrainedOptimizerBase> 
  createOptimizer(const string& filename,
    int verb=0);
};
}

#endif
