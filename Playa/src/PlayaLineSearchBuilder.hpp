#ifndef PLAYA_LINE_SEARCH_BUILDER_H
#define PLAYA_LINE_SEARCH_BUILDER_H

#include "PlayaLineSearchBase.hpp"
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
 * Builder class for line search objects
 *
 * @author Kevin Long
 */
class LineSearchBuilder
{
public:
  /** Construct a line search object from a parameter specification */
  static RCP<LineSearchBase> createLineSearch(const ParameterList& params,
    int verb=0);
};
}

#endif
