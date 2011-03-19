#ifndef PLAYA_CONVERGENCE_MONITOR_H
#define PLAYA_CONVERGENCE_MONITOR_H

#include "Teuchos_Array.hpp"

namespace Playa
{
using Teuchos::Array;

/** 
 *
 */
class ConvergenceMonitor
{
public:
  /** */
  ConvergenceMonitor();

  /** */
  void addRecord(int iter, const Array<double>& record) const ;

  /** */
  void reset() ;

  /** */
  void write(std::ostream& os) const ;

private:
  mutable Array<int> iters_;
  mutable Array<Array<double> > records_;
};


}


#endif
