#ifndef PYSUNDANCENOXSOLVERHANDLE_H
#define PYSUNDANCENOXSOLVERHANDLE_H

#include "TSFNOXSolver.H"


namespace TSFExtended
{
  class NOXSolverHandle

  {
  public:
    /** */
    NOXSolverHandle();

    /** */
    NOXSolverHandle(const Teuchos::RefCountPtr<NOXSolver>& ptr)
      : ptr_(ptr) {;}

    /** */
    NOX::StatusTest::StatusType solve() const {return ptr_->solve();}

  private:
    Teuchos::RefCountPtr<NOXSolver> ptr_;
  };
}
#endif // 
