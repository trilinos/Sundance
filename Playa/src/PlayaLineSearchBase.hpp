#ifndef PLAYA_LINESEARCHBASE_H
#define PLAYA_LINESEARCHBASE_H


#include "PlayaObjectiveBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Playa
{
  
using Teuchos::RCP;
using Teuchos::ParameterList;

/** */
enum LineSearchStatus {LS_Success, LS_ExceededMaxiters, LS_StepTooSmall, 
                       LS_Crashed};


/**
 * Base class for line search methods. 
 * \author Paul Boggs and Kevin Long
 */
class LineSearchBase : public ObjectWithVerbosity,
                       public Describable, 
                       public Printable
  
{
public:
  /** */
  LineSearchBase(const ParameterList& params); 


  /** */
  virtual LineSearchStatus search(const RCP<ObjectiveBase>& obj,
    const Vector<double>& x0,
    const double& f0,
    const Vector<double>& direction,
    const double& alphaMax,
    Vector<double>& xn, 
    Vector<double>& gradF,
    double& fVal) const = 0 ;
    
  /** */
  virtual double minStepSize() const {return minStepSize_;}

  /** */
  virtual int maxSteps() const {return maxSteps_;}

  /** */
  const ParameterList& params() const {return params_;}

private:
  ParameterList params_;

  int maxSteps_;

  double minStepSize_;
};
}

#endif
