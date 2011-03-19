#ifndef PLAYA_SIMPLEBACKTRACKING_H
#define PLAYA_SIMPLEBACKTRACKING_H


#include "PlayaLineSearchBase.hpp"

namespace Playa
{
/**
 * Simple backtracking line search
 */
class SimpleBacktracking : public LineSearchBase
{
public:
  /** */
  SimpleBacktracking(const ParameterList& params); 

  /** */
  LineSearchStatus search(const RCP<ObjectiveBase>& obj,
    const Vector<double>& x0,
    const double& f0,
    const Vector<double>& direction,
    const double& alphaMax,
    Vector<double>& xn, 
    Vector<double>& gradF,
    double& fVal) const ;
    
  /** */
  std::string description() const ;

  /** */
  void print(std::ostream& os) const 
    {os << description();}
private:
    
};
}

#endif
