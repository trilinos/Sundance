#ifndef PLAYA_BASIC_LMBFGS_OPT_H
#define PLAYA_BASIC_LMBFGS_OPT_H


#include "PlayaLineSearchBasedOptBase.hpp"
#include <deque>

namespace Playa
{
/** 
 * Implements the L-BFGS algorithm
 *
 * @author Kevin Long
 */
class BasicLMBFGS : public LineSearchBasedOptBase
{
public:
  /** */
  BasicLMBFGS(const ParameterList& params);
        
  /** */
  std::string description() const {return "LM-BFGS";}

  /** */
  RCP<DirectionGeneratorBase> makeDirectionGenerator() const ; 

private:
  int memSize_;
};

/** */
class BasicLMBFGSDirection : public DirectionGeneratorBase
{
public:
  /** */
  BasicLMBFGSDirection(int memSize);

  /** */
  ~BasicLMBFGSDirection() {}

  /** */
  bool generateDirection(const RCP<ObjectiveBase>& obj,
    const Vector<double>& xCur,
    const Vector<double>& gradCur,
    const double& fCur,
    Vector<double>& p) ;

private:
  int memSize_;
  Vector<double> xPrev_;
  Vector<double> gradPrev_;
  std::deque<Vector<double> > sMem_;
  std::deque<Vector<double> > yMem_;
};

}

#endif
