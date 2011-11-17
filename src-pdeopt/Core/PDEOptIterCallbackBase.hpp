#ifndef PDEOPT_ITER_CALLBACK_BASE_H
#define PDEOPT_ITER_CALLBACK_BASE_H

#include "SundanceDefs.hpp"
#include <string>

namespace Sundance
{

class PDEConstrainedObjBase;

/**
 * IterCallbackBase provides an abstract interface for callbacks to be
 * performed at the start of each optimization iteration.
 */
class IterCallbackBase 
{
public:
  /** */
  virtual ~IterCallbackBase(){}

  /** */
  virtual void call(const PDEConstrainedObjBase* obj, int iter) const = 0 ;
};



/**
 * Default callback writes all variables to disk at a specified frequency
 */
class DefaultIterCallback : public IterCallbackBase
{
public:

  /** */
  DefaultIterCallback(
    const std::string& filename, 
    const std::string& type,
    int frequency=1);

  /** */
  virtual ~DefaultIterCallback(){}

  /** */
  void call(const PDEConstrainedObjBase* obj, int iter) const ;

private:
  std::string type_;
  std::string filename_;
  int frequency_;
};




  
}

#endif
