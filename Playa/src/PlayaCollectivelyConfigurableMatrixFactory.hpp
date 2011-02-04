/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_COLLECTIVELYCONFIGURABLEMATRIXFACTORY_HPP
#define PLAYA_COLLECTIVELYCONFIGURABLEMATRIXFACTORY_HPP

#include "PlayaDefs.hpp"
#include <vector>

namespace Playa
{
  /** 
   * Class CollectivelyConfigurableMatrixFactory provides an abstract 
   * interface for all-at-once matrix structure configuration
   */
  class CollectivelyConfigurableMatrixFactory
  {
  public:
    /** Virtual dtor */
    virtual ~CollectivelyConfigurableMatrixFactory(){;}


    /** */
    virtual void configure(int lowestRow,
                           const std::vector<int>& rowPtrs,
                           const std::vector<int>& nnzPerRow,
                           const std::vector<int>& data) = 0 ;
  private:
    
    
  };

}

#endif
