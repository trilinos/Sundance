/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_IFPACKOPERATOR_HPP
#define PLAYA_IFPACKOPERATOR_HPP


#include "PlayaEpetraMatrix.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"


namespace Playa
{
/**
 *
 */
class IfpackOperator : 
  public LinearOpWithSpaces<double>
{
public:
  /** */
  IfpackOperator(const EpetraMatrix* A,
    int fillLevels,
    int overlapFill,
    double relaxationValue,
    double relativeThreshold,
    double absoluteThreshold);

  /** 
   * Apply the operator. 
   */
  virtual void apply(
    Teuchos::ETransp applyType,
    const Vector<double>& in,
    Vector<double> out) const ;


private:
  RCP<Ifpack_IlukGraph> precondGraph_;

  RCP<Ifpack_CrsRiluk> precond_;

};


}

#endif 
