/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_IFPACKILUOPERATOR_HPP
#define PLAYA_IFPACKILUOPERATOR_HPP


#include "PlayaEpetraMatrix.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"


namespace Playa
{
/**
 *
 */
class IfpackILUOperator : 
  public LinearOpWithSpaces<double>
{
public:
  /** */
  IfpackILUOperator(const EpetraMatrix* A,
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
