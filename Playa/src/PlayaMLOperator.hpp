/* @HEADER@ */
//   
 /* @HEADER@ */

#ifndef PLAYA_MLOPERATOR_HPP
#define PLAYA_MLOPERATOR_HPP


#include "PlayaEpetraMatrix.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"

#include "ml_MultiLevelPreconditioner.h"
#include "EpetraPlayaOperator.hpp"

namespace Playa
{
/**
 *
 */
class MLOperator :
  public LinearOpWithSpaces<double>
{
public:
  /** */
  MLOperator(
    const LinearOperator<double>& op,
    const ParameterList& mlParams);

  /** 
   * Apply the operator. 
   * 
   * \param applyType Indicates whether to apply the operator, its transpose,
   * or its conjugate transpose. 
   * \param in The vector on which the operator is to act
   * \param out The vector into which the result of the operation 
   * is to be written. This vector should already be initialized by the
   * appropriate space.
   **/
  virtual void apply(
    Teuchos::ETransp applyType,
    const Vector<double>& in,
    Vector<double> out) const ;


private:
  RCP<ML_Epetra::MultiLevelPreconditioner> mlPrec_;

};
}

#endif /* PlayaIFPACKOPERATOR_HPP */
