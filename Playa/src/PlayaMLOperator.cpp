/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaMLOperator.hpp"
#include "Teuchos_Array.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaEpetraVector.hpp"



#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Playa;
using namespace Teuchos;

MLOperator::MLOperator(
  const LinearOperator<double>& op,
  const ParameterList& mlParams)
  : LinearOpWithSpaces<double>(op.domain(), op.range()),
    mlPrec_()
{
	Epetra_CrsMatrix& A = EpetraMatrix::getConcrete(op);
  
  
  mlPrec_ = rcp(new ML_Epetra::MultiLevelPreconditioner(A, mlParams));
}


void MLOperator::apply(
  Teuchos::ETransp applyType,
  const Vector<double>& in,
  Vector<double> out) const
{
  /* grab the epetra vector objects underlying the input and output vectors */
  const Epetra_Vector& epIn = EpetraVector::getConcrete(in);
  Epetra_Vector& epOut = EpetraVector::getConcrete(out);


  int ierr;

  /* do the solve (or transpose solve) */
  if (applyType==NO_TRANS)
    {
      ierr = mlPrec_->ApplyInverse(epIn, epOut);
    }
  else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(applyType != NO_TRANS, std::runtime_error,
        "ML preconditioner does not support transposes");
    }
  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
			     "Playa::MLOperator apply: "
			     "mlPrec_->ApplyInvers() failed with ierr="
			     << ierr);
}
  
