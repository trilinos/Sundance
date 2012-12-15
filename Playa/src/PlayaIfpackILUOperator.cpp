/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaIfpackILUOperator.hpp"
#include "Teuchos_Array.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaEpetraVector.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

namespace Playa
{

using namespace Teuchos;

IfpackILUOperator::IfpackILUOperator(const EpetraMatrix* A,
  int fillLevels,
  int overlapFill,
  double relaxationValue,
  double relativeThreshold,
  double absoluteThreshold)
  : LinearOpWithSpaces<double>(A->domain(), A->range()),
    precondGraph_(),
    precond_()
{
  const Epetra_CrsMatrix* matrix = A->crsMatrix();

  const Epetra_CrsGraph& matrixGraph = matrix->Graph();
			
  Ifpack_IlukGraph* precondGraph 
    = new Ifpack_IlukGraph(matrixGraph, fillLevels, overlapFill);
  precondGraph_ = rcp(precondGraph);

  int ierr = precondGraph->ConstructFilledGraph();

  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
    "IfpackOperator ctor: "
    "precondGraph->ConstructFilledGraph() failed with ierr="
    << ierr);

  Ifpack_CrsRiluk* precond = new Ifpack_CrsRiluk(*precondGraph);
  precond_ = rcp(precond);

  precond->SetRelaxValue(relaxationValue);
  precond->SetRelativeThreshold(relativeThreshold);
  precond->SetAbsoluteThreshold(absoluteThreshold);

  ierr = precond->InitValues(*matrix);

  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
    "IfpackOperator ctor: "
    "precond->InitValues() failed with ierr="
    << ierr);

  ierr = precond->Factor();

  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
    "IfpackOperator ctor: "
    "precond->Factor() failed with ierr="
    << ierr);
}

void IfpackILUOperator::apply(Teuchos::ETransp transApplyType,
  const Vector<double>& in,
  Vector<double> out) const
{
  /* grab the epetra vector objects underlying the input and output vectors */
  const Epetra_Vector& epIn = EpetraVector::getConcrete(in);
  Epetra_Vector& epOut = EpetraVector::getConcrete(out);

  /* ifpack's solve is logically const but declared non-const because
   *  internal data changes. So, do a const_cast. */
  Ifpack_CrsRiluk* p = const_cast<Ifpack_CrsRiluk*>(precond_.get());

  int ierr;

  /* do the solve (or transpose solve) */
  if (transApplyType==NO_TRANS)
  {
    ierr = p->Solve(false, epIn, epOut);
  }
  else
  {
    ierr = p->Solve(true, epIn, epOut);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error,
			     "Playa::IfpackILUOperator apply: "
			     "p->Solve() (transpose state="
			     << transApplyType << ") failed with ierr="
			     << ierr);
}
  
}
