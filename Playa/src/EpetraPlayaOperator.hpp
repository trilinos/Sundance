/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef EPETRA_PLAYA_OPERATOR_HPP
#define EPETRA_PLAYA_OPERATOR_HPP

#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "Epetra_Operator.h"

namespace Epetra
{
using namespace Teuchos;
using namespace Playa;
  
  

/** */
class Epetra_PlayaOperator : public Epetra_Operator
{
public:
  /** */
  Epetra_PlayaOperator(const LinearOperator<double>& A,
    const LinearSolver<double>& solver=LinearSolver<double>());
    
  /** */
  int SetUseTranspose(bool useTrans) {useTranspose_ = useTrans; return 0;}

  /** */
  int Apply(const Epetra_MultiVector& in, Epetra_MultiVector& out) const ;

  /** */
  int ApplyInverse(const Epetra_MultiVector& in, Epetra_MultiVector& out) const ;

  /** */
  double NormInf() const ;

  /** */
  const char* Label() const ;

  /** */
  bool UseTranspose() const {return useTranspose_;}

  /** */
  bool HasNormInf() const {return false;}

  /** */
  const Epetra_Comm& Comm() const {return *comm_;}

  /** */
  const Epetra_Map& OperatorDomainMap() const {return *domain_;}

  /** */
  const Epetra_Map& OperatorRangeMap() const {return *range_;}

    

private:
  LinearOperator<double> A_;
  LinearSolver<double> solver_;
  bool useTranspose_;
  RCP<Epetra_Comm> comm_;
  RCP<const Epetra_Map> domain_;
  RCP<const Epetra_Map> range_;
  bool isNativeEpetra_;
  bool isCompoundEpetra_;
  std::string label_;
};
}

#endif
