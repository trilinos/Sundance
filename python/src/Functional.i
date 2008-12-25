// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceFunctional.hpp"

  %}


namespace SundanceStdFwk
{
class Functional
  {
  public:
    /** */
    Functional(){;}

    /** */
    Functional(const SundanceStdMesh::Mesh& mesh, 
      const SundanceCore::Expr& integral, 
      const TSFExtended::VectorType<double>& vecType);

    /** */
    Functional(const SundanceStdMesh::Mesh& mesh, 
      const SundanceCore::Expr& integral, 
      const SundanceCore::Expr& essentialBC,
      const TSFExtended::VectorType<double>& vecType);

    /** */
    TSFExtended::NonlinearOperator<double>
    nonlinearVariationalProb(const SundanceCore::Expr& var,
                             const SundanceCore::Expr& varEvalPts,
                             const SundanceCore::Expr& unk,
                             const SundanceCore::Expr& unkEvalPts,
                             const SundanceCore::Expr& fixed,
                             const SundanceCore::Expr& fixedEvalPts) const ;


    /** */
    FunctionalEvaluator evaluator(const SundanceCore::Expr& var,
                                  const SundanceCore::Expr& varEvalPts,
                                  const SundanceCore::Expr& fixed,
                                  const SundanceCore::Expr& fixedEvalPts) const ;


    /** */
    FunctionalEvaluator evaluator(const SundanceCore::Expr& var,
                                  const SundanceCore::Expr& varEvalPts) const ;

    /** */
    const SundanceStdMesh::Mesh& mesh() const ;
};

class FunctionalEvaluator 
  : public TSFExtended::ParameterControlledObjectWithVerbosity<FunctionalEvaluator>
{
public:
  /** */
  FunctionalEvaluator();

  /** */
  FunctionalEvaluator(const SundanceStdMesh::Mesh& mesh, 
    const SundanceCore::Expr& integral,
    const Teuchos::ParameterList& verbParams = *defaultVerbParams());
  /** */
  FunctionalEvaluator(const SundanceStdMesh::Mesh& mesh, 
    const SundanceCore::Expr& integral,
    const SundanceCore::Expr& bcs,
    const SundanceCore::Expr& var,
    const SundanceCore::Expr& varEvalPts,
    const TSFExtended::VectorType<double>& vectorType,
    const Teuchos::ParameterList& verbParams = *defaultVerbParams());
  /** */
  FunctionalEvaluator(const SundanceStdMesh::Mesh& mesh, 
    const SundanceCore::Expr& integral,
    const SundanceCore::Expr& bcs,
    const SundanceCore::Expr& vars,
    const SundanceCore::Expr& varEvalPts,
    const SundanceCore::Expr& fields,
    const SundanceCore::Expr& fieldValues,
    const TSFExtended::VectorType<double>& vectorType,
    const Teuchos::ParameterList& verbParams = *defaultVerbParams());


  /** */
  double evaluate() const ;

  /** */
  SundanceCore::Expr evalGradient(double& value) const ;

  /** */
  double fdGradientCheck(double h) const ;
};


}
