#ifndef SUNDANCE_OBJECTBUILDER_H
#define SUNDANCE_OBJECTBUILDER_H

#include "Sundance.hpp"
#include "SundanceVarMap.hpp"

namespace SundanceXML
{
  class ObjectBuilder
  {
  public:

    void interpret(const XMLObject& xml);

    const Mesh& mesh(const string& name) const {return mesh_.lookup(name);}

    const Expr& expr(const string& name) const {return expr_.lookup(name);}

    const QuadratureFamily& quad(const string& name) const 
    {return quadrature_.lookup(name);}

    const BasisFamily& basis(const string& name) const 
    {return basis_.lookup(name);}

    const DiscreteSpace& discreteSpace(const string& name) const 
    {return space_.lookup(name);}

    const CellFilter& cellFilter(const string& name) const 
    {return filter_.lookup(name);}

    const LinearProblem& linearProblem(const string& name) const 
    {return linearProb_.lookup(name);}

    Mesh createMesh(const XMLObject& xml)  ;

    Expr createExpr(const XMLObject& xml)  ;

    BasisFamily createBasis(const XMLObject& xml)  ;

    QuadratureFamily createQuadrature(const XMLObject& xml)  ;

    LinearProblem createLinearProblem(const XMLObject& xml)  ;

    DiscreteSpace createDiscreteSpace(const XMLObject& xml)  ;

    CellFilter createCellFilter(const XMLObject& xml)  ;

  private:
    
    VarMap<Mesh> mesh_;

    VarMap<Expr> expr_;

    VarMap<BasisFamily> basis_;

    VarMap<QuadratureFamily> quadrature_;

    VarMap<CellFilter> filter_;

    VarMap<DiscreteSpace> space_;

    VarMap<LinearProblem> linearProb_;

  };

}
#endif
