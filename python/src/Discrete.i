// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceL2Projector.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceStdFwk
{
  class DiscreteSpace
  {
  public:
    /* */
    DiscreteSpace(const SundanceStdMesh::Mesh& mesh, const BasisFamily& basis,
                  const TSFExtended::VectorType<double>& vecType);
    /* */
    ~DiscreteSpace();

    
  };


  class L2Projector
  {
  public:
    /* */
    L2Projector(const DiscreteSpace& space, 
                const SundanceCore::Expr& expr);

    /* */
    ~L2Projector();

    /* */
    SundanceCore::Expr project() const ;
    
  };
}

%rename(DiscreteFunction) makeDiscreteFunction;

%inline %{
  /* Create a discrete function */
  SundanceCore::Expr makeDiscreteFunction(const SundanceStdFwk::DiscreteSpace& space,
                                          const TSFExtended::Vector<double>& vec)
  {
    return new SundanceStdFwk::DiscreteFunction(space, vec);
  }
  %}

%inline %{
  /* Create a discrete function */
  SundanceCore::Expr makeDiscreteFunction(const SundanceStdFwk::DiscreteSpace& space,
                                          const double& val)
  {
    return new SundanceStdFwk::DiscreteFunction(space, val);
  }
  %}


%inline %{
  /* Create a discrete function */
  SundanceCore::Expr makeDiscreteFunction(const SundanceStdFwk::DiscreteSpace& space,
                                          const TSFExtended::Vector<double>& vec,
                                          const std::string& name)
  {
    return new SundanceStdFwk::DiscreteFunction(space, vec, name);
  }
  %}

%inline %{
  /* Create a discrete function */
  SundanceCore::Expr makeDiscreteFunction(const SundanceStdFwk::DiscreteSpace& space,
                                          const double& val,
                                          const std::string& name)
  {
    return new SundanceStdFwk::DiscreteFunction(space, val, name);
  }
  %}




