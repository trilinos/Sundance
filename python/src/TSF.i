// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFVector.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFNOXSolver.H"
#include "TSFLinearSolverBuilder.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PyTeuchos_Utils.hpp"
#include "PySundanceNOXSolverHandle.hpp"
  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

%rename(Vector) Vec;
%rename(VectorType) VecType;
%rename(VectorSpace) VecSpace;
%rename(LinearOperator) LinOp;
%rename(LinearSolver) LinSol;
%rename(NonlinearOperator) NonlinOp;
%rename(NOXSolver) makeNOXSolver;



 /* --------- vector space ------------ */
namespace TSFExtended
{
  template <class Scalar> class Vector;
  template <class Scalar>
  class VectorSpace
  {
  public:
    Vector<Scalar> createMember();
  };

  %template(VecSpace) VectorSpace<double>;
}

/* --------- vector ------------ */
namespace TSFExtended
{
  template <class Scalar> class Vector
  {
  public:
    Vector();
    ~Vector();

    void setElement(int globalIndex, const Scalar& x); 

    VectorSpace<Scalar> space() const ;

    %extend 
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(Vec) Vector<double>;

}

 /* --------- vector space ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class VectorSpace
  {
  public:
    VectorSpace();
    ~VectorSpace();

    Vector<Scalar> createMember();

    %extend 
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(VecSpace) VectorSpace<double>;
}



/* --------- vector type ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class VectorType
  {
  public:
    VectorType();
    ~VectorType();

    VectorSpace<Scalar> createEvenlyPartitionedSpace(int nLocal) const ;

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(VecType) VectorType<double>;

}

%rename(EpetraVectorType) makeEpetraVectorType;

%inline %{
  /* Create an epetra vector type */
  TSFExtended::VectorType<double> makeEpetraVectorType()
  {
    return TSFExtended::VectorType<double>(new TSFExtended::EpetraVectorType());
  }
  %}


/* --------- linear operator ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class LinearOperator
  {
  public:
    LinearOperator();
    ~LinearOperator();

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(LinOp) LinearOperator<double>;

}

/* --------- linear solver ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class LinearSolver
  {
  public:
    LinearSolver();
    ~LinearSolver();

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(LinSol) LinearSolver<double>;

}

%inline %{
  /* Read a linear solver from an XML file */
  TSFExtended::LinearSolver<double> readSolver(const std::string& filename)
  {
    Teuchos::ParameterXMLFileReader reader(filename);
    Teuchos::ParameterList solverParams = reader.getParameters();
    TSFExtended::LinearSolver<double> solver 
      = TSFExtended::LinearSolverBuilder::createSolver(solverParams);
    return solver;
  }
  %}

%inline %{
  /* Read a linear solver from a parameter list */
  TSFExtended::LinearSolver<double> buildSolver(const Teuchos::ParameterList& params)
  {
    TSFExtended::LinearSolver<double> solver ;
    try
      {
        solver = TSFExtended::LinearSolverBuilder::createSolver(params);
      }
    catch(std::exception& e)
      {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        std::cerr << "detected exception "
                  << e.what() << " in buildSolver()" << std::endl;
      }
    return solver;
  }
  %}




/* --------- nonlinear operator ------------ */
namespace TSFExtended
{
  template <class Scalar>
  class NonlinearOperator
  {
  public:
    NonlinearOperator();
    ~NonlinearOperator();

    %extend
    {
      std::string __str__() 
      {
        std::string rtn; 
        std::stringstream os;
        self->print(os);
        rtn = os.str();
        return rtn;
      }
    }
  };

  %template(NonlinOp) NonlinearOperator<double>;

}

namespace NOX
{
  namespace StatusTest
  {
    enum StatusType {Unevaluated, Unconverged, Converged, Failed};
  }
}


/* --------- NOX solver ------------ */
namespace TSFExtended
{
  class NOXSolverHandle
  {
  public:
    ~NOXSolverHandle();
    NOXSolverHandle();

    NOX::StatusTest::StatusType solve() const ;

  };
}





%inline %{
  /* Read a linear solver from an XML file */
  TSFExtended::
    NOXSolverHandle makeNOXSolver(const Teuchos::ParameterList& params,
                                  const TSFExtended::NonlinearOperator<double>& F)
  {

    NOXSolverHandle rtn;
    try
      {
        rtn = rcp(new TSFExtended::NOXSolver(params, F));
      }
    catch(std::exception& e)
      {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        std::cerr << "detected exception "
                  << e.what() << " in makeNOXSolver()" << std::endl;
      }
    return rtn;
  }
  %}
