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

    int dim() const ;
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

    Vector<Scalar> copy() const ;

    Vector<Scalar> acceptCopyOf(const Vector<Scalar>& x);

    Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

    Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

    Vector<Scalar> reciprocal() const ;

    Vector<Scalar> abs() const ;

    void setToConstant(const Scalar& alpha) ;

    Scalar norm1() const ;

    Scalar norm2() const ;

    Scalar normInf() const ;

    void zero();

    Scalar max() const;

    Scalar max(int& index)const;

    Scalar max(const Scalar& bound, int& index)const;

    Scalar min()const;

    Scalar min(int& index)const;

    Scalar min(const Scalar& bound, int& index)const;
    
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

      Vector<Scalar> __add__(const Vector<Scalar>& other) 
      {
        return (*self) + other;
      }

      Vector<Scalar> __sub__(const Vector<Scalar>& other) 
      {
        return (*self) - other;
      }

      Vector<Scalar> __mul__(const Scalar& other) 
      {
        return (*self) * other;
      }

      Vector<Scalar> __div__(const Scalar& other) 
      {
        return (*self) *(1.0/ other);
      }

      Vector<Scalar> __rmul__(const Scalar& other) 
      {
        return (*self) * other;
      }

      Scalar __mul__(const Vector<Scalar>& other) 
      {
        return (*self) * other;
      }

      
      Scalar __getitem__(int globalIndex) const 
      {
        return self->getElement(globalIndex);
      }
      
      void __setitem__(int globalIndex, const Scalar& value)
      {
        return self->setElement(globalIndex, value);
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
