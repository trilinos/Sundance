// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceLinearProblem.hpp"
#include "PySundanceLinearSolver.hpp"

  %}

namespace SundanceStdFwk
{
class Block
  {
  public:
    /** */
    Block(const SundanceCore::Expr& expr, const TSFExtended::VectorType<double>& vecType);

    /** */
    const SundanceCore::Expr& expr() const ;

    /** */
    const TSFExtended::VectorType<double>& vecType() const ;
  };

  class BlockArray
  {
  public:
    BlockArray(int n);
  };

  %extend Block
  {
    using namespace std;
    string __str__() 
    {
      string rtn; 
      stringstream os;
      os << *self;
      rtn = os.str();
      return rtn;
    }
  }
    
  %extend BlockArray
  {
    using namespace std;
    string __str__() 
    {
      string rtn; 
      stringstream os;
      os << *self;
      rtn = os.str();
      return rtn;
    }
  }
}


%inline %{
/* */
  SundanceStdFwk::BlockArray 
    BlockList(const SundanceStdFwk::Block& a)
  {
    return tuple(a);
  }

  /* */
  SundanceStdFwk::BlockArray 
    BlockList(const SundanceStdFwk::Block& a,
              const SundanceStdFwk::Block& b)
  {
    return tuple(a,b);
  }

  /* */
  SundanceStdFwk::BlockArray 
    BlockList(const SundanceStdFwk::Block& a,
              const SundanceStdFwk::Block& b,
              const SundanceStdFwk::Block& c)
  {
    return tuple(a,b,c);
  }

  /* */
  SundanceStdFwk::BlockArray 
    BlockList(const SundanceStdFwk::Block& a,
              const SundanceStdFwk::Block& b,
              const SundanceStdFwk::Block& c,
              const SundanceStdFwk::Block& d)
  {
    return tuple(a,b,c,d);
  }

  /* */
  SundanceStdFwk::BlockArray 
    BlockList(const SundanceStdFwk::Block& a,
              const SundanceStdFwk::Block& b,
              const SundanceStdFwk::Block& c,
              const SundanceStdFwk::Block& d,
              const SundanceStdFwk::Block& e)
  {
    return tuple(a,b,c,d,e);
  }
                                      
  %}

// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace SundanceStdFwk
{
  
  typedef RefCountPtr<DOFMapBase> DOFMap;

  class LinearProblem
  {
  public:
    LinearProblem(const SundanceStdMesh::Mesh& mesh, 
                  const SundanceCore::Expr& eqn,
                  const SundanceCore::Expr& bc,
                  const SundanceCore::Expr& v, 
                  const SundanceCore::Expr& u,
                  const TSFExtended::VectorType<double>& vecType);

    LinearProblem(const SundanceStdMesh::Mesh& mesh, 
                  const SundanceCore::Expr& eqn,
                  const SundanceCore::Expr& bc,
                  const SundanceStdFwk::BlockArray& v, 
                  const SundanceStdFwk::BlockArray& u);

    TSFExtended::Vector<double> getSingleRHS() const ;

    TSFExtended::LinearOperator<double> getOperator() const ;

    SundanceCore::Expr solve(const TSFExtended::LinearSolver<double>& solver) const ;

  };

  %extend LinearProblem {
    SundanceCore::Expr solve(PyObject* pySolver)
    {
      TSFExtended::PySundanceLinearSolver* tmp  
        = new TSFExtended::PySundanceLinearSolver(pySolver);
      
      Teuchos::RefCountPtr<TSFExtended::LinearSolverBase<double> > 
        r = rcp(tmp);

      TSFExtended::LinearSolver<double> cxxSolver = r;
      return self->solve(cxxSolver);
    }
  }

  %extend LinearProblem {
    void printRowMap() const 
    {
      for (int b=0; b<self->numBlockRows(); b++)
        {
          self->rowMap(b)->print(cout);
        }
    }
  }

  %extend LinearProblem {
    void printColMap() const 
    {
      for (int b=0; b<self->numBlockCols(); b++)
        {
          self->colMap(b)->print(cout);
        }
    }
  }
 
}
