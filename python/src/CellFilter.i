// -*- c++ -*-


%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceCellFilter.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "SundanceBoundaryCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "PySundanceCellPredicate.hpp"

  %}


// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"

namespace SundanceStdFwk
{

  class CellPredicate
  {
  public:
    CellPredicate();
    ~CellPredicate();


  };

  %extend CellPredicate
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

  namespace Internal
  {
    class CellSet
    {
    public:
      CellSet();
      ~CellSet();

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
  }

  class CellFilter
  {
  public:
    CellFilter();
    ~CellFilter();

    
    Internal::CellSet getCells(const SundanceStdMesh::Mesh& mesh) const ;
    int dimension(const SundanceStdMesh::Mesh& mesh) const ;

    CellFilter labeledSubset(int label) const ;

    CellFilter intersection(const CellFilter& other) const ;
  };

  %extend CellFilter
  {
    using namespace std;
    string __str__() 
    {
      string rtn; 
      stringstream os;
      self->print(os);
      rtn = os.str();
      return rtn;
    }

    CellFilter subset(PyObject* functor) const
    {
      Teuchos::RefCountPtr<SundanceStdFwk::CellPredicateFunctorBase> f 
        = Teuchos::rcp(new SundanceStdFwk::PySundanceCellPredicate(functor));
      CellPredicate p = new SundanceStdFwk::PositionalCellPredicate(f);
      return self->subset(p);
    }

    CellFilter __add__(const CellFilter& other) const
    {
      return self->operator+(other);
    }

    CellFilter __sub__(const CellFilter& other) const
    {
      return self->operator-(other);
    }
  }


  class CellFilterArray
  {
  public:
    CellFilterArray();
    ~CellFilterArray();

    unsigned int size() const ;

    void append(const CellFilter& x);
  };

}

%rename(MaximalCellFilter) makeMaximalCellFilter;
%rename(BoundaryCellFilter) makeBoundaryCellFilter;
%rename(DimensionalCellFilter) makeDimensionalCellFilter;
%rename(PositionalCellPredicate) makePyFunctorCellPredicate;

%inline %{
  /* Create a maximal cell filter */
  SundanceStdFwk::CellFilter makeMaximalCellFilter()
  {
    return new SundanceStdFwk::MaximalCellFilter();
  }
  %}


%inline %{
  /* Create a boundary cell filter */
  SundanceStdFwk::CellFilter makeBoundaryCellFilter()
  {
    return new SundanceStdFwk::BoundaryCellFilter();
  }
  %}

%inline %{
  /* Create a dimensional cell ftiler */
  SundanceStdFwk::CellFilter makeDimensionalCellFilter(int i)
  {
    return new SundanceStdFwk::DimensionalCellFilter(i);
  }
  %}

%inline %{
  /*  */
  SundanceStdFwk::CellPredicate makePyFunctorCellPredicate(PyObject* functor)
  {
    Teuchos::RefCountPtr<SundanceStdFwk::CellPredicateFunctorBase> f 
      = Teuchos::rcp(new SundanceStdFwk::PySundanceCellPredicate(functor));
    return new SundanceStdFwk::PositionalCellPredicate(f);
  }
  %}



%inline %{
  /*  */
  SundanceStdFwk::CellFilterArray CellFilterList()
  {
    return CellFilterArray();
  }
  /*  */
  SundanceStdFwk::CellFilterArray CellFilterList(const SundanceStdFwk::CellFilter& a)
  {
    return SundanceStdFwk::List(a);
  }

  /*  */
  SundanceStdFwk::CellFilterArray CellFilterList(const SundanceStdFwk::CellFilter& a,
                                 const SundanceStdFwk::CellFilter& b)
  {
    return SundanceStdFwk::List(a, b);
  }

  /*  */
  SundanceStdFwk::CellFilterArray CellFilterList(const SundanceStdFwk::CellFilter& a,
                                 const SundanceStdFwk::CellFilter& b,
                                 const SundanceStdFwk::CellFilter& c)
  {
    return SundanceStdFwk::List(a, b, c);
  }

  /*  */
  SundanceStdFwk::CellFilterArray CellFilterList(const SundanceStdFwk::CellFilter& a,
                                 const SundanceStdFwk::CellFilter& b,
                                 const SundanceStdFwk::CellFilter& c,
                                 const SundanceStdFwk::CellFilter& d)
  {
    return SundanceStdFwk::List(a, b, c, d);
  }

  /*  */
  SundanceStdFwk::CellFilterArray CellFilterList(const SundanceStdFwk::CellFilter& a,
                                 const SundanceStdFwk::CellFilter& b,
                                 const SundanceStdFwk::CellFilter& c,
                                 const SundanceStdFwk::CellFilter& d,
                                 const SundanceStdFwk::CellFilter& e)
  {
    return SundanceStdFwk::List(a, b, c, d, e);
  }

  %}


