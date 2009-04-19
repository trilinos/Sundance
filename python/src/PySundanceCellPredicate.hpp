#ifndef PYSUNDANCECELLPREDICATE_H
#define PYSUNDANCECELLPREDICATE_H

#include "PySundanceCallback.hpp"
#include "SundancePositionalCellPredicate.hpp"


namespace SundanceStdFwk
{
  class PySundanceCellPredicate : public CellPredicateFunctorBase

  {
  public:
    /** */
    PySundanceCellPredicate(PyObject* functor = Py_None);
    /** */
    ~PySundanceCellPredicate();

    /* */
    GET_RCP(CellPredicateFunctorBase);

    /** */
    bool operator()(const SundanceUtils::Point& x) const ;

    /** */
    string description() const ;

  protected:
    PyObject* setEvalOp(PyObject* pyClass);
    PyObject* setDescr(PyObject* pyClass);

  private:
    // Private and not implemented
    //PyInterface();
    PySundanceCellPredicate(const PySundanceCellPredicate &);
    PySundanceCellPredicate & operator=(const PySundanceCellPredicate &);

  private:
    PyObject* py_functor_;
    mutable RefCountPtr<PySundanceCallback>  evalOpCallback_;
    mutable RefCountPtr<PySundanceCallback>  descrCallback_;
  };
}
#endif // 
