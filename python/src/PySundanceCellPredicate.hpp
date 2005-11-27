#ifndef PYSUNDANCECELLPREDICATE_H
#define PYSUNDANCECELLPREDICATE_H

#include "Python.h"
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

  protected:
    PyObject * setEvalOp(PyObject *);

  private:
    // Private and not implemented
    //PyInterface();
    PySundanceCellPredicate(const PySundanceCellPredicate &);
    PySundanceCellPredicate & operator=(const PySundanceCellPredicate &);

  private:
    PyObject* py_functor_;
    mutable PySundanceCallback  py_evalOp_;
  };
}
#endif // 
