#include "PySundanceCellPredicate.hpp"

using namespace SundanceStdFwk;
using namespace SundanceUtils;

PySundanceCellPredicate::PySundanceCellPredicate(PyObject* functor) 
  : py_functor_(functor), py_evalOp_()

{
  // Increment the reference count
  Py_XINCREF(py_functor_);

  // If the python object has a "evalOp" attribute, set it
  // to be the PySundanceCellPredicate computeF callback function
  if (PyObject_HasAttrString (py_functor_,
			      "evalOp")) {
    setEvalOp(PyObject_GetAttrString(py_functor_,
				       "evalOp"));
  }
}

PySundanceCellPredicate::~PySundanceCellPredicate() {
  // Decrement the reference count
  Py_XDECREF(py_functor_);
}


bool PySundanceCellPredicate::operator()(const Point& x) const 
{
  PyObject * arglist;
  PyObject * result;

  switch(x.dim())
    {
    case 1:
      arglist = Py_BuildValue("(d)", x[0]);
      break;
    case 2:
      arglist = Py_BuildValue("(dd)", x[0], x[1]);
      break;
    case 3:
      arglist = Py_BuildValue("(ddd)", x[0], x[1], x[2]);
      break;
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "point dimension = " << x << " not supported");
    }
  result = PyEval_CallObject(py_evalOp_.getFunction(), arglist);
  Py_DECREF(arglist);  // All done with argument list

  if (0 == result) {
    PyErr_Print();
    return false;
  }
  Py_DECREF(result); // All done with returned result object
  return (bool) PyObject_IsTrue(result);
}


PyObject * PySundanceCellPredicate::setEvalOp(PyObject * p_pyObject)
{
  return py_evalOp_.setFunction(p_pyObject);
}

