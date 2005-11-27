// -*- c++ -*-



%{

#ifdef HAVE_MPI
#include "mpi.h"
  PyObject* Init_Argv(PyObject *args) 
    {  
      int i, error, myid, size;
      int argc = 0;  
      char **argv;   

      /* Reconstruct C-commandline */     
      /*                           */ 
      argc = PyList_Size(args); //Number of commandline arguments
      argv = (char**) malloc((argc+1)*sizeof(char*)); 
  
      for (i=0; i<argc; i++)  
        argv[i] = PyString_AsString( PyList_GetItem(args, i) );
    
      argv[i] = NULL; //Lam 7.0 requires last arg to be NULL  
  
      error = MPI_Init(&argc, &argv); 
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
      MPI_Comm_size(MPI_COMM_WORLD, &size);    

      if (error != 0) {
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);    
        PyErr_SetString(PyExc_RuntimeError, "error");   
        return NULL;
      }  

      return Py_BuildValue("");  
    } 
  
  PyObject* Finalize() {  
    int error, myid;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);  
  
    error = MPI_Finalize();
    if (error != 0) {
      PyErr_SetString(PyExc_RuntimeError, "error");    //raise ValueError, errmsg
      return NULL;
    }  
    
    return Py_BuildValue("");
  }
#else
  PyObject* Init_Argv(PyObject *args) 
    {
      return Py_BuildValue("");
    }

  PyObject* Finalize()
    {
      return Py_BuildValue("");
    }
#endif
 
  %}

// MPI stuff
PyObject* Init_Argv(PyObject *args);
PyObject* Finalize();

// Python code.  Here we set the __version__ string, call MPI_Init()
// and arrange for MPI_Finalize to be called (if appropriate), and
// declare classes that inherit both from Epetra objects and
// UserArrays, to give users additional functionality.
%pythoncode %{

# Call MPI_Init if appropriate
import sys
Init_Argv(sys.argv)
del sys

# Arrange for MPI_Finalize to be called at exit, if appropriate
import atexit
atexit.register(Finalize)


    %}
