#include "SundanceFIATScalarAdapter.hpp"

namespace SundanceStdFwk
{
  FIATScalarAdapter::FIATScalarAdapter( PyObject *pyfamilyclass ) :
    order_( order )
  {
    for (int i=0;i<3;i++) {
      PyObject *arglist = Py_BuildValue( "(dd)" , i+1 , order_ );
      bases_[i] = PyEval_CallObject( pyfamilyclass , arglist );
      Py_DECREF( arglist );
    }
  }
  
  FIATScalarAdapter::~FIATScalarAdapter()
  {
    for (int i=0;i<3;i++) {
      Py_DECREF( bases_[i] );
    }
  }
  
  void FIATScalarAdapter::getLocalDOFs(const CellType& cellType,
				       Array<Array<Array<int> > >& dofs ) const
  {
    if (PointCell == cellType) {
    }
    else {
      const int celldim = dimension( cellType );
      PyObject *basis = bases_[ celldim - 1 ];
      
      PyObject *dual_basis = PyObject_CallMethod( basis , "dual_basis" ,
						  NULL );
      
      dofs.resize( celldim );
      
      
      for (int i=0;i<celldim;i++) {
	PyObject *nodes_per_dim = 
	  PyObject_CallMethod( dual_basis , "getNodeIds" , 
			       "i" , i ); 
	
	int num_facets = PyObject_Length( nodes_per_dim );
	
	dofs[i].resize( num_facets );
	
	for (int j=0;j<num_facets;j++) {
	  // convert j to a Python integer
	  PyObject *pyj = PyInt_FromLong( (long) j );
	  
	  // Extract the Python list of integers associated
	  // with facet j of dimension i
	  PyObject *nodes_per_facet =
	    PyObject_GetItem( nodes_per_dim , pyj );
	  
	  int num_nodes_cur = PyObject_Length( nodes_per_facet );
	  
	  
	  // resize the innermost array and copy the data from
	  // the Python list to the Techos array.
	  dofs[i][j].resize( num_nodes_cur );
	  for (int k=0;k<num_nodes_cur;k++) {
	    PyObject *pyk = PyInt_FromLong( (long) k );
	    PyOjbect *pynodecur = PyObject_GetItem( nodes_per_facet , pyk );
	    dofs[i][j][k] = PyInt_AsLong( pynodecur );
	    Py_DECREF( pynodecur );
	    Py_DECREF( pyk );
	  }
	  
	  Py_DECREF( pyj );
	  Py_DECREF( nodes_per_facet );
	}
	
	Py_DECREF( nodes_per_dim );
	Py_DECREF( dual_basis );
      }
      
    }
  }
  

  // sum over spatial dimensions up to and including spatialDim
  // for the given cellType
  int FIATScalarAdapter::nNodes( int spatialDim,
				 const CellType& cellType ) const 
  {
    if (PointCell == cellType) {
    }
    else {
      int celldim = dimension( cellType );
      // throw a Sundance exception if spatialDim > cellDim
      PyObject *basis = bases[celldim-1];
      PyObject *dual_basis = PyObject_CallMethod( basis , 
						  "dual_basis" ,
						  NULL );
      int nnod = 0;
      for (int i=0;i<=spatialDim;i++) {
	PyObject *nodes_per_dim =
	  PyObject_CallMethod( dual_basis , "getNodeIDs" ,
			       "i" , i );

	PyObject *nodes_on_facet0 =
	  PyObject_GetItem( nodes_per_dim , 0 );
	
	nnod += PyObject_Length( nodes_on_facet0 );

	Py_DECREF( nodes_on_facet0 );
	Py_DECREF( nodes_per_dim );
      }

      Py_DECREF( dual_basis );
    }

    return nnod;
  }

  void FIATScalarAdapter::refEval( int spatialDim ,
				   const CellType& cellType ,
				   const Array<Point>& pts,
				   const MultiIndex& deriv ,
				   Array<Array<double> >& result) const
  {
    if (PointCell == cellType) {
      // what about derivatives?
      result = tuple(tuple(1.0));
    }
    else {
      // Get list of points to put into FIAT
      int cellDim = dimension( cellType );
      PyObject *py_list_of_points = PyList_New(pts.size());

      for (int i=0;i<pts.size();i++) {
	// Create a Python tuple for the point, converting from
	// Sundance (0,1)-based coordinates to FIAT (-1,1)-based coordinates
	PyObject *py_pt_cur = PyTuple_New( spatialDim );
	for (int j=0;j<cellDim;j++) {
	  double coord_cur = 2.0 * ( pts[i][j] - 0.5 );
	  PyObject *py_coord = PyFloat_FromDouble( coord_cur );
	  PyTuple_SetItem( py_pt_cur , j , py_coord );
	}
	for (int j=cellDim;j<spatialDim;j++) {
	  PyObject *py_coord = PyFloat_FromDouble( -1.0 );
	  PyTuple_SetItem( py_pt_cur , j , py_coord );
	}
	PyObject_SetItem( py_list_of_points , i , py_pt_cur );
      }

      // Extract the function space from the basis
      PyObject *py_basis = bases_[spatialDim-1];
      PyObject *py_function_space = 
	PyObject_CallMethod( py_basis ,
			     "function_space" , NULL );

      // convert the multiindex from Sundance into a Python tuple
      PyObject *py_alpha = PyTuple_New( spatialDim );
      for (int i=0;i<spatialDim;i++) {
	PyObject *py_alpha_i = PyInt_FromLong( (long) alpha[i] );
	PyObject_SetItem( py_alpha , i , py_alpha_i );
      }

      // Call multi_deriv_all on the function space; requires
      PyObject *py_deriv_space =
	PyObject_CallMethod( py_function_space , 
			     "multi_deriv_all" , 
			     "O" , py_alpha );
      
      
			     

      // constructing an appropriate Python tuple out of "deriv"


      //
      // Call trace_tabulate on the polynomial set of derivatives
      //
      // Iterate over the resulting Numeric array and move
      // everything into "result"
      //
      // Decrement all references, including a loop over the list of points
      // with a loop over the doubles inside.
      Py_DECREF( py_deriv_space );
      Py_DECREF( py_function_space );


    }
  }

  int dim() const
  {
  }




  
}
