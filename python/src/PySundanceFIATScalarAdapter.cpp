#include "SundanceFIATScalarAdapter.hpp"

#ifdef HAVE_PY_FIAT

static int line_sdvert_to_fvert[] = {0,1};
static int line_sdline_to_fline[] = {0};
static int *line_sd_to_fiat[] = {line_sdvert_to_fvert,line_sdline_to_fline};
static int tri_sdvert_to_fvert[] = {0,1,2};
static int tri_sdline_to_fline[] = {2,0,1};
static int tri_sdtri_to_ftri[] = {0};
static int *tri_sd_to_fiat[] = {tri_sdvert_to_vert,
				tri_sdline_to_fline,
				tri_sdtri_to_ftri};
static int tet_sdvert_to_fvert[] = {0,1,2,3};
static int tet_sdline_to_fline[] = {2,0,1,3,4,5};
static int tet_sdtri_to_ftri[] = {0,1,2,3};
static int tet_sdtet_to_ftet[] = {0};

static int *tet_sd_to_fiat[] = {tet_sdvert_to_fvert,
				tet_sdline_to_fline,
				tet_sdtri_to_ftri,
				tet_sdtet_to_ftet};

static int **sd_to_fiat[] = {NULL,
			     line_sd_to_fiat,
			     tri_sd_to_fiat,
			     tet_sd_to_fiat};

static CellType sdim_to_cellType[] = 
  {PointCell,LineCell,TriangleCell,TetCell};

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
      // keep a stack of every python object we create.  This makes
      // DECREF'ing easy to do
      std::stack<PyObject *> to_decref;

      // Get list of points to put into FIAT
      int cellDim = dimension( cellType );
      PyObject *py_list_of_points = PyList_New(pts.size());
      to_decref.push( py_list_of_points );

      for (int i=0;i<pts.size();i++) {
	// Create a Python tuple for the point, converting from
	// Sundance (0,1)-based coordinates to FIAT (-1,1)-based coordinates
	PyObject *py_pt_cur = PyTuple_New( spatialDim );
	to_decref.push( py_pt_cur );
	for (int j=0;j<cellDim;j++) {
	  double coord_cur = 2.0 * ( pts[i][j] - 0.5 );
	  PyObject *py_coord = PyFloat_FromDouble( coord_cur );
	  PyTuple_SetItem( py_pt_cur , j , py_coord );
	  to_decref.push( py_coord );
	}
	for (int j=cellDim;j<spatialDim;j++) {
	  PyObject *py_coord = PyFloat_FromDouble( -1.0 );
	  PyTuple_SetItem( py_pt_cur , j , py_coord );
	  to_decref.push( py_coord );
	}
	PyObject_SetItem( py_list_of_points , i , py_pt_cur );
      }

      // Extract the function space from the basis
      PyObject *py_basis = bases_[spatialDim-1];
      PyObject *py_function_space = 
	PyObject_CallMethod( py_basis ,
			     "function_space" , NULL );
      to_decref.push( py_function_space );

      // convert the multiindex from Sundance into a Python tuple
      PyObject *py_alpha = PyTuple_New( spatialDim );
      to_decref.push( py_alpha );
      for (int i=0;i<spatialDim;i++) {
	PyObject *py_alpha_i = PyInt_FromLong( (long) alpha[i] );
	to_decref.push( py_alpha_i );
	PyObject_SetItem( py_alpha , i , py_alpha_i );

      }

      // Call multi_deriv_all on the function space; requires
      PyObject *py_deriv_space =
	PyObject_CallMethod( py_function_space , 
			     "multi_deriv_all" , 
			     "O" , py_alpha );
      to_decref.push( py_deriv_space );

      // This should be a Numeric.array object.  I can
      // go to the low-level Numeric API and grab the C pointer to
      // the data directly if we have a speed problem.
      // 
      // This function should give a 2d array that is num_pts by
      // num_bf
      PyObject *py_tabulation =
	PyObject_CallMethod( py_deriv_space , 
			     "tabulate" ,
			     "O" , py_list_of_points );
      to_decref.push( py_tabulation );

      PyObject *py_zero = PyInt_FromLong( (long) 0 );
      to_decref.push( py_zero );

      PyObject *py_tabulation_first_row =
	PyObject_GetItem( py_tabulation , py_zero );
      to_decref.push( py_tabulation_first_row );

      int num_bf = PyObject_Length( py_tabulation_first_row );

      result.resize( pts.size() );
      for (int i=0;i<pts.size();i++) {
	result.resize( num_bf );
      }

      CellType ct = sdim_to_cellType[ sdatialDim ];
      int cd = dimension[ct];

      
      Array<Array<Array<int> > > dofs;
      getLocalDOFs( ct , dofs );


      int cur = 0;
      int **sd_to_fiat_spd = sd_to_fiat[spatialDim];
      for (int d=0;d<=pd;d++) {
	int *sd_to_fiat_spd_d = sd_to_fiat_spd[d];
	for (int e=0;e<=dofs[e].size();e++) {
	  int fiat_e = sd_to_fiat_spd_d[e];
	  for (int n=0;n<(int)dofs[e][d].size();n++) {
	    for (int p=0;p<(int)pts.length();p++) {
	      PyObject *py_ij_tuple = 
		Py_BuildObject( "(ii)" , dofs[d][fiat_e][n] , p );
	      to_decref( py_ij_tuple );
	      PyObject *py_tab_cur_ij = 
		PyObject_GetItem( py_tabulation , py_ij_tuple );
	      to_decref( py_tab_cur_ij );
	      result[p][cur] = PyFloat_AsDouble( py_tab_cur_ij );
	    }
	    cur++;
	  }
	}
      }

      while (!to_decref.empty()) {
	PyObject *cur = to_decref.top();
	Py_DECREF( cur );
	to_decref.pop();
      }

    }
  }
  
}

#endif
