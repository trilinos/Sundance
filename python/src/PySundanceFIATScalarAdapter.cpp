#include "PySundanceFIATScalarAdapter.hpp"
#include <stack>
#include <iostream>
using namespace std;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;


static int line_sdvert_to_fvert[] = {0,1};
static int line_sdline_to_fline[] = {0};
static int *line_sd_to_fiat[] = {line_sdvert_to_fvert,line_sdline_to_fline};
static int tri_sdvert_to_fvert[] = {0,1,2};
static int tri_sdline_to_fline[] = {2,0,1};
static int tri_sdtri_to_ftri[] = {0};
static int *tri_sd_to_fiat[] = {tri_sdvert_to_fvert,
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

#define HAVE_PY_FIAT
#ifdef HAVE_PY_FIAT

namespace SundanceStdFwk
{
	FIATScalarAdapter::FIATScalarAdapter( PyObject *pyfamilyclass ,
										  int order ) :
		order_( order )
	{
    	stack<PyObject *> to_decref;

		cout << PyCallable_Check( pyfamilyclass ) << endl;

	    /* instantiate basis for each shape */
		bases_.resize( 3 );
		for (int i=0;i<3;i++) {
			PyObject *arglist = Py_BuildValue( "(ii)" , i+1 , order_ );
			if (!arglist) cout << "barf" << endl;
			bases_[i] = PyObject_CallObject( pyfamilyclass , arglist  );
			if (!bases_[i]) cout << "puke" << endl;
			to_decref.push( arglist );
	    }

    	/* pretabulate all the dof in the constructor */
    	/* one dof array for each cellType I instantiate */
		dof_.resize( 3 );
		for (int cd=1;cd<=3;cd++) {
			dof_[cd-1].resize( cd + 1 );
			PyObject *basis = bases_[cd-1];
			PyObject *dual_basis = PyObject_CallMethod( basis , "dual_basis" ,
						 						    	NULL );
			to_decref.push( dual_basis );
			for (int i=0;i<=cd;i++) {
		  		PyObject *nodes_per_dim =
		  		PyObject_CallMethod( dual_basis , "getNodeIDs" ,
					  	         	 "i" , i );
				to_decref.push( nodes_per_dim );
				int num_facets = PyObject_Length( nodes_per_dim );
				dof_[cd-1][i].resize( num_facets );
				for (int j=0;j<num_facets;j++) {
					cout << "facet number " << j << endl;
		  			PyObject *pyj = PyInt_FromLong( (long) j );
		  			to_decref.push( pyj );
		  			PyObject *nodes_per_facet =
		    		PyObject_GetItem( nodes_per_dim , pyj );
		  			to_decref.push( nodes_per_facet );
		  			int num_nodes_this_facet = PyObject_Length( nodes_per_facet );
	  
		  			dof_[cd-1][i][j].resize( num_nodes_this_facet );
		  			for (int k=0;k<num_nodes_this_facet;k++) {
		    			PyObject *pyk = PyInt_FromLong( (long) k );
		    			to_decref.push( pyk );
		    			PyObject *pynodecur = PyObject_GetItem( nodes_per_facet , pyk );
		    			to_decref.push( pynodecur );
		    			dof_[cd-1][i][j][k] = (int) PyInt_AsLong( pynodecur );
					}
				}
			}
		}

    	while (!to_decref.empty()) {
      		PyObject *foo = to_decref.top();
      		Py_DECREF( foo );
      		to_decref.pop();
    	}	

	}
  
	FIATScalarAdapter::~FIATScalarAdapter()
	{
		for (int i=0;i<3;i++) {
			Py_DECREF( bases_[i] );
    	}
	}
  
  	/* loops over pretabulated data structures and copies into dofs */
	void FIATScalarAdapter::getLocalDOFs(const CellType& cellType,
										 Array<Array<Array<int> > >& dofs ) const
	{
    	if (PointCell == cellType) {
			dofs.resize(1);
			dofs[0] = tuple(tuple(0));
			return;
	    }
    	else {
			int cd = dimension( cellType );
			Array<Array<Array<int> > >& dofs_cur = dofs_[cd-1];
			dofs.resize( dofs_cur.size() );
			for (unsigned i=0;i<dofs_cur.size();i++) {
				dofs[i].resize( dofs_cur[i].size() );
				for (unsigned j=0;j<dofs_cur[i].size();j++) {
					dofs[i][j].resize( dofs_cur[i][j].size() );
					for (unsigned k=0;k<dofs_cur[i][j].size();k++) {
						dofs[i][j][k] = dofs_cur[i][j][k];
					}
				}
			}
			return;
  		}
	}

	// sum over spatial dimensions up to and including spatialDim
	// for the given cellType
	int FIATScalarAdapter::nNodes( int spatialDim,
								   const CellType& cellType ) const 
	{
		if (PointCell == cellType) {
			return 1;
		}
		else {
			int celldim = dimension( cellType );
			int nn = 0;
			Array<Array<Array<int> > >& dofs_cur = dof_[celldim-1];
			for (int i=0;i<=spatialDim;i++) {
				nn += dofs_cur[i][0].size();
			}

			return nn;
		}
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
      stack<PyObject *> to_decref;

      // Get list of points to put into FIAT
      int cellDim = dimension( cellType );
      PyObject *py_list_of_points = PyList_New(pts.size());
      to_decref.push( py_list_of_points );

      for (int i=0;i<(int)pts.size();i++) {
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
	PyObject *py_i = PyInt_FromLong( (long) i );
	to_decref.push( py_i );
	PyObject_SetItem( py_list_of_points , py_i , py_pt_cur );
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
	PyObject *py_alpha_i = PyInt_FromLong( (long) deriv[i] );
	to_decref.push( py_alpha_i );
	PyObject *py_i = PyInt_FromLong( (long) i );
	to_decref.push( py_i );
	PyObject_SetItem( py_alpha , py_i , py_alpha_i );

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
      for (int i=0;i<(int)pts.size();i++) {
	result.resize( num_bf );
      }

      CellType ct = sdim_to_cellType[ spatialDim ];
      int cd = dimension(ct);

      
      Array<Array<Array<int> > > dofs;
      getLocalDOFs( ct , dofs );


      int cur = 0;
      int **sd_to_fiat_spd = sd_to_fiat[spatialDim];
      for (int d=0;d<=cd;d++) {
	int *sd_to_fiat_spd_d = sd_to_fiat_spd[d];
	for (int e=0;e<=(int)dofs[e].size();e++) {
	  int fiat_e = sd_to_fiat_spd_d[e];
	  for (int n=0;n<(int)dofs[e][d].size();n++) {
	    for (int p=0;p<(int)pts.length();p++) {
	      PyObject *py_ij_tuple = 
		Py_BuildValue( "(ii)" , dofs[d][fiat_e][n] , p );
	      to_decref.push( py_ij_tuple );
	      PyObject *py_tab_cur_ij = 
		PyObject_GetItem( py_tabulation , py_ij_tuple );
	      to_decref.push( py_tab_cur_ij );
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
