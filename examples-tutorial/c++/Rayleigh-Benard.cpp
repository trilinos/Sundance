#include "Sundance.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceEvaluator.hpp"

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})


Expr flowEquation( Expr flow , 
		   Expr lagFlow ,
		   Expr varFlow , Expr temp ,
		   Expr rayleigh, Expr inv_prandtl ,
		   QuadratureFamily quad )
{
  CellFilter interior = new MaximalCellFilter();
  /* Create differential operators  */
  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);
  Expr grad = List(dx, dy);

  Expr lagUx = lagFlow[0]; Expr lagUy = lagFlow[1];
  Expr lagU = List( lagUx , lagUy ); 
  Expr ux = flow[0]; Expr uy = flow[1]; 
  Expr u = List( ux , uy );
  Expr vx = varFlow[0]; Expr vy = varFlow[1];
  Expr p = flow[2]; Expr q = varFlow[2];
  Expr temp0 = temp;
  return Integral(interior,
		  (grad*vx)*(grad*ux) + (grad*vy)*(grad*uy)
		  + vx*(lagU*grad)*ux + vy*(lagU*grad)*uy 
		  - p*(dx*vx+dy*vy) - q*(dx*ux+dy*uy) 
		  - temp0*rayleigh*inv_prandtl*vy,quad);
}

Expr flowBoundaryCondition( Expr flow , Expr varFlow ,
			    QuadratureFamily quad )
{
  CellFilter bdry = new BoundaryCellFilter();
  Expr ux = flow[0]; Expr uy = flow[1];
  Expr vx = varFlow[0]; Expr vy = varFlow[1];
  return EssentialBC( bdry , ux*vx + uy*vy , quad );
}


Expr tempEquation( Expr temp , Expr varTemp , Expr flow ,
		   Expr inv_prandtl , 
		   QuadratureFamily quad )
{
  CellFilter interior = new MaximalCellFilter();
  Expr dx = new Derivative(0); Expr dy = new Derivative(1);
  Expr grad = List(dx, dy);  
  
  return Integral( interior ,
		   inv_prandtl * (grad*temp)*(grad*varTemp)
		   + (flow[0]*(dx*temp)+flow[1]*(dy*temp))*varTemp ,
		   quad );
}

Expr tempBoundaryCondition( Expr temp , Expr varTemp , QuadratureFamily quad )
{
  CellFilter edges = new DimensionalCellFilter(1);
  CellFilter top = edges.subset( new TopPointTest() );
  CellFilter bottom = edges.subset( new BottomPointTest() );  

  return EssentialBC( top , temp*varTemp , quad )
    + EssentialBC( bottom , (temp-1.0)*varTemp , quad );
}

using SundanceCore::List;
/** 
 * \example iterative.cpp
 * 
 * Solves the coupled NSE/thermal convection 
 */



int main(int argc, char** argv)
{
  try {
    Sundance::init(&argc, &argv);
    int np = MPIComm::world().getNProc();
    
    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    ParameterXMLFileReader reader( "bouss.xml" );
    ParameterList params = reader.getParameters();
    
    const int nx = params.get( "nx" , 1 );

    MeshType meshType = new BasicSimplicialMeshType();
    
    MeshSource mesher 
      = new PartitionedRectangleMesher( 0.0 , 1.0 , nx , np ,
					0.0 , 1.0 , nx , 1 ,
					meshType );

    Mesh mesh = mesher.getMesh();
    CellFilter interior = new MaximalCellFilter();
    
    /* Create unknown and test functions using Taylor-Hood for fluids
     * and linears for temperature */
    BasisFamily L2 = new Lagrange(2);
    BasisFamily L1 = new Lagrange(1);

    /* fluid test and trial functions */
    Expr ux = new UnknownFunction(L2, "u_x");
    Expr uy = new UnknownFunction(L2, "u_y");
    Expr u = List(ux, uy);
    Expr p = new UnknownFunction(L1, "p");
    Expr flow = List( ux , uy , p );
    Expr vx = new TestFunction(L2, "v_x");
    Expr vy = new TestFunction(L2, "v_y");
    Expr v = List( vx , vy );
    Expr q = new TestFunction(L1, "q");
    Expr varFlow = List( vx , vy , q );

    /* temperature test and trial functions */
    Expr T = new UnknownFunction(L1,"T");
    Expr w = new TestFunction(L1,"w");

  
    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad2 = new GaussianQuadrature(2);
    
    /* parameter expression for the Raleigh and Prandtl number */
    Expr rayleigh =
      new SundanceCore::Parameter( params.get( "Ra" , 0.0 ) );
    Expr inv_prandtl
      = new SundanceCore::Parameter( params.get( "InvPr" , 0.0 ) );


    /* create Discrete Spaces for the flow and temperature variables */
    DiscreteSpace flowSpace( mesh , List( L2 , L2 , L1 ) , vecType );
    DiscreteSpace tempSpace( mesh , L1 , vecType );

    /* create DiscreteFunction, inititialized to zero, to store the
       initial guess and solution for the nonlinear systems */
    Expr flow0 = new DiscreteFunction( flowSpace , 0.0 , "u0" );
    Expr temp0 = new DiscreteFunction( tempSpace , 0.0 , "T" );

    /* These DiscreteFunction objects hold the previous iterates for
       comparison */ 
    Expr flowold = new DiscreteFunction( flowSpace , 0.0 , "uold" );
    Expr tempold = new DiscreteFunction( tempSpace , 0.0 , "told" );
    
    /* call the functions above to set up the split system.
       The flow equation uses the previous iterate of the temperature */
    Expr flowEqn = flowEquation( flow , flow0 ,  
				 varFlow , temp0 ,
				 rayleigh , inv_prandtl , quad2 );
    
    Expr flowBC = flowBoundaryCondition( flow , varFlow , quad2 );

    /* equation and BC for temp update */
    Expr tempEqn = tempEquation( T , w , flow0 , inv_prandtl , quad2 );

    Expr tempBC = tempBoundaryCondition( T , w , quad2 );


    LinearProblem flowProb( mesh , flowEqn , flowBC , varFlow , flow ,
			    vecType );

    LinearProblem tempProb( mesh , tempEqn , tempBC , w , T , vecType );


    /* set up nonlinear solver.  We will use the same parameters for
       both the flow and temperature equation */
    const string nonlinearReaderFileName =
      params.get( "NonlinearSolverFile" , "" );
    ParameterXMLFileReader nonlinearReader( nonlinearReaderFileName );
    ParameterList noxParams = nonlinearReader.getParameters();

    /* set up linear solver */
    const string linearReaderFileName =
      params.get( "LinearSolverFile" , "" );
    ParameterXMLFileReader linearReader( linearReaderFileName );
    ParameterList linParams = linearReader.getParameters();
    LinearSolver<double> linSolver 
      = LinearSolverBuilder::createSolver( linParams );

    const int maxiter = params.get( "MaxIters" , 1 );

    int i=0;
    while( i<maxiter ) {
      // back up vector data for comparison later
      DiscreteFunction::discFunc( tempold )->setVector( DiscreteFunction::discFunc( temp0 )->getVector().copy() );
      DiscreteFunction::discFunc( flowold )->setVector( DiscreteFunction::discFunc( flow0 )->getVector().copy() );

      Expr flownew = flowProb.solve( linSolver );

      DiscreteFunction::discFunc( flow0 )->setVector( DiscreteFunction::discFunc( flownew )->getVector().copy() );

      Expr tempnew = tempProb.solve( linSolver );

      DiscreteFunction::discFunc( temp0 )->setVector( DiscreteFunction::discFunc( tempnew )->getVector().copy() );

      // compare new and old flow and temperature
      Expr udiff_sq_expr =
	Integral( interior ,
		  pow( flowold[0] - flow0[0] , 2)
		  + pow( flowold[1] - flow0[1] , 2 ) ,
		  quad2 );
      Expr tempdiff_sq_expr =
	Integral( interior ,
		  pow( tempold - temp0 , 2 ) ,
		  quad2);
      double udiff = sqrt( evaluateIntegral( mesh , udiff_sq_expr ) );
      double tdiff = sqrt( evaluateIntegral( mesh , tempdiff_sq_expr ) );
      
      cout << i << endl;
      cout << "vel diff" << udiff  << endl;
      cout << "temp diff" << tdiff << endl;

      // output this iterature to watch the evolution of the solver
      FieldWriter writer = new VTKWriter( "bouss" + Teuchos::toString(i) );
      writer.addMesh( mesh );
      writer.addField( "ux" , new ExprFieldWrapper( flow0[0] ) );
      writer.addField( "uy" , new ExprFieldWrapper( flow0[1] ) );
      writer.addField( "p" , new ExprFieldWrapper( flow0[2] ) );
      writer.addField( "T" , new ExprFieldWrapper( temp0 ) );
      writer.write();

      i++;
    }

    /* now, we are going to do full Newton starting from the last
       iterate of the split method.  First, we create the full
       discrete space and project the iterate onto it */
    DiscreteSpace fullSpace( mesh , List( L2 , L2 , L1 , L1 ) , vecType );    
    L2Projector fullProj( fullSpace ,
			  List( flow0[0] , flow0[1] , flow0[2] , temp0 ) );
    Expr full0 = fullProj.project();

    /* now we will set up the fully coupled equations.  Notice this
       puts an UnknownFunction in for the temperature in the flow
       equations and the flow variables for the temperature equations,
       whereas before we put in a DiscreteFunction. */
    Expr fullEqn = flowEquation( flow , flow , varFlow , T ,
				 rayleigh , inv_prandtl , quad2 )
      + tempEquation( T , w , flow , inv_prandtl , quad2 );
    Expr fullBC = flowBoundaryCondition( flow , varFlow , quad2 )
      + tempBoundaryCondition( T , w , quad2 );

    NonlinearOperator<double> fullOperator
      = new NonlinearProblem( mesh , fullEqn , fullBC ,
			      List( vx , vy , q , w ) ,
			      List( ux , uy , p , T ) ,
			      full0 , vecType );

    NOXSolver fullSolver( noxParams , fullOperator );
    NOX::StatusTest::StatusType fullStatus = fullSolver.solve();
    
    
    FieldWriter writer = new VTKWriter( "bouss-newt" );
    writer.addMesh( mesh );
    writer.addField( "ux" , new ExprFieldWrapper( full0[0] ) );
    writer.addField( "uy" , new ExprFieldWrapper( full0[1] ) );
    writer.addField( "p" , new ExprFieldWrapper( full0[2] ) );
    writer.addField( "T" , new ExprFieldWrapper( full0[3] ) );
    writer.write();

    
    
  }
  catch(exception& e) {
      cerr << e.what() << endl;
    }
  Sundance::finalize();
}
