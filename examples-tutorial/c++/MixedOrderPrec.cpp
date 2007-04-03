// by Robert C. Kirby
// Texas Tech University
// This code is designed to show how to create full-rank preconditioners for
// high order discretizations using low-order discretizations in Sundance
// and Trilinos.


#include "Sundance.hpp"
#include "Thyra_SingleRhsLinearOpBase.hpp"
#include "Thyra_DefaultRealLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_TestingTools.hpp"
#include "EpetraExt_CrsMatrixIn.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_getConst.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})

Expr residual( const Expr & , const Expr & ,
               const QuadratureFamily & );

class MultiOrderPC;
class ResidPlusBC;


// This pure virtual class encapsulates the symbol of an operator:
// its residaul and bc.  Both of these are abstracted over (lists of) test
// and trial functions.  The point of this is I define the symbol of
// an operator in one place, and then put in different orders of test and
// trial functions, so as to reuse code and have a single point of entry
// for the problem I'm solving.
class ResidPlusBC
{
public:
  ResidPlusBC() {}
  virtual ~ResidPlusBC() {}

  virtual Expr residual( const Expr &v , const Expr &u , const QuadratureFamily &q ) const = 0;
  virtual Expr bc( const Expr &v , const Expr &u , const QuadratureFamily &q ) const = 0;
};

// This is an example of filling in the ResidPlusBC class with a (low
// Peclet-number) advection-diffusion operator.  Should be self-explanatory.

class AdvectionDiffusion : public ResidPlusBC
{
public:
  AdvectionDiffusion() {}
  Expr residual(const Expr &v , const Expr &u,
                const QuadratureFamily &q ) const
  {
    double eps = 1.0;
    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    CellFilter interior = new MaximalCellFilter();
    return Integral(interior,eps*((dx*u)*(dx*v)+(dy*u)*(dy*v)) 
                                 + (dx*u)*v - 1.0 * v,q);
  } 

  Expr bc( const Expr &v , const Expr &u , const QuadratureFamily &q) const
  {
    CellFilter points = new DimensionalCellFilter( 0 );
    CellFilter leftPoint = points.subset(new LeftPointTest());
    CellFilter rightPoint = points.subset(new RightPointTest());
    CellFilter topPoint = points.subset( new TopPointTest());
    CellFilter bottomPoint = points.subset( new BottomPointTest() );
    CellFilter allBoundaryPoints = leftPoint+rightPoint+topPoint+bottomPoint;
    return EssentialBC( allBoundaryPoints , v * u , q );
  }
};

// This class turns a particular mesh and a ResidPlusBC object and
// creates a factory for creating linear problem objects given
// the test, trial functions and quadrature rules.
class LinProbMaker
{
public:
  LinProbMaker( const Mesh &mesh ,
                const ResidPlusBC &builder ,
                const VectorType<double> &vecType ) :
    mesh_( mesh ) , 
    builder_( builder ) , 
    vecType_( vecType ) {}

  ~LinProbMaker( ) {}

  Teuchos::RefCountPtr<LinearProblem> 
  getProblem( const Expr &test , const Expr &trial , 
	      const QuadratureFamily &q )
  {
    return Teuchos::rcp( new LinearProblem( mesh_ , 
					    builder_.residual(test,trial,q) , 
					    builder_.bc(test,trial,q) , 
					    test , trial , vecType_ ) );
  }

private:
  const Mesh &mesh_;
  const ResidPlusBC &builder_;
  const VectorType<double> &vecType_;
};


//
// Abstract class takes two BasisFamily object and implements the correct
// interface to use the lower order method to precondition the higher one Ross
// will need to turn this somehow into a Thyra::LinearOpBase.  so we can use
// it as a preconditioner.  I also need to fill in how to apply the transpose
// of the operator.
//
// The PC is an approximation to the inverse of the high order operator, so
// its domain and range are the range and domain of the high order operator,
// respectively These two will be needed to make this a
// Thyra::LinearOpBase<double>
//
class MultiOrderPC : public Thyra::SingleRhsLinearOpBase<double>
{
public:

  // 

  MultiOrderPC(
    DiscreteSpace &lowOrder, DiscreteSpace &highOrder,
    const Teuchos::RefCountPtr<const Thyra::LinearOpBase<double> > &highOrderOp,
    const Teuchos::RefCountPtr<const Thyra::LinearOpBase<double> > &lowOrderPC
    )
    : lowOrderSpace_(lowOrder), highOrderSpace_(highOrder),
      highOrderOp_(highOrderOp) , lowOrderPC_(lowOrderPC)
    {}

  // Overridden from SingeRhsLinearOpBase

  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> >
  range() const
    {
      return highOrderOp_->domain();
    }

  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<double> >
  domain() const {
    return highOrderOp_->range();
  }

  bool opSupported(Thyra::ETransp M_trans) const
    {
      return ( real_trans(M_trans) == Thyra::NOTRANS );
      // Can't handle transpose yet!
    }

  void apply(
    const Thyra::ETransp M_trans,
    const Thyra::VectorBase<double> &x_in,
    Thyra::VectorBase<double> *y_out,
    const double alpha,
    const double beta
    ) const
  {
    typedef Teuchos::ScalarTraits<double> ST;
    TEST_FOR_EXCEPT( 0 == y_out );
    TEST_FOR_EXCEPT( alpha != ST::one() ); // Can't handle yet!
    TEST_FOR_EXCEPT( beta != ST::zero() ); // Can't handle yet!
    if ( real_trans(M_trans) != Thyra::NOTRANS )
      applyTranspose( x_in , *y_out );
    else
      applyNoTranspose( x_in, *y_out );
  }

private:

  DiscreteSpace &lowOrderSpace_ , &highOrderSpace_;

  Teuchos::RefCountPtr<const Thyra::LinearOpBase<double> > highOrderOp_;
  Teuchos::RefCountPtr<const Thyra::LinearOpBase<double> > lowOrderPC_;

  void applyNoTranspose(
    const Thyra::VectorBase<double> &x_in,
    Thyra::VectorBase<double> &y_out
    ) const;

  void applyTranspose(
    const Thyra::VectorBase<double> &x_in,
    Thyra::VectorBase<double> &y_out
    ) const;

};


//
// Implementations
//


void MultiOrderPC::applyNoTranspose(
  const Thyra::VectorBase<double> &x_in,
  Thyra::VectorBase<double> &y_out
  ) const
{

  // Cast the input vector to TSFExtended::Vector<double> to put
  // inside a Sundance DiscreteFunction
  const TSFExtended::Vector<double> x_in_tsf( 
    Teuchos::rcp( const_cast<Thyra::VectorBase<double> * >(&x_in) , false) );

  TSFExtended::Vector<double> y_out_tsf( Teuchos::rcp( &y_out , false ) );

  // Later, these could be moved to member variables to avoid the
  // allocation/deallocation at each application
  Expr dfHigh_ = new DiscreteFunction( highOrderSpace_ , 0.0 , "dfhigh" );
  Expr dfLow_ = new DiscreteFunction( lowOrderSpace_ , 0.0 , "dflow" );

  // cast input vector as a discrete function
  DiscreteFunction::discFunc( dfHigh_ )->setVector( x_in_tsf.copy() );

  // project the input function into the lower order space
  L2Projector proj1( lowOrderSpace_ , dfHigh_ );
  Expr projX = proj1.project();

  // turn the projected input function into a vector and apply the low
  // order PC 
  Vector<double> projXVec = DiscreteFunction::discFunc( projX )->getVector();

  Thyra::apply(
    *lowOrderPC_, NONCONJ_ELE, 
    *(projXVec.ptr()), 
    &*(DiscreteFunction::discFunc(dfLow_)->getVector().ptr())
    );
  
  // put the result back into the high order space.
  // this corresponds to applying the identity to the part orthogonal
  // to the low order space and the preconditioner to the projection
  // into the low order space.  If PC is the low order preconditioner 
  // u is the input function
  // P1 is the projection into the low order space, 
  // P2 is the projection into the high order space
  // this whole operation corresponds to
  // P2 * ( PC * P1 * u + (I - P1) * u )
  L2Projector proj2( highOrderSpace_ , ( dfHigh_ - projX ) + dfLow_ );
  Expr projResult = proj2.project();
  y_out_tsf.acceptCopyOf( DiscreteFunction::discFunc( projResult)->getVector() );

}


void MultiOrderPC::applyTranspose(
  const Thyra::VectorBase<double> &x_in,
  Thyra::VectorBase<double> &y_out
  ) const
{
  TEST_FOR_EXCEPT("MultiOrderPC::applyTranspose(...) not implemented yet!");
}


//
// Main program that creates the physics-based preconditioner and solves the
// lienar system.
//

int main(int argc, char** argv)
{

  using Teuchos::describe;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::getConst;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  typedef ParameterList::PrintOptions PLPrintOptions;
  using Teuchos::sublist;
  using Thyra::inverse;
  using Thyra::prec;
  typedef RefCountPtr<const Thyra::LinearOpBase<double> > LinearOpPtr;
  typedef RefCountPtr<Thyra::VectorBase<double> > VectorPtr;

  bool success = true;

  try {

    Sundance::init(&argc, &argv);

    Teuchos::RefCountPtr<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    *out << "\nA) Setting up the basic operators using Sundance ...\n";

    // needed for Sundance
    VectorType<double> vecType = new EpetraVectorType();
    
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_MEDIUM;

    // get the mesh
    MeshType meshType = new BasicSimplicialMeshType();
    int nx = 100;
    MeshSource mesher = new PartitionedRectangleMesher(0.0,1.0,nx,1,
                                                       0.0,1.0,nx,1,meshType);
    Mesh mesh = mesher.getMesh();

    // create basis objects, quadrature rules and discrete spaces.
    BasisFamily L1 = new Lagrange( 1 );
    BasisFamily L2 = new Lagrange( 2 );

    QuadratureFamily Q1 = new GaussianQuadrature( 2 );
    QuadratureFamily Q2 = new GaussianQuadrature( 4 );

    DiscreteSpace P1Space( mesh , L1 , vecType );
    DiscreteSpace P2Space( mesh , L2 , vecType );

    // these will be the test and trial functions for two
    // LinearProblem objects for the same equation but different
    // orders of discretization.
    Expr u1 = new UnknownFunction( L1 , "u1" );
    Expr u2 = new UnknownFunction( L2 , "u2" );
    Expr v1 = new TestFunction( L1 , "v1" );
    Expr v2 = new TestFunction( L2 , "v2" );

    // create the factory for building linear problems.
    AdvectionDiffusion AdvDiff;
    LinProbMaker LPM( mesh , AdvDiff , vecType );

    // get the Sundance LinearProblem objects for each order.
    Teuchos::RefCountPtr<LinearProblem> LP1 = LPM.getProblem( v1 , u1 , Q1 );
    Teuchos::RefCountPtr<LinearProblem> LP2 = LPM.getProblem( v2 , u2 , Q2 );

    // extract the matrices from each LinearProblem

    LinearOpPtr P1 = LP1->getOperator().ptr();
    LinearOpPtr P2 = LP2->getOperator().ptr();

    // now, set up a preconditioner for P1

    RefCountPtr<ParameterList> paramList = rcp(new ParameterList);
    Teuchos::updateParametersFromXmlFile( "solver.xml", &*paramList );

    const double solveTol =
      getParameter<double>(getConst(*paramList),"Solve Tol");

    Thyra::DefaultRealLinearSolverBuilder P1_linsolve_strategy_builder;
    P1_linsolve_strategy_builder.setParameterList(
      sublist(paramList,"P1 Solver",true) );

    Thyra::DefaultRealLinearSolverBuilder P2_linsolve_strategy_builder;
    P2_linsolve_strategy_builder.setParameterList(
      sublist(paramList,"P2 Solver",true) );

    *out << "\nB) Creating prec(P1) as just an algebraic preconditioner ...\n";

    RefCountPtr<Thyra::PreconditionerFactoryBase<double> > P1_prec_strategy;
    P1_prec_strategy = P1_linsolve_strategy_builder.createPreconditioningStrategy("");
    RefCountPtr<Thyra::PreconditionerBase<double> >
      precP1 = prec<double>(*P1_prec_strategy,P1);
    *out << "\nprecP1 = " << describe(*precP1,verbLevel) << "\n"; 
    LinearOpPtr precP1Op = precP1->getUnspecifiedPrecOp();

    *out << "\nC) Creating the physics-based precondtioiner for P2 ...\n";
    LinearOpPtr precP2Op = rcp(
      new MultiOrderPC( P1Space, P2Space, P2, precP1Op )
      );
    *out << "\nprecP2Op = " << describe(*precP2Op,verbLevel) << "\n"; 

    *out << "\nD) Setting up the solver P2 ...\n";

    RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<double> >
      P2_linsolve_strategy
      = P2_linsolve_strategy_builder.createLinearSolveStrategy("");

    const double usePhysicsPrec = getParameter<bool>(getConst(*paramList),"Use Physics Prec");

    RefCountPtr<Thyra::LinearOpWithSolveBase<double> >
      P2_lows = P2_linsolve_strategy->createOp();

    if (usePhysicsPrec) {
      Thyra::initializePreconditionedOp<double>(
	*P2_linsolve_strategy,
	P2,
	Thyra::unspecifiedPrec(precP2Op),
	&*P2_lows
	);
    }
    else {
      Thyra::initializeOp<double>(
	*P2_linsolve_strategy,
	P2,
        &*P2_lows
	);
    }

    *out << "\nE) Solve P2 for a random RHS ...\n";

    VectorPtr x = createMember(P2->domain());
    VectorPtr b = createMember(P2->range());
    Thyra::randomize(-1.0,+1.0,&*b);
    Thyra::assign(&*x,0.0); // Must give an initial guess!

    Thyra::SolveStatus<double>
      solveStatus = Thyra::solve( *P2_lows, Thyra::NOTRANS, *b, &*x );

    *out << "\nSolve status:\n" << solveStatus;

    *out << "\nSolution ||x|| = " << Thyra::norm(*x) << "\n";
    //
    *out << "\nF) Checking the error in the solution of r=b-P2*x ...\n";
    //
    
    VectorPtr P2x = Thyra::createMember(b->space());
    Thyra::apply( *P2, Thyra::NOTRANS, *x, &*P2x );
    VectorPtr r = Thyra::createMember(b->space());
    Thyra::V_VmV(&*r,*b,*P2x);
    
    double
      P2x_nrm = Thyra::norm(*P2x),
      r_nrm = Thyra::norm(*r),
      b_nrm = Thyra::norm(*b),
      r_nrm_over_b_nrm = r_nrm / b_nrm;

    bool result = r_nrm_over_b_nrm <= solveTol;
    if(!result) success = false;
    
    *out
      << "\n||P2*x|| = " << P2x_nrm << "\n";
    
    *out
      << "\n||P2*x-b||/||b|| = " << r_nrm << "/" << b_nrm
      << " = " << r_nrm_over_b_nrm << " <= " << solveTol
      << " : " << Thyra::passfail(result) << "\n";

    *out << "\nParameter list after use:\n\n";
    paramList->print(*out,PLPrintOptions().indent(2).showTypes(true));
    
  }
  catch(exception &e) {
    Sundance::handleException(e);
    success = false;
  }

  Sundance::finalize();

  return ( success ? 0 : 1 );

}
