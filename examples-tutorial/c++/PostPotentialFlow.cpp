/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"

/** 
 * \example PostPotentialFlow.cpp
 * Solves the Laplace equation for potential flow past an elliptical 
 * post in a wind tunnel.  This example will serve as a master example
 * and outline all the major components commen to nearlt all Sundance 
 * simulators.
 * 
 */


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);

      /** 
       * \subsubsection{Boilerplate}
       * A dull but essential first step is to show 
       * the boilerplate C++ common to nearly every Sundance code:
       * \begin{verbatim}
       * #include "Sundance.hpp"

       * int main(int argc, void** argv)
       * {
       * try
       *     {
       *       Sundance::init(argc, argv);
       *
       *    }
       *  catch(exception& e)
       *    {
       *       Sundance::handleException(__FILE__, e);
       *     }
       *  Sundance::finalize();
       * }
       * \end{verbatim}
       * These lines control initialization and result gathering for
       * profiling timers, initializing and
       * finalizing MPI if MPI is being used, and other administrative
       * tasks. 
       * We will use Epetra as our linear algebra representation, which is
       * specified by selecting the corresponding \verb+VectorType+ subtype,
       */
      VectorType<double> vecType = new EpetraVectorType();

      /** 
       * \subsubsection{Getting the mesh}
       *
       * Sundance uses a {\tt Mesh} object to represent 
       * a tesselation of the problem domain.
       * There are many ways of getting a mesh, all abstracted with the 
       * \verb+MeshSource+ interface. Sundance is designed to work with
       * different mesh underlying implementations, the choice of which is
       * done by specifying a \verb+MeshType+ object. In this example,
       * we use the \verb+BasicSimplicialMeshType+ which is a lightweight parallel
       * simplicial mesh. The mesh was created using 
       * Cubit~(\verb+http://cubit.sandia.gov+) and
       * saved to a NetCDF file. 
       */

      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusNetCDFMeshReader("post.ncdf", meshType);
      Mesh mesh = mesher.getMesh();
      /**
       * If you know a little C++ -- just enough to be dangerous -- you might 
       * think it odd
       * that the result of the {\tt new} operator, which returns a pointer, 
       * is being assigned
       * to a {\tt MeshType} object which is -- apparently -- not a pointer. 
       * That's not
       * a typo: the {\tt MeshType} object is a {\bf handle} class that stores and
       * manages the pointer to the {\tt BasicSimplicialMeshType} object. 
       * Handle classes are used
       * throughout user-level Sundance code, and among other things 
       * relieve you of the need to worry about memory management. 
       *
       * \subsubsection{Defining the domains of integration}
       *
       * We've already read a mesh. We need a way to specify {\it where} on the mesh
       * equations or boundary conditions are to be applied. Sundance 
       * uses a {\tt CellFilter} object
       * to represent subregions of a geometric 
       * domain. A {\tt CellFilter} can be any collection of mesh cells,
       * for example a block of maximal cells, a set of boundary edges, or a set of points.
       *
       * We need do nothing for the wall BCs because they are natural and drop
       * out of the weak form. We need only identify the interior, the
       * inlet, and the outlet.
       * We first create a cell filter object for the entire domain,
       */

      CellFilter interior = new MaximalCellFilter();
      /**
       * which will identify all maximal cells. We next create a cell
       * filter to identify the boundary cells,
       */
      CellFilter boundary = new BoundaryCellFilter();

      /**
       * and then we find the subsets of the boundary corresponding to the
       * inlet and outlet. These have been labeled by Cubit as ``1'' and ``2,''
       * respectively, so we can grab the subsets having those labels:
       */

      CellFilter in = boundary.labeledSubset(1);
      CellFilter out = boundary.labeledSubset(2);
      /**
       * \subsubsection{Defining unknown and test functions}
       *
       * We'll use first order piecewise Lagrange interpolation to represent our unknown
       * solution $\phi$. With a Galerkin method we use a test function ${\hat\phi}$
       * defined using the same basis as the unknown. Expressions representing 
       * the test and unknown functions are defined easily:
       */

      Expr phi = new UnknownFunction(new Lagrange(1), "u");
      Expr phiHat = new TestFunction(new Lagrange(1), "v");
      /**
       * \subsubsection{Creating the gradient operator}
       *
       * The gradient operator is formed by making a {\tt List} containing the partial
       * differentiation operators in the $x$, $y$, and $z$ directions. 
       */

      Expr x = new CoordExpr(0);
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);

      /**
       * Notice that we always number directions (and all other indices)
       * starting from zero as in C, C++, and Java,
       * as opposed to from one as in Fortran and Matlab.
       * The gradient thus defined 
       * is treated as a vector with respect to the overloaded multiplication 
       * operator used to apply the gradient, so that an operation such as {\tt grad*u}
       * expands correctly to {\tt \{dx*u, dy*u, dz*u\}}.
       *
       * \subsubsection{Writing the weak form}
       * 
       * We will need to define expressions for $L$ and $\phi_{out}$. Such constant
       * expressions are defined easily as
       */

      QuadratureFamily quad2 = new GaussianQuadrature(2);

      double L = 1.0;
      /**
       * The weak form for the interior, walls, and outlet is
       * \begin{equation}
       * \int_\text{interior} \nabla{\hat\phi}\cdot\nabla\phi 
       * -\int_\text{outlet} \frac{1}{L}{\hat\phi}\left(\phi_{out} - \phi\right) = 0.
       * \end{equation}
       * All of these integrations can be done exactly so there is no 
       * need to specify a quadrature rule, however, there is no harm in doing so;
       * if an exactly integrable term is detected any quadrature rule will be ignored. 
       * For generality in the event we later wish to change one of the coefficients
       * to something non-constant, we explicitly specify a quadrature rule. We'll
       * use second-order Gaussian quadrature. The
       * weak form with a quadrature specification is written in Sundance as
       */ 

      Expr eqn = Integral(interior, (grad*phiHat)*(grad*phi), quad2)
        + Integral(in, phiHat*(x-phi)/L, quad2) ;

	/**
	* \subsubsection{Writing the essential BCs}
	* 
	* The weak form above contains the physics in the body of the domain 
	* plus the Neumann BCs on the walls and
	* the Robin boundary conditions on the outlet. We still need to apply the 
	* Dirichlet boundary condition 
	* on the inlet, which we do with an {\tt EssentialBC} object
	*/

      Expr bc = EssentialBC(out, phiHat*phi/L, quad2);

      /**
       * \subsubsection{Creating the linear problem object}
       * 
       * A {\tt LinearProblem} object contains everything that is needed to assemble
       * a discrete approximation to our PDE: a mesh, a weak form, boundary conditions, 
       * specification of test and unknown functions, and a specification 
       * of the low-level matrix and vector representation to be used. 
       * 
       * We can now create a problem object
       */ 

      LinearProblem prob(mesh, eqn, bc, phiHat, phi, vecType);
      /**
       * It may seem unnecessary to provide {\tt phiHat} and {\tt phi} as
       * constructor arguments here; after all, the test and unknown functions
       * could be deduced from the weak form.  In more complex problems with
       * vector-valued unknowns, however, we will want to specify the {\it
       * order} in which the different unknowns and test functions appear, and
       * we may want to group unknowns and test functions into blocks to create
       * a block linear system.  Such considerations can make a great
       * difference in the performance of linear solvers for some problems.
       * The test and unknown slots in the linear problem constructor are used
       * to pass information about the function ordering and blocking to the
       * linear problem; these features will be used to effect in subsequent
       * examples.
       *
       * \subsubsection{Specifying the solver}
       * 
       * A good choice of solver for this problem is BICGSTAB with ILU preconditioning.
       * \footnote{The Laplacian operator is symmetric positive definite,
       * however, the imposition of essential BCs destroys the symmetry making
       * conjugate gradients unsuitable. Symmetry may be restored by doing
       * block manipulations or by using Robin BCs as an 
       * approximation to Dirichlet BCs.}
       * We'll use level 1 preconditioning, and ask for a convergence 
       * tolerance of $10^-14$
       * within $1000$ iterations. The TSF BICGSTAB solver is configured using
       * a Trilinos \verb+ParameterList+ object, read in from an XML file
       */

      ParameterXMLFileReader reader("./bicgstab.xml");
      ParameterList solverParams = reader.getParameters();
      /**
       * \subsubsection{Solving the problem}
       *
       * The syntax of Sundance makes the next step look simpler than it really is:
       */

      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(solverParams);


      /* solve the problem */
      Expr soln = prob.solve(linSolver);
      /**
       * What is happening under the hood is that the problem object {\tt prob} 
       * builds a stiffness
       * matrix and load vector, feeds that matrix and vector into 
       * the linear solver {\tt solver}.
       * If all goes well, a solution vector is returned from the solver, 
       * and that solution
       * vector is captured into a discrete function wrapped in the expression object 
       * {\tt soln}. 
       * 
       * \subsubsection{Computing Velocity Fields}
       * 
       * We have computed the potential $\phi$, however, we will often want to obtain
       * the velocity field $\uvec=\nabla\phi$. The gradient of the solution
       * $\phi$ is not a continuous function and its components do not exist in the same
       * function space as used for $\phi$. However, we can project the velocity field
       * onto any vector space by solving the least-squares problem of minimizing the 
       * $L^2$ norm of the difference between the velocity field and its projection
       * onto a discrete space. In Sundance, this is accomplished as follows:
       */ 

      DiscreteSpace discreteSpace(mesh, 
                                  List(new Lagrange(1), 
                                       new Lagrange(1), 
                                       new Lagrange(1)),
                                  vecType);
      L2Projector projector(discreteSpace, grad*soln);
      Expr velocity = projector.project();
      /** 
       * We now have the velocity field as a discrete function.
       * 
       * \begin{note}
       * The derivatives of a $C^0$ function are not 
       * defined at element boundaries. Derivatives of $C^0$ functions can safely be
       * used within weak forms because values on element boundaries are not
       * required. 
       * However, 
       * in cases such as visualization where
       * a derivative of a $C^0$ discrete function is to be used in such a way
       * that its nodel values are required, it is essential to project the derivative
       * onto a discrete space. 
       * \end{note}
       * 
       * \subsubsection{Viewing the solution}
       * 
       * We next write the solution in a form suitable for viewing by VTK.
       */

      FieldWriter w = new VTKWriter("Post3d");
      w.addMesh(mesh);
      w.addField("phi", new ExprFieldWrapper(soln[0]));
      w.addField("ux", new ExprFieldWrapper(velocity[0]));
      w.addField("uy", new ExprFieldWrapper(velocity[1]));
      w.addField("uz", new ExprFieldWrapper(velocity[2]));
      w.write();
      /**
      * The resulting VTK file, \verb+Post3d.vtu+, can then be visualized with
      * any VTK viewer such as Paraview. A plot of potential and particle streamlines
      * is shown in figure~\ref{PostFlow}.
      * \begin{figure}[p]
      * \epsfig{file=postFlow.ps, width=6.0in}
      * \caption{Solution for potential flow past a post in a wind tunnel. The 
      * coloring of the planar slice near the bottom of the domain corresponds 
      * to velocity potential. The curves above are streamlines; the coloring along
      * the streamline indicates the $z$ component of velocity.}
      * \label{PostFlow}
      * \end{figure}
      */ 
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
