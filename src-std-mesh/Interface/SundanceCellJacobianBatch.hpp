/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_CELLJACOBIANBATCH_H
#define SUNDANCE_CELLJACOBIANBATCH_H


#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
  using namespace Teuchos;
  using std::string;
  using std::ostream;

  namespace Internal
    {
      /**
       * A CellJacobianBatch is a collection of Jacobian matrices 
       * for many quadrature
       * points distributed over batch of cells. All cells 
       * must have the same dimension
       * and number of quadrature points. Affine cells have constant Jacobians,
       * so in that case the quadrature point index can be ignored.
       * See the ReferenceCellBase documentation
       * for definitions of the coordinate systems used.
       *
       * <H4> Data layout </H4>
       * All Jacobian elements for all points and cells are 
       * packed into a single vector.
       * The length of the vector is 
       * \f$ N_{cell} \times N_{quad} \times D^2, \f$ where
       * \f$D\f$ is the spatial dimension.
       * The indices into the vector cycle in the following order 
       * (from slowest to fastest)
       * <ol>
       * <li> cell number
       * <li> quadrature point number
       * <li> physical coordinate direction
       * <li> reference coordinate direction
       * </ol>
       * Thus, the jacobian values for the \f$q\f$-th quadrature point on the
       * \f$c\f$-th cell are the \f$D^2\f$ entries following 
       * the \f$(q + cN_{quad})D^2\f$-th
       * element.
       *
       */
      class CellJacobianBatch
        {
        public:
          /** empty ctor */
          CellJacobianBatch();

          /** get the spatial dimension */
          int dim() const {return dim_;}

          /** get the number of cells in the batch */
          int numCells() const {return numCells_;}

          /** get the number of quad points per cell */
          int numQuadPoints() const {return numQuad_;}

          /** resize the batch */
          void resize(int numCells, int numQuad, int dim);

          /** resize the batch, using one quadrature point per cell 
           * (appropriate for affine elements) */
          void resize(int numCells, int dim);

          /** set the Jacobian values at the q-th quadrature 
           * point on the c-th cell.
           * @param cell the index of the cell in the batch
           * @param q the index of the quadrature point
           * @param jVals jacobian values 
           * <b> (overwritten during the method call) </b>
           */
          void setJacobian(int cell, int q, Array<double>& jVals);

          /** set the Jacobian values on the c-th cell,
           * assuming cell jacobians are constant (appropriate for
           * affine cells)
           * @param cell the index of the cell in the batch
           * @param jVals jacobian values 
           * <b> (overwritten during the method call) </b>
           */
          void setJacobian(int cell, Array<double>& jVals)
            {setJacobian(cell, 0, jVals);}

          /** get the vector of determinant values */
          const Array<double>& detJ() const {return detJ_;}

          /** get the vector of jacobian values */
          const double* invJ() const
            {return &(invJ_[0]);}

          /** get the inverse jacobian elements at the given (cell, quad) pair */
          const double* invJPtr(int cell, int q) const
            {return &(invJ_[dim_*dim_*(q + numQuad_*cell)]);}

          /** get the inverse jacobian elements at the given cell */
          const double* invJPtr(int cell) const
            {return &(invJ_[dim_*dim_*cell]);}

          /** */
          static bool& verbose() {static bool rtn = false; return rtn;}

        private:
          /** write J, J^T, or J^-T to a stream */
          void viewMatrix(ostream& os, const string& header, const double* vals) const ;

        private:
          int dim_;
          int numCells_;
          int numQuad_;
          Array<double> invJ_;
          Array<double>  detJ_;

        };

    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
