/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_PARTITIONEDRECTANGLEMESHER_H
#define SUNDANCE_PARTITIONEDRECTANGLEMESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

namespace SundanceStdMesh
{
  using namespace TSFExtended;
  using namespace Teuchos;
using namespace SundanceUtils;
  using namespace Internal;
  /**
   * PartitionedRectangleMesher meshes the rectangle 
   * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
   * with \f$ n_x \otimes n_y \f$ elements per processor. The 
   * rectangle is partitioned among processors, with \f$np_x\f$
   * equal sized 
   * subdomains in the \f$x\f$ direction and \f$np_y\f$ in the \f$y\f$
   * direction.
   */
  class PartitionedRectangleMesher : public MeshSourceBase
  {
  public:
    /** 
     * Set up meshing of the rectangle 
     * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
     * with \f$ n_x \otimes n_y \f$ elements per processor. The 
     * rectangle is partitioned among processors, with \f$np_x\f$
     * equal sized 
     * subdomains in the \f$x\f$ direction and \f$np_y\f$ in the \f$y\f$
     * direction.
     */
    PartitionedRectangleMesher(double ax, double bx, int nx, int npx,
                               double ay, double by, int ny, int npy,
                               const MeshType& meshType,
                               const MPIComm& comm = MPIComm::world())
      : 
      MeshSourceBase(meshType, comm),
      ax_(ax), bx_(bx), nx_(nx), npx_(npx),
      ay_(ay), by_(by), ny_(ny), npy_(npy) {;}

    
    /** */
    virtual ~PartitionedRectangleMesher() {;}

    /** Print a short descriptive string */
    virtual string describe() const 
    {return "PartitionedRectangleMesher[ax=" + Teuchos::toString(ax_)
       + ", bx=" + Teuchos::toString(bx_)
       + ", nx=" + Teuchos::toString(nx_) +
       + ", ay=" + Teuchos::toString(ay_)
       + ", by=" + Teuchos::toString(by_)
       + ", ny=" + Teuchos::toString(ny_) + "]";}
      

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<MeshSourceBase> getRcp() {return rcp(this);}

  protected:

    /** */
    virtual Mesh fillMesh() const ;

  private:

    /** */
    double ax_;
    /** */
    double bx_;
    /** */
    int nx_;
    /** */
    int npx_;

    /** */
    double ay_;
    /** */
    double by_;
    /** */
    int ny_;
    /** */
    int npy_;

#endif  /* DOXYGEN_DEVELOPER_ONLY */   
  };
}

#endif
