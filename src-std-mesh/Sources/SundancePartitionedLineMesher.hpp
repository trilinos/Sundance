/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_PARTITIONEDLINEMESHER_H
#define SUNDANCE_PARTITIONEDLINEMESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

namespace SundanceStdMesh
{
  using namespace TSFExtended;
  using namespace Teuchos;
using namespace SundanceUtils;
  using namespace Internal;
  /**
   * PartitionedLineMesher meshes the one-dimensional interval 
   * \f$\left[a_x, b_x\right]\f$
   * with \f$n_x\f$ elements per processor. 
   */
  class PartitionedLineMesher : public MeshSourceBase
  {
  public:
    /** 
     * Set up a mesher for the interval \f$\left[a_x, b_x\right]\f$
     * with \f$n_x\f$ elements per processor. 
     */
    PartitionedLineMesher(double ax, double bx, int nx,
                          const MeshType& meshType,
                          const MPIComm& comm = MPIComm::world())
      : 
      MeshSourceBase(meshType, comm),
      ax_(ax), bx_(bx), nx_(nx) {;}

    
    /** */
    virtual ~PartitionedLineMesher() {;}

    /** Print a short descriptive string */
    virtual string describe() const 
    {return "PartitionedLineMesher[ax=" + Teuchos::toString(ax_)
       + ", bx=" + Teuchos::toString(bx_)
       + ", nx=" + Teuchos::toString(nx_) + "]";}
      

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<MeshSourceBase> getRcp() {return rcp(this);}

  private:

    /** */
    virtual Mesh fillMesh() const ;

    /** */
    double ax_;
    /** */
    double bx_;
    /** */
    int nx_;

#endif  /* DOXYGEN_DEVELOPER_ONLY */   
  };
}

#endif
