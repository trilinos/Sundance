/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXTRUSIONMESHFILTER_H
#define SUNDANCE_EXTRUSIONMESHFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTransformationBase.hpp"

namespace SundanceStdMesh
{
  using namespace TSFExtended;
  using namespace Teuchos;
using namespace SundanceUtils;
  using namespace Internal;
  /**
   * ExtrusionMeshTransformation extrudes a 2D mesh to 3D. 
   */
  class ExtrusionMeshTransformation : public MeshTransformationBase
  {
  public:
    /** Construct a filter to extrude a 2D mesh from
     * the plane \f$z=z_0\f$ to the plane \f$z=z_1\f$ in 
     * \f$n_z\f$ steps. */
    ExtrusionMeshTransformation(double z0, double z1, int nzLevels,
                        const MeshType& meshType)
      : MeshTransformationBase(meshType),
        z0_(z0), z1_(z1), nzLevels_(nzLevels) {;}

    /** virtual dtor */
    virtual ~ExtrusionMeshTransformation() {;}

    /** Apply the filter to an input mesh, returning an output mesh */
    virtual Mesh apply(const Mesh& inputMesh) const ;

    /** Print a short descriptive string */
    virtual string describe() const 
    {return "ExtrusionMeshTransformation[z0=" + Teuchos::toString(z0_)
       + ", z1=" + Teuchos::toString(z1_)
       + ", nz=" + Teuchos::toString(nzLevels_) + "]";}
      

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<MeshTransformationBase> getRcp() {return rcp(this);}

  private:
    
    /** */
    double z0_;
    
    /** */
    double z1_;

    /** */
    int nzLevels_;
    
    

#endif  /* DOXYGEN_DEVELOPER_ONLY */   
  };
}

#endif
