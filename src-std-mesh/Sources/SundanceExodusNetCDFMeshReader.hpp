/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXODUSNETCDFMESHREADER_H
#define SUNDANCE_EXODUSNETCDFMESHREADER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Array.hpp"

namespace SundanceStdMesh
{
  using namespace TSFExtended;
  using namespace Teuchos;
  using namespace SundanceUtils;
  using namespace Internal;
  /**
   * ExodusNetCDFMeshReader reads a mesh from a NetCDF dump of an Exodus file.
   * This will often be less efficient than reading from an exodus file
   * directly, but does not require any proprietary libraries.
   * 
   * Utilities to dump exodus to NetCDF are available from 
   * <A HREF="http://my.unidata.ucar.edu/content/software/netcdf/index.html"> 
   * </A>
   */
  class ExodusNetCDFMeshReader : public MeshReaderBase
  {
  public:
    /** */
    ExodusNetCDFMeshReader(const string& filename, 
                           const MeshType& meshType,
                           const MPIComm& comm = MPIComm::world());

    /** virtual dtor */
    virtual ~ExodusNetCDFMeshReader(){;}


    /** Create a mesh */
    virtual Mesh fillMesh() const ;

    /** Print a short descriptive string */
    virtual string describe() const 
    {return "ExodusNetCDFMeshReader[file=" + filename() + "]";}
      

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<MeshSourceBase> getRcp() {return rcp(this);}

  private:

    
                      
#endif  /* DOXYGEN_DEVELOPER_ONLY */   
  };
}

#endif
