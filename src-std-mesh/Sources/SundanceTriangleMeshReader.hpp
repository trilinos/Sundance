/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TRIANGLEMESHREADER_H
#define SUNDANCE_TRIANGLEMESHREADER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"

namespace SundanceStdMesh
{
  using namespace TSFExtended;
  using namespace Teuchos;
using namespace SundanceUtils;
  using namespace Internal;
  /**
   * TriangleMeshReader reads a mesh stored in Shewchuk's Triangle format.
   * This format is documented at 
   * <A HREF="http://www-2.cs.cmu.edu/~quake/triangle.html"> 
   * the Triangle homepage. 
   * </A>
   * This reader expects to find node information in <tt>.node</tt> files
   * and element information in <tt>.ele</tt> files. The <tt> filename </tt>
   * constructor argument is the stem of the filenames, and so that 
   * a reader constructed with filename <tt>joe</tt> will look for node and
   * element data in <tt>joe.node</tt> and <tt>joe.ele</tt> respectively.
   * Node and element
   * attributes are read if present, and can be accessed with the 
   * <tt>getAttributes()</tt> method of <tt>MeshSource.</tt>
   * 
   * <h4> Parallel extensions </h4>
   * We have extended the Triangle format to deal with distributed meshes.
   * A TriangleMeshReader is constructed with an MPIComm object, and if
   * that communicator has more than one processor the mesh is assumed
   * to be split into files, one for each processor. Data
   * on mesh "joe" for the <i>nnn</i>-th processor will be found in the files
   * <ul>
   * <li> <tt>joe.node.<tt><i>nnn</i>
   * <li> <tt>joe.ele.<tt><i>nnn</i>
   * <li> <tt>joe.par.<tt><i>nnn</i>
   * </ul>
   * The <tt>.node.</tt><i>nnn</i> and <tt>.ele.<tt><i>nnn</i> files contain the
   * node and element data for the part of the mesh seen 
   * by the <i>nnn</i>-th
   * processor. The node and element 
   * numberings given in those two files are <b>local</b> indexes.
   * The <tt>.par.</tt><i>nnn</i> contains node and element 
   * ownership information for the part of the mesh seen 
   * by the <i>nnn</i>-th
   * processor. 
   *
   * <br> 
   *
   * A <tt>.par</tt> file is formatted as follows:
   * <ul>
   * <li> First line: <tt> rank numProcs </tt>
   * <li> Second line: <tt> numPoints </tt>
   * <li> Next <i> nPoints </i> lines: <tt> ptLID ptGID ptOwnerRank </tt>
   * <li> Next line: <tt> numElems </tt>
   * <li> Next <i> nElems </i> lines: <tt> elemLID elemGID elemOwnerRank </tt>
   * </ul>
   * 
   */
  class TriangleMeshReader : public MeshReaderBase
  {
  public:
    /** */
    TriangleMeshReader(const string& filename, 
                       const MeshType& meshType,
                       const MPIComm& comm = MPIComm::world());

    /** virtual dtor */
    virtual ~TriangleMeshReader(){;}


    /** Create a mesh */
    virtual Mesh getMesh() const ;

    /** Print a short descriptive string */
    virtual string describe() const 
    {return "TriangleMeshReader[file=" + filename() + "]";}
      

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RefCountPtr<MeshSourceBase> getRcp() {return rcp(this);}

  private:
    /** */
    void readParallelInfo(Array<int>& ptGID, Array<int>& ptOwner,
                          Array<int>& elemGID, Array<int>& elemOwner) const ;

    /** */
    void readNodes(const Array<int>& ptGID,
                   const Array<int>& ptOwner) const ;

    /** */
    void readElems(const Array<int>& elemGID,
                   const Array<int>& elemOwner) const ;
    

    /** */
    string nodeFilename_;

    /** */
    string elemFilename_;

    /** */
    string parFilename_;
#endif  /* DOXYGEN_DEVELOPER_ONLY */   
  };
}

#endif