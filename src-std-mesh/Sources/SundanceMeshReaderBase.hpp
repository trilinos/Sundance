/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MESHREADERBASE_H
#define SUNDANCE_MESHREADERBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY


#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "Teuchos_StrUtils.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    /**
     * MeshReaderBase is a base class for mesh sources that get a mesh
     * from a file. It provides several utilities for parsing lines
     * from mesh files. 
     */
    class MeshReaderBase : public MeshSourceBase
    {
    public:
      /** Construct with a filename */
      MeshReaderBase(const std::string& filename,
                     const MeshType& meshType,
                     const MPIComm& comm)
        : MeshSourceBase(meshType, comm), filename_(filename)
      {}

      /** */
      virtual ~MeshReaderBase(){;}

    protected:
      /** access to the filename */
      const std::string& filename() const {return filename_;}

      /** convert a string to its integer value */
      int atoi(const std::string& x) const ;

      /** convert a string to its double value */
      double atof(const std::string& x) const ;

      /** Determine whether a line is empty */
      bool isEmptyLine(const std::string& x) const ;

      /** Open a file "fname" and check for success.
       * @param fname name of the file to be opened
       * @param description a description of the file, e.g., "node file",
       * to be included in any error messages generated.  
       **/
      RefCountPtr<ifstream> openFile(const string& fname, 
                                     const string& description) const ;

      /** 
       * Read the next non-empty, non-comment line from a stream
       * @param is the stream from which to get the line
       * @param line upon return, filled in with the line that was read
       * @param tokens array of space-separated tokens in the line
       * @param comment a character indicating that everything after it
       * is a comment
       */
      bool getNextLine(istream& is, string& line,
                       Array<string>& tokens,
                       char comment) const ;
    private:
      std::string filename_;
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
