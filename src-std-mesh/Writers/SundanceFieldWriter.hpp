/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FIELDWRITER_H
#define SUNDANCE_FIELDWRITER_H

#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"
#include "TSFHandle.hpp"

namespace SundanceStdMesh
{
  /**
   * FieldWriter is the user level object for writing fields and meshes
   * to output stream. 
   *
   * <h4> Example: </h4> Write fields u0 and w0 to a VTK file "results.vtu"
   * \code
   * FieldWriter vtkWriter = new VTKWriter("results");
   * vtkWriter.addField(u0);
   * vtkWriter.addField(w0);
   * vtkWriter.write();
   * \endcode
   *
   * <h4> Example: </h4> Write verbose mesh information to cout
   * \code
   * FieldWriter writer = new VerboseFieldWriter();
   * writer.addMesh(mesh);
   * writer.write();
   * \endcode
   */
  class FieldWriter : public TSFExtended::Handle<FieldWriterBase>
  {
  public:
    /* Boilerplate handle ctors */
    HANDLE_CTORS(FieldWriter, FieldWriterBase);

    /** add a mesh to the list of things to be written */
    void addMesh(const Mesh& mesh) const ;

    /** add a field, tagging it with the given string as a name */
    void addField(const string& name, 
                  const Handle<FieldBase>& field) ;
    

    /** write to stream */
    void write() const ;
  };
}

#endif
