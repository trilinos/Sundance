// -*- c++ -*-

%{
  // System includes
#include <Python.h>

  // Sundance includes
#include "SundanceFieldWriter.hpp"
#include "SundanceVTKWriter.hpp"
#include "SundanceTriangleWriter.hpp"
#include "SundanceMatlabWriter.hpp"
#include "SundanceExprFieldWrapper.hpp"

  %}




// SWIG library includes
%include "std_string.i"
%include "std_vector.i"
%include "exception.i"


namespace SundanceStdMesh
{
  class FieldWriter
  {
  public:
    FieldWriter();
    ~FieldWriter();

    void addMesh(const Mesh& mesh);

    %extend 
    {
      void addField(const std::string& name, const SundanceCore::Expr& f)
      {
        self->addField(name, new SundanceStdFwk::Internal::ExprFieldWrapper(f));
      }
    }

    void write();

    void setUndefinedValue(const double& x);

  };

}

%rename(VTKWriter) makeVTKWriter;
%rename(MatlabWriter) makeMatlabWriter;
%rename(TriangleWriter) makeTriangleWriter;


%inline %{
  /* Create a VTK writer */
  SundanceStdMesh::FieldWriter makeVTKWriter(const std::string& filename)
  {
    return new SundanceStdMesh
      ::VTKWriter(filename);
  }
  %}

%inline %{
  /* Create a Triangle writer */
  SundanceStdMesh::FieldWriter makeTriangleWriter(const std::string& filename)
  {
    return new SundanceStdMesh
      ::TriangleWriter(filename);
  }
  %}

%inline %{
  /* Create a Matlab writer */
  SundanceStdMesh::FieldWriter makeMatlabWriter(const std::string& filename)
  {
    return new SundanceStdMesh
      ::MatlabWriter(filename);
  }
  %}



