/* @HEADER@ */
/* @HEADER@ */

#include "SundanceVTKWriter.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "Teuchos_XMLObject.hpp"



using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace Teuchos;
using namespace TSFExtended;



void VTKWriter::write() const 
{
  lowLevelWrite(filename(), false);
  if (nProc() > 1 && myRank()==0) lowLevelWrite(filename(), true);
}

void VTKWriter::lowLevelWrite(const string& filename, bool isPHeader) const 
{
  string PHeader = "";
  if (isPHeader) PHeader="P";

  string f = filename;
  
  if (isPHeader) f = f + ".pvtu";
  else if (nProc() > 1) 
    {
      f = f + Teuchos::toString(myRank()) + ".vtu";
    }
  else
    {
      f = f + ".vtu";
    }
  
  SUNDANCE_OUT(verbosity() > VerbSilent, "writing VTK file " << f);

  ofstream os(f.c_str());

  XMLObject head("VTKFile");
  head.addAttribute("type", PHeader + "UnstructuredGrid");
  head.addAttribute("version", "0.1");
  
  os << head.header() << endl;

	for (int i=0; i<comments().length(); i++)
		{
			os << "<!-- " << comments()[i] << " -->" << endl;
		}

  XMLObject ug(PHeader + "UnstructuredGrid");
  os << ug.header() << endl;

  if (isPHeader)
    {
      writePoints(os, isPHeader);
      writePointData(os, isPHeader);
      writeCellData(os, isPHeader);
      for (int p=0; p<nProc(); p++)
        {
          XMLObject pc("Piece");
          string pfile = filename + Teuchos::toString(p) + ".vtu";
          pc.addAttribute("Source", pfile);
          os << pc << endl;
        }
    }
  else
    {
      XMLObject pc("Piece");
      pc.addAttribute("NumberOfPoints", Teuchos::toString(mesh().numCells(0)));
      pc.addAttribute("NumberOfCells", Teuchos::toString(mesh().numCells(mesh().spatialDim())));

      os << pc.header() << endl;

      writePoints(os, false);
      writeCells(os);
      writePointData(os, false);

      os << pc.footer() << endl;
    }

	os << ug.footer() << endl;
	os << head.footer() << endl;
}

void VTKWriter::writePoints(ostream& os, bool isPHeader) const 
{
  string PHeader = "";
  if (isPHeader) PHeader="P";
  XMLObject pts(PHeader + "Points");

  os << pts.header() << endl;

  XMLObject xml(PHeader + "DataArray");
  xml.addAttribute("NumberOfComponents", "3");
  xml.addAttribute("type", "float");
  xml.addAttribute("format", "ascii");

  os << xml.header() << endl;

  /* write the points, unless this call is for the dummy header on the root proc */
  if (!isPHeader)
    {
      int np = mesh().numCells(0);
      int dim = mesh().spatialDim();
      
      for (int i=0; i<np; i++)
        {
          const Point& x = mesh().nodePosition(i);
          
          for (int d=0; d<dim; d++)
            {
              os << x[d] << " ";
            }
          for (int d=dim; d<3; d++)
            {
              os << "0.0 ";
            }
          os << endl;
        }
    }

  os << xml.footer() << endl;

  os << pts.footer() << endl;
}


void VTKWriter::writeCells(ostream& os) const 
{
  XMLObject cells("Cells");
  os << cells.header() << endl;

  XMLObject conn("DataArray");
  conn.addAttribute("type", "vtkIdType");
  conn.addAttribute("Name", "connectivity");
  conn.addAttribute("format", "ascii");

  int dim = mesh().spatialDim();
  int nc = mesh().numCells(dim);

  os << conn.header() << endl;
  
  for (int c=0; c<nc; c++)
    {
      int nNodes = dim+1;
      
      for (int i=0; i<nNodes; i++)
        {
          os << " " << mesh().facetLID(dim,c,0,i);
        }
      os << endl;
    }
  
  os << conn.footer() << endl;


  XMLObject offsets("DataArray");
  offsets.addAttribute("type", "vtkIdType");
  offsets.addAttribute("Name", "offsets");
  offsets.addAttribute("format", "ascii");
  
  os << offsets.header() << endl;

  int count = 0;
  for (int c=0; c<nc; c++)
    {
			count += mesh().numFacets(dim, c, 0);
      os << count << endl;
    }

  os << offsets.footer() << endl;

  XMLObject types("DataArray");
  types.addAttribute("type", "unsigned char");
  types.addAttribute("Name", "types");
  types.addAttribute("format", "ascii");

  os << types.header() << endl;

  CellType cellType = mesh().cellType(dim);
  for (int c=0; c<nc; c++)
    {

			int vtkCode = 0;
			switch(cellType)
				{
				case TriangleCell:
					vtkCode = 5;
					break;
				case QuadCell:
					vtkCode = 9;
					break;
				case TetCell:
					vtkCode = 10;
					break;
				}
			os << vtkCode << endl;
    }

  os << types.footer() << endl;

  os << cells.footer() << endl;
}

void VTKWriter::writePointData(ostream& os, bool isPHeader) const 
{
  string PHeader = "";
  if (isPHeader) PHeader="P";

  XMLObject xml(PHeader + "PointData");

  if (pointVectorNames().length() > 0) xml.addAttribute("Vectors", pointVectorNames()[0]);
  if (pointScalarNames().length() > 0) xml.addAttribute("Scalars", pointScalarNames()[0]);

  os << xml.header() << endl;

  for (int i=0; i<pointScalarNames().length(); i++)
    {
      writeDataArray(os, pointScalarNames()[i], pointScalarFields()[i], isPHeader, true);
    }

  for (int i=0; i<pointVectorNames().length(); i++)
    {
      writeDataArray(os, pointVectorNames()[i], pointVectorFields()[i], isPHeader, true);
    }

  os << xml.footer() << endl;
}

void VTKWriter::writeCellData(ostream& os, bool isPHeader) const 
{
  string PHeader = "";
  if (isPHeader) PHeader="P";

  XMLObject xml(PHeader + "CellData");

  if (cellVectorNames().length() > 0) xml.addAttribute("Vectors", cellVectorNames()[0]);
  if (cellScalarNames().length() > 0) xml.addAttribute("Scalars", cellScalarNames()[0]);

  os << xml.header() << endl;

  for (int i=0; i<cellScalarNames().length(); i++)
    {
      writeDataArray(os, cellScalarNames()[i], cellScalarFields()[i], isPHeader, false);
    }

  for (int i=0; i<pointVectorNames().length(); i++)
    {
      writeDataArray(os, cellVectorNames()[i], cellVectorFields()[i], isPHeader, false);
    }

  os << xml.footer() << endl;
}


void VTKWriter::writeDataArray(ostream& os, const string& name, 
                               const RefCountPtr<FieldBase>& expr, bool isPHeader, bool isPointData) const 
{
  string PHeader = "";
  if (isPHeader) PHeader="P";

  XMLObject xml(PHeader + "DataArray");
  xml.addAttribute("type", "float");
  xml.addAttribute("Name", name);
  xml.addAttribute("format", "ascii");
  
  if (expr->numElems() > 1)
    {
      xml.addAttribute("NumberOfComponents", 
                       Teuchos::toString(expr->numElems()));
    }
  
  os << xml.header() << endl;

  /* write the point|cell data, unless this is a parallel header */
  if (!isPHeader)
    {

      if (isPointData)
        {
          int numNodes = mesh().numCells(0);
          
          for (int i=0; i<numNodes; i++)
            {
              for (int j=0; j<expr->numElems(); j++)
                {
                  os << (float) expr->getData(0, i, j) << endl;
                }
            }
        }
      else
        {
          int dim = mesh().spatialDim();
          int nc = mesh().numCells(dim);
          
          for (int c=0; c<nc; c++)
            {
              for (int j=0; j<expr->numElems(); j++)
                {
                  os << (float) expr->getData(dim, c, j) << endl;
                }
            }
        }
    }

  os << xml.footer() << endl;
}

