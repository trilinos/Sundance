#include "TriangleMeshReader.hpp"

using namespace Sundance;
using namespace Sundance::StdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;


TriangleMeshReader::TriangleMeshReader(const string& filename,
                                       const MeshType& meshType,
                                       const MPIComm& comm)
  : MeshReaderBase(filename, meshType, comm),
    nodeFilename_(),
    elemFilename_(),
    parFilename_()
{
  string nodeFilename_ = filename_;
  string elemFilename_ = filename_;
  string parFilename_ = filename_;

  if (nProc() > 1)
    {
      nodeFilename_ = nodeFilename_ + Teuchos::toString(myRank());
      parFilename_ = parFilename_ + Teuchos::toString(myRank());
      elemFilename_ = elemFilename_ + Teuchos::toString(myRank());
    }
  nodeFilename_ = nodeFilename_ + ".node";
  elemFilename_ = elemFilename_ + ".ele";
  parFilename_ = parFilename_ + ".par";
}


Mesh TriangleMeshReader::getMesh() const 
{
  Mesh mesh;

  Array<int> ptGID;
  Array<int> ptOwner;
  Array<int> cellGID;
  Array<int> cellOwner;

  readParallelInfo(ptGID, ptOwner, elemGID, elemOwner);

  readNodes(ptGID, ptOwner);

  readElems(elemGID, elemOwner);

  mesh.assignGlobalIndices();

	return mesh;
}

void TriangleMeshReader::readParallelInfo(Array<int>& ptGID, 
                                          Array<int>& ptOwner,
                                          Array<int>& elemGID, 
                                          Array<int>& elemOwner,
                                          int& nPoints,
                                          int& nCells) const
{
  try
    {
      ptGID.resize(0);
      ptOwner.resize(0);
      elemGID.resize(0);
      elemOwner.resize(0);

      /* if we're running in parallel, read the info on processor 
       * distribution */
      if (nProc() > 1)
        {
          ifstream parStream = openFile(parFilename_, "parallel info");
     
          /* read the number of processors and the processor rank in 
           * the file. These must be consistent with the current number of
           * processors and the current rank */
          getNextLine(parStream, line, tokens, '#');
      
          TEST_FOR_EXCEPTION(tokens.length() != 2, RuntimeError,
                             "TriangleMeshReader::getMesh() expects 2 entries "
                             "on the first line of .par file. In " 
                             << parFilename_ << " I found \n[" << line << "]\n");

          int np = atoi(tokens[1]);
          int pid = atoi(tokens[0]);

          /* check consistency with the current number of
           * processors and the current rank */
      
          TEST_FOR_EXCEPTION(np != nProc(), RuntimeError,
                             "TriangleMeshReader::getMesh() found "
                             "a mismatch between the current number of processors="
                             << nProc() << 
                             "and the number of processors=" << np
                             << "in the file " << parFilename_);

          TEST_FOR_EXCEPTION(pid != myRank(), RuntimeError,
                             "TriangleMeshReader::getMesh() found "
                             "a mismatch between the current processor rank="
                             << myRank() << "and the processor rank="
                             << pid << " in the file " << parFilename_);

          /* read the number of points */
          getNextLine(parStream, line, tokens, '#');

          TEST_FOR_EXCEPTION(tokens.length() != 1, RuntimeError,
                             "TriangleMeshReader::getMesh() requires 1 entry "
                             "on the second line of .par file. Found line \n[" 
                             << line << "]\n in file " << parFilename_);
      
          nPoints = StrUtils::atoi(tokens[0]);

          /* read the global ID and the owner PID for each point */
          ptGID.resize(nPoints);
          ptOwner.resize(nPoints);

          for (int i=0; i<nPoints; i++)
            {
              getNextLine(parStream, line, tokens, '#');

              TEST_FOR_EXCEPTION(tokens.length() != 3, RuntimeError,
                                 "TriangleMeshReader::getMesh() requires 3 "
                                 "entries on each line of the point section in "
                                 "the .par file. Found line \n[" << line
                                 << "]\n in file " << parFilename_);

              ptGID[i] = StrUtils::atoi(tokens[1]);
              ptOwner[i] = StrUtils::atoi(tokens[2]);
            }


          /* Read the number of elements */

          getNextLine(parStream, line, tokens, '#');

          TEST_FOR_EXCEPTION(tokens.length() != 1, RuntimeError,
                             "TriangleMeshReader::getMesh() requires 1 entry "
                             "on the cell count line of .par file. Found line \n[" 
                             << line << "]\n in file " << parFilename_);

          nElems = StrUtils::atoi(tokens[0]);

          if (verbose())
            {
              Out::println("read nElem=" + Teuchos::toString(nElems));
            }

          /* read the global ID and the owner PID for each element */

          cellGID.resize(nElems);
          cellOwner.resize(nElems);
          for (int i=0; i<nElems; i++)
            {
              getNextLine(parStream, line, tokens, '#');

              TEST_FOR_EXCEPTION(tokens.length() != 3, RuntimeError,
                                 "TriangleMeshReader::getMesh() requires 3 "
                                 "entries on each line of the element section in "
                                 "the .par file. Found line \n[" << line
                                 << "]\n in file " << parFilename_);

              cellGID[i] = StrUtils::atoi(tokens[1]);
              cellOwner[i] = StrUtils::atoi(tokens[2]);
            }
        }

      nPoints = ptGID.length();
      nElems = cellGID.length();
    }
  catch(std::exception& e)
    {
      SUNDANCE_TRACE(e);
    }
}

void TriangleMeshReader::readNodes(const Array<int>& ptGID,
                                   const Array<int>& ptOwner) const 
{
  
  /* Open the node file so we can read in the nodes */
	
	ifstream nodeStream = openFile(nodeFile_, "node info");
	
  /* read the header line */
  getNextLine(parStream, line, tokens, '#');
  TEST_FOR_EXCEPTION(tokens.length() != 4, RuntimeError,
                     "TriangleMeshReader::getMesh() requires 4 "
                     "entries on the header line in "
                     "the .node file. Found line \n[" << line
                     << "]\n in file " << nodeFilename_);
	string headerLine = line;
  if (verbose())
    {
      Out::println("read pt header " + line);
    }

  
	if (nProc()==1)
    {
      nPoints = atoi(tokens[0]);
      ptGID.resize(nPoints);
      ptOwner.resize(nPoints);
      for (int i=0; i<nPoints; i++)
        {
          ptGID[i] = i;
          ptOwner[i] = 0;
        }
    }
  else
    {
      /* If we're running in parallel, we'd better have consistent numbers
       * of points in the .node and .par file. */
      TEST_FOR_EXCEPTION(atoi(tokens[0]) != nPoints, RuntimeError,
                         "TriangleMeshReader::getMesh() found inconsistent "
                         "numbers of points in .node file and par file. Node "
                         "file " << nodeFilename_ << " had nPoints=" 
                         << atoi(tokens[0]) << " but .par file " 
                         << parFilename_ << " had nPoints=" << nPoints);
    }

  if (verbose())
    {
      Out::println("expecting to read " + Teuchos::toString(nPoints) + " points");
    }
  
	int dimension = atoi(tokens[1]);
	int nAttributes = atoi(tokens[2]);
	int nBdryMarkers = atoi(tokens[3]);

  /* now that we know the dimension, we can build the mesh object */
	Mesh mesh = meshType.createEmptyMesh(dimension);

  /* size point-related arrays */
	Array<int> ptIndices(nPoints);
	Array<int> rawIndices(nPoints);
	Array<bool> usedPoint(nPoints);
	pointAttributes_.resize(nPoints);

	int offset=0;
	bool first = true;

  /* read all the points */
  for (int count=0; count<nPoints; count++)
    {
      getNextLine(parStream, line, tokens, '#');
      
      TEST_FOR_EXCEPTION(tokens.length() 
                         != (1 + dimension + nAttributes + nBdryMarkers),
                         RuntimeError,
                         "TriangleMeshReader::getMesh() found bad node input "
                         "line. Expected " 
                         << (1 + dimension + nAttributes + nBdryMarkers)
                         << " entries but found line \n[" << 
                         line << "]\n in file " << nodeFilename_);
      /* Triangle files can use either 0-offset or 1-offset numbering. We'll
       * inspect the first node line to decide which numbering to use. */
      if (first)
        {
          offset = atoi(tokens[0]);
          TEST_FOR_EXCEPTION(offset < 0 || offset > 1, RuntimeError,
                             "TriangleMeshReader::getMesh() expected "
                             "either 0-offset or 1-offset numbering. Found an "
                             "initial offset of " << offset << " in line \n["
                             << line << "]\n of file " << nodeFilename_);
        }
      
      /* now we can add the point to the mesh */
      rawIndices[atoi(tokens[0]) - offset] = count;
			double x = atof(tokens[1]); 
			double y = atof(tokens[2]);
      double z = 0.0;
      Point pt;

      if (dimension==3)
				{
					z = atof(tokens[3]);
          pt = Point(x,y,z);
				}
			else 
				{
          pt = Point(x,y);
					ptIndices[count] = mesh.addPoint(Point(x, y), ptOwner[count], ptGID[count]);
				}

      ptIndices[count] 
        = mesh.addVertex(ptGID[count], pt, ptOwner[count], ptLabel);

			pointAttributes_[count].resize(nAttributes);
			for (int i=0; i<nAttributes; i++)
				{
					pointAttributes_[count][i] = atof(tokens[dimension+1+i]);
				}
		}

}

void TriangleMeshReader::readElems(const Array<int>& elemGID,
                                   const Array<int>& elemOwner) const 
{
	/* Open the element file */
	
	ifstream elemStream = openFile(elemFilename_, "element info");

  getNextLine(parStream, line, tokens, '#');

  TEST_FOR_EXCEPTION(tokens.length() != 3, RuntimeError,
                     "TriangleMeshReader::getMesh() requires 3 "
                     "entries on the header line in "
                     "the .ele file. Found line \n[" << line
                     << "]\n in file " << elemFilename_);
                   
  int nElems = elemGID.length();

	if (nProc()==1)
    {
      nElems = atoi(tokens[0]);
      elemGID.resize(nElems);
      elemOwner.resize(nElems);
      for (int i=0; i<nElems; i++)
        {
          elemGID[i] = i;
          elemOwner[i] = 0;
        }
    }
  else
    {
      /* If we're running in parallel, we'd better have consistent numbers
       * of points in the .node and .par file. */
      TEST_FOR_EXCEPTION(atoi(tokens[0]) != nElems, RuntimeError,
                         "TriangleMeshReader::readElems() found inconsistent "
                         "numbers of elements in .ele file and par file. Elem "
                         "file " << elemFilename_ << " had nElems=" 
                         << atoi(tokens[0]) << " but .par file " 
                         << parFilename_ << " had nElems=" << nElems);
    }

	int ptsPerElem = atoi(tokens[1]);

  TEST_FOR_EXCEPTION(ptsPerElem != dimension+1, RuntimeError,
                     "TriangleMeshReader::readElems() found inconsistency "
                     "between number of points per element=" << ptsPerElem 
                     << " and dimension=" << dimension << ". Number of pts "
                     "per element should be dimension + 1");

	nAttributes = atoi(tokens[2]);
	cellAttributes_.resize(nElems);

  Array<int> nodes(dimension);

  for (int count=0; count<nElems; count++)
    {
      getNextLine(parStream, line, tokens, '#');
      
      TEST_FOR_EXCEPTION(tokens.length() 
                         != (1 + ptsPerElem + nAttributes),
                         RuntimeError,
                         "TriangleMeshReader::readElems() found bad elem "
                         "input line. Expected " 
                         << (1 + ptsPerElem + nAttributes)
                         << " entries but found line \n[" << 
                         line << "]\n in file " << elemFilename_);

      for (int d=0; d<dimension; d++)
        {
          nodes[d] = ptGID[rawIndices[atoi(tokens[d+1])-offset]]
        }
      
      mesh.addElement(elemGID[count], nodes, elemOwner[count], elemLabel);
			
			cellAttributes_[count].resize(nAttributes);
			for (int i=0; i<nAttributes; i++)
				{
					cellAttributes_[count][i] = atof(tokens[1+ptsPerElem+i]);
				}
    }

}
