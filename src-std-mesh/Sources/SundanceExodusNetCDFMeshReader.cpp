#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceOut.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace TSFExtended;
using namespace Teuchos;
using namespace SundanceUtils;


ExodusNetCDFMeshReader::ExodusNetCDFMeshReader(const string& fname,
                                               const MeshType& meshType,
                                               const MPIComm& comm)
  : MeshReaderBase(fname, meshType, comm)
{
  TEST_FOR_EXCEPTION(nProc() > 1, RuntimeError,
                     "ExodusNetCDFMeshReader not useable with parallel meshes");
}


Mesh ExodusNetCDFMeshReader::fillMesh() const 
{
  Mesh mesh;

  RefCountPtr<ifstream> is = openFile(filename(), "NetCDF");

  string line;
  Array<string> tokens;

  /* read the header line */
  getNextLine(*is, line, tokens, '#');

  TEST_FOR_EXCEPTION(tokens[0] != "netcdf", RuntimeError,
                     "ExodusNetCDF reader expected to find [netcdf] as first "
                     "token, found " << tokens[0]);

  /* read the list of dimensions */
  getNextLine(*is, line, tokens, '#');


  TEST_FOR_EXCEPTION(tokens[0] != "dimensions:", RuntimeError,
                     "ExodusNetCDF reader expected to find [dimension:] as first "
                     "token, found " << tokens[0]);
  
  int nElem = 0;
  int nNodes = 0;
  int nElemBlocks = 0;
  int nSideSets = 0;
  int dimension = 0 ;
  Array<int> blockSizes;
  Array<int> sideSetSizes;

  while (true)
    {
      getNextLine(*is, line, tokens, '#');

      if (tokens[0] == "variables:") break;
      string keyword = tokens[0];
      string equals = tokens[1];
      TEST_FOR_EXCEPTION(equals!="=", RuntimeError, "ExodusNetCDF reader "
                         "expected [=] as second token, found "
                         << equals);

      TEST_FOR_EXCEPTION(tokens.size() < 4, RuntimeError,
                         "ExodusNetCDF reader found a dimension line with "
                         "fewer than 4 tokens");

      int val = atoi(tokens[2]);
      if (keyword=="num_dim")
        {
          dimension = val;
        }
      else if (keyword=="num_nodes")
        {
          nNodes = val;
        }
      else if (keyword=="num_elem")
        {
          nElem = val;
        }
      else if (keyword=="num_el_blk")
        {
          nElemBlocks = val;
          blockSizes.resize(nElemBlocks);
        }
      else if (keyword=="num_side_sets")
        {
          nSideSets = val;
          sideSetSizes.resize(nSideSets);
        }
      else
        {
          for (int i=0; i<nElemBlocks; i++)
            {
              if (keyword=="num_el_in_blk" + Teuchos::toString(i+1))
                {
                  blockSizes[i] = val;
                }
              break;
            }
          for (int i=0; i<nSideSets; i++)
            {
              if (keyword=="num_sides_ss" + Teuchos::toString(i+1))
                {
                  sideSetSizes[i] = val;
                }
              break;
            }
        }
    }

  

  /* skip until we find [data:] */
  while (true)
    {
      getNextLine(*is, line, tokens, '#');

      if (tokens[0] == "data:") break;
    }

  /* read the data */
  

  Array<double> coords;
  coords.reserve(nNodes*dimension);
  Array<Array<int> > connect(nElemBlocks);
  Array<Array<int> > sideSetElems(nSideSets);
  Array<Array<int> > sideSetFacets(nSideSets);

  bool doneWithData = false;
  while(!doneWithData)
    {
      if (!getNextLine(*is, line, tokens, '#')) 
        {
          doneWithData = true;
          break;
        }

      if (tokens.size()==0) 
        {

          doneWithData=true;
          break;
        }
      
      if (tokens[0]=="coord")
        {
          cerr << "coord data" << endl;
          bool done = false;
          for (int i=1; i<tokens.size(); i++)
            {
              if (tokens[i] == "=") continue;
              if (tokens[i] == ";") 
                {
                  done = true;
                  break;
                }
              coords.append(atof(StrUtils::before(tokens[i], ",")));
            }
          while (!done)
            {
              if (!getNextLine(*is, line, tokens, '#'))
                {
                  doneWithData = true;
                  break;
                }

              for (int i=0; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  coords.append(atof(StrUtils::before(tokens[i], ",")));
                }
            }
          cerr << "coords = " << coords << endl;
          continue;
        }
      /* see if the line lists connectivity data for a block */
      for (int b=0; b<nElemBlocks; b++)
        {
          if (tokens[0] == "connect" + Teuchos::toString(b+1))
            {
              cerr << "conn data" << endl;
              connect[b].reserve(blockSizes[b]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  connect[b].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }

                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      connect[b].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              cerr << "block " << b << " connect = " << connect[b] << endl;
              continue;
            }
        }

      /* see if the line lists side set element numbers */
      for (int s=0; s<nSideSets; s++)
        {
          if (tokens[0] == "elem_ss" + Teuchos::toString(s+1))
            {
              cerr << "ss elem data" << endl;
              sideSetElems[s].reserve(sideSetSizes[s]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  sideSetElems[s].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }

                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      sideSetElems[s].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              cerr << "side set " << s << " elems = " << sideSetElems[s] << endl;
              continue;
            }
        }

      /* see if the line lists side set facet numbers */
      for (int s=0; s<nSideSets; s++)
        {
          if (tokens[0] == "side_ss" + Teuchos::toString(s+1))
            {
              cerr << "ss facet data" << endl;
              sideSetFacets[s].reserve(sideSetSizes[s]);
              bool done = false;
              for (int i=1; i<tokens.size(); i++)
                {
                  if (tokens[i] == "=") continue;
                  if (tokens[i] == ";") 
                    {
                      done = true;
                      break;
                    }
                  sideSetFacets[s].append(atoi(StrUtils::before(tokens[i], ",")));
                }
              while (!done)
                {
                  if (!getNextLine(*is, line, tokens, '#'))
                    {
                      doneWithData = true;
                      break;
                    }
                  for (int i=0; i<tokens.size(); i++)
                    {
                      if (tokens[i] == "=") continue;
                      if (tokens[i] == ";") 
                        {
                          done = true;
                          break;
                        }
                      sideSetFacets[s].append(atoi(StrUtils::before(tokens[i], ",")));
                    }
                }
              cerr << "side set " << s 
                   << " facets = " << sideSetFacets[s] << endl;
              continue;
            }
        }
    }

  TEST_FOR_EXCEPTION(true, RuntimeError, "mesh generation not done");

	return mesh;
}

