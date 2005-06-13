#include "SundanceIntHashSet.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_Time.hpp"

using namespace SundanceUtils;
using namespace Teuchos;

int main(int argc, void** argv)
{
  try
    {
      MPISession::init(&argc, &argv);
      
      int nRow = 40000;
      int nReps = 4;
      
      Time tSet("STL set time");
      Time tHash("hash set time");
      
      for (int nData=20; nData<360; nData+=20)
      {
        {
        tSet.start();
        Array<Set<int> > sets(nRow);


        for (int r=0; r<nRow; r++)
          {
            Set<int>& s = sets[r];
            for (int rep=0; rep<nReps; rep++)
              {
                for (int x=0; x<nData; x++)
                  {
                    int y = rand() % nData;
                    s.put(y);
                  }
              }
          }
        tSet.stop();
        }
        {
          tHash.start();
          Array<IntHashSet> sets(nRow+1);
          for (int r=0; r<nRow; r++)
            {
              sets[r].setCapacity(nData+1);
            }


          for (int r=0; r<nRow; r++)
            {
            IntHashSet& s = sets[r];
            for (int rep=0; rep<nReps; rep++)
              {
                for (int x=0; x<nData; x++)
                  {
                    int y = rand() % nData;
                    s.put(y);
                  }
              }
          }
        tHash.stop();
        }

        
        cerr << nData << "\t set=" << tSet.totalElapsedTime()
             << "\t hash=" << tHash.totalElapsedTime() 
             << "\t ratio=" 
             << tHash.totalElapsedTime()/tSet.totalElapsedTime()
             << endl;
      }

      
    }
  catch(std::exception& e)
    {
      cerr << "caught exception " << e.what() << endl;
    }
  MPISession::finalize();
}
