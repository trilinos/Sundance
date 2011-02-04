/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_Array.hpp"

using namespace Playa;
using namespace Teuchos;
using namespace std;


namespace Playa
{

void writeTable(std::ostream& os, const Tabs& tab,
  const Array<double>& a, int cols)
{
  int rows = a.size() / cols;

  for (int i=0; i<rows; i++)
  {
    os << tab << setw(10) << i << ":";
    for (int j=0; j<cols; j++) 
      os << setw(12) << setprecision(6) << a[i*cols+j];
    os << std::endl;
  }
  int n = a.size() - rows * cols;
  if (n==0) return ;
  os << tab << setw(10) << rows << ":" ;
  for (int j=0; j<n; j++) 
    os << setw(12) << setprecision(6) << a[rows*cols+j];
  os << std::endl;
}


void writeTable(std::ostream& os, const Tabs& tab,
  const Array<int>& a, int cols)
{
  int rows = a.size() / cols;

  for (int i=0; i<rows; i++)
  {
    os << tab << setw(10) << i << ":";
    for (int j=0; j<cols; j++) 
      os << setw(10) << a[i*cols+j];
    os << std::endl;
  }
  int n = a.size() - rows * cols;
  if (n==0) return ;
  os << tab << setw(10) << rows << ":" ;
  for (int j=0; j<n; j++) 
    os << setw(10) << a[rows*cols+j];
  os << std::endl;
}


}







