#include "PlayaConvergenceMonitor.hpp"


namespace Playa
{

using std::endl;
using std::setw;

ConvergenceMonitor::ConvergenceMonitor()
  : iters_(),
    records_()
{}

void ConvergenceMonitor::addRecord(int iter, const Array<double>& record) const
{
  iters_.append(iter);
  records_.append(record);
}

void ConvergenceMonitor::write(std::ostream& os) const
{
  for (int i=0; i<iters_.size(); i++)
  {
    os << setw(5) << i ;
    for (int j=0; j<records_[i].size(); j++)
    {
      os << setw(16) << records_[i][j];
    }
    os << endl;
  }
}



}
