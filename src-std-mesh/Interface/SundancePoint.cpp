#include "SundancePoint.hpp"



using namespace Sundance;

using namespace Teuchos;



void Point::boundsCheck(int i) const
{
  TEST_FOR_EXCEPTION(i < 0 || i>=dim_, RuntimeError,
                     "Point::boundsCheck()");
}

namespace Sundance
{
ostream& operator<<(ostream& os, const Point& point)
{
	os << "(";
	for (int i=0; i<point.dim(); i++)
		{
			os << point[i];
			if (i<(point.dim()-1)) os << ", ";
		}
	os << ")";
	return os;
}
}

