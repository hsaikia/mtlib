/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "TrackInfo.h"
#include <algorithm>
#include <queue>
#include <stack>

namespace mtlib {

	void TrackInfo::init(int n1, int n2)
	{
		numFeatures1 = n1;
		numFeatures2 = n2;
		histogramDist.resize(n1, std::vector<double>(n2, 0));
		normOverlap.resize(n1, std::vector<double>(n2, 0.0));
		normOverlapOld.resize(n1, std::vector<double>(n2, 0.0));
		volDist.resize(n1, std::vector<double>(n2, 0.0));
		cmDist.resize(n1, std::vector<double>(n2, 0.0));
	}
}