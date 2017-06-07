/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "globals.h"

namespace mtlib {

	class TrackInfo {

	public:
		void init(int n1, int n2);
		int numFeatures1;
		int numFeatures2;
		std::vector<std::vector<double> > histogramDist;
		std::vector<std::vector<double> > volDist;
		std::vector<std::vector<double> > normOverlap;
		std::vector<std::vector<double> > normOverlapOld;
		std::vector<std::vector<double> > cmDist;
	};

}