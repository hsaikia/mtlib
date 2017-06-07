/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once
#include <vector>
#include <iomanip>
#include <limits>
#include <sstream>
#include <fstream>
#include <ctime>
#include "MergeTree.h"
#include "TrackInfo.h"
#include "Dag.h"
#include "FileIO.h"
#include "Histogram4d.h"

namespace mtlib {

	class Tracking {
	public:
		static void findOverlapRecursive(
			const std::vector<FeatureIdx>& belongs1,
			const std::vector<FeatureIdx>& belongs2,
			const std::map <FeatureIdx, std::vector<FeatureIdx> >& invHierarchyTree1,
			const std::map <FeatureIdx, std::vector<FeatureIdx> >& invHierarchyTree2,
			const FeatureIdx& invTreeRoot1,
			const FeatureIdx& invTreeRoot2,
			TrackInfo& tri,
			clock_t& time_taken
			
		);

		static void findOverlapFlat(
			const std::vector<FeatureIdx>& belongs1,
			const std::vector<FeatureIdx>& belongs2,
			TrackInfo& tri,
			clock_t& time_taken
		);

		static void findHistDiff(
			const std::vector < Hist >& hists1,
			const std::vector < Hist >& hists2,
			const HistDiffType& type,
			const int& normFactor,
			TrackInfo& tri,
			clock_t& time_taken
			);

		static void createDAG(
			const std::vector< std::vector< Hist > >& histsAll,
			const std::vector<TrackInfo>& tracks,
			Dag& dag,
			bool ignoreLastFeature = false
		);

	};

}