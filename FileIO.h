/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once
#include <vector>
#include <map>
#include <iomanip>
#include <limits>
#include <sstream>
#include <fstream>
#include "TrackInfo.h"
#include "MergeTree.h"
#include "Dag.h"

namespace mtlib {

	class FileIO {
	public:
		static void readRegionsFromCSV(const std::string& filename, std::map<dagkey, Region>& dagNodes);
		static void readMergeTreesFromCSV(const std::string& filename, std::vector<std::map<int, int> >& unAugTrees);
		static void writeCostMatrixToDisk(const std::string& filename, const std::vector<std::vector<double> >& costMatrix);
		static void readCostMatrixFromDisk(const std::string& filename, std::vector<std::vector<double> >& costMatrix);
		static void drawMatches(const std::map<size_t, size_t>& matches, const std::string& prefix1, const std::string& prefix2, std::ofstream& file);
		static void addHeaderToStatisticsFile(std::ofstream & fs, const int hist_buckets);
		static void addToStatistics(
			std::ofstream& fs,
			const int& timestep,
			const ConciseSMT& smt,
			const std::vector < Hist >& hists,
			const double globalRange
		);
	};
}