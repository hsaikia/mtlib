/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once

#include "globals.h"
#include "MergeTree.h"

namespace mtlib {

	enum HistDiffType { L1_NORM, L2_NORM, INTERSECTION, CHI_SQUARED, EMD_HAT, EMD_HAT_GD_METRIC };

	class Histogram4d {
	public:
		static void createHistograms(
			const std::vector<double>& sf,
			const Dimensions& dims,
			const double& min4dglobal,
			const double& max4dglobal,
			const FeatureIdx& maxFeatureIdx,
			const std::vector<FeatureIdx>& belongsTo,
			std::vector<Hist>& hists,
			const int hist_buckets
		);

		static void accumulateHistograms(
			std::vector<Hist>& hists,
			const std::map<FeatureIdx, FeatureIdx>& hierarchyTree
		);

		static void DTWres(
			const std::vector<std::vector<Region> >& sequences,
			const std::vector< std::vector<Hist> >& histsAll,
			std::vector<std::vector<double> >& costMatrix,
			const int totVoxels,
			const HistDiffType typ
		);

		static double DTW(const std::vector<Hist> & hist1, const std::vector<Hist> & hist2, const int totVoxels, const HistDiffType typ);

		static void DTW_matrix(const std::vector<Hist>& arr1, const std::vector<Hist>& arr2, std::vector< std::vector<double> >& matrix, const int totVoxels, const HistDiffType typ);

		//static double histDiff(const Hist& h1, const Hist& h2, const int totVoxels); // the one used for all purposes;

		static double histDiff(const Hist& h1, const Hist& h2, const int totVoxels, const HistDiffType typ); // the one used for all purposes;

		static std::string getStringIdentifier(HistDiffType type);

	private:

		static double histDiffL1(const Hist& h1, const Hist& h2, const int totVoxels);
		static double histDiffL2SizeRatio(const Hist& h1, const Hist& h2, const int totVoxels);
		static double histDiffL2Cumul(const Hist& h1, const Hist& h2, const int totVoxels);
		static double histDiffL2(const Hist& h1, const Hist& h2, const int totVoxels);
		static double histDiffIntersection(const Hist& h1, const Hist& h2);
		static double histDiffChiSquared(const Hist& h1, const Hist& h2, const int totVoxels);
		static double histDiffEMDHat(const Hist& h1, const Hist& h2, const int totVoxels);
		static double histDiffEMDHatGdMetric(const Hist& h1, const Hist& h2, const int totVoxels);

	};

}