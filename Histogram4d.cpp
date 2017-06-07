/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "Histogram4d.h"
#include <math.h>

#include "FastEMD/emd_hat.hpp"
#include "FastEMD/emd_hat_signatures_interface.hpp"

namespace mtlib {

	void Histogram4d::createHistograms(
		const std::vector<double>& sf, 
		const Dimensions & dims, 
		const double & min4dglobal, 
		const double & max4dglobal, 
		const FeatureIdx& maxFeatureIdx,
		const std::vector<FeatureIdx>& belongsTo, 
		std::vector<Hist>& hists, 
		const int hist_buckets)
	{

		hists.clear();
		hists.resize(maxFeatureIdx + 1);

		for (auto i = 0; i <= maxFeatureIdx; i++) {
			hists[i].init(hist_buckets);
		}

		for (auto i = 0; i < belongsTo.size(); i++) {
			
			Vertex3D v = ScalarField::getVertex3(i, dims);
			hists[belongsTo[i]].add_elem(min4dglobal, max4dglobal, sf[i],
				Vertex3Df(static_cast<double>(v.x) / dims.nx, static_cast<double>(v.y) / dims.ny, static_cast<double>(v.z) / dims.nz));
		}

	}

	void Histogram4d::accumulateHistograms(std::vector<Hist>& hists, const std::map<FeatureIdx, FeatureIdx>& hierarchyTree)
	{
		//make the histograms cumulative

		std::vector<Hist> histsOld = hists;
		for (auto a = 0; a < hists.size(); a++) {
			auto b = hierarchyTree.find(a)->second;
			while (b >= 0 && b < hists.size()) {
				hists[b].add_hist(histsOld[a]);
				b = hierarchyTree.find(b)->second;
			}
		}
	}

	//void Histogram4d::createHistograms(
	//	const ScalarField& sf,
	//	const double& min4dglobal,
	//	const double& max4dglobal,
	//	const SimplifiedMergeTree& smt,
	//	std::vector<Hist>& hists,
	//	const int hist_buckets
	//) {

	//	hists.clear();
	//	hists.resize(smt.subtrees.size());

	//	for (auto i = 0; i < smt.subtrees.size(); i++) {
	//		hists[i].init(hist_buckets);
	//	}

	//	for (auto i = 0; i < smt.levelSetReprIdxMap.size(); i++) {
	//		Vertex3D v = ScalarField::getVertex3(i, sf.getDims());
	//		hists[smt.levelSetReprIdxMap[i]].add_elem(min4dglobal, max4dglobal, sf.F(i), 
	//			Vertex3Df(static_cast<double>(v.x) / sf.getDims().nx, static_cast<double>(v.y) / sf.getDims().ny, static_cast<double>(v.z) / sf.getDims().nz));
	//	}

	//	//make the histograms cumulative

	//	std::vector<Hist> histsOld = hists;
	//	for (const auto& edge : smt.unAugIdxTree) {
	//		auto a = edge.first;
	//		auto b = edge.second;
	//		while (b != -1) {
	//			hists[b].add_hist(histsOld[a]);
	//			b = smt.unAugIdxTree.find(b)->second;
	//		}
	//	}
	//}

	double Histogram4d::histDiffL1(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		//L1 - norm

		long long ret = 0;

		for (auto i = 0; i < h1.bins.size(); i++) {
			long long dist = std::abs(h1.bins[i] - h2.bins[i]);
			ret += dist;
		}

		if (totVoxels < 1) {
			return ret;
		}

		return  static_cast<double>(ret) / totVoxels;
	}

	double Histogram4d::histDiffL2SizeRatio(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		double SizeRatio = static_cast<double>(std::max(h1.vol, h2.vol)) / std::min(h1.vol, h2.vol);

		return histDiffL2(h1, h2, totVoxels) * SizeRatio;
	}

	double Histogram4d::histDiffL2Cumul(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		//first create the cumul hists

		Hist h1_cumul = h1;
		Hist h2_cumul = h2;

		for (int i = (int)h1.bins.size() - 1; i >= 0; i--) {
			for (int j = i + 1; j >= 0; j--) {
				h1_cumul.bins[j] += h1.bins[i];
				h2_cumul.bins[j] += h2.bins[i];
			}
		}

		return histDiffL2(h1_cumul, h2_cumul, totVoxels);

	}

	double Histogram4d::histDiffL2(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		//L2 - norm

		long long ret = 0;

		for (auto i = 0; i < h1.bins.size(); i++) {
			long long dist = std::abs(h1.bins[i] - h2.bins[i]);
			ret += (dist * dist);
		}

		if (ret < 0) {
			std::cout << "Error! Negative L2 norm of " << ret << "! \n";
		}

		if (totVoxels < 1) {
			return ret;
		}

		return  std::sqrt(static_cast<double>(ret) / 2) / totVoxels;

	}

	double Histogram4d::histDiffIntersection(const Hist & h1, const Hist & h2)
	{
		double ret = 0.0;

		for (auto i = 0; i < h1.bins.size(); i++) {
			ret += std::min(h1.bins[i], h2.bins[i]);
		}

		return ret / std::min(h1.vol, h2.vol);

	}

	double Histogram4d::histDiffChiSquared(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		double ret = 0.0;

		for (auto i = 0; i < h1.bins.size(); i++) {
			if (h1.bins[i] + h2.bins[i] == 0) {
				continue;
			}
			double dist = std::fabs(h1.bins[i] - h2.bins[i]);
			ret += ((dist * dist) / (h1.bins[i] + h2.bins[i]));
		}

		if (totVoxels < 1) {
			return ret;
		}

		return ret / totVoxels;

	}

	double Histogram4d::histDiffEMDHat(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		std::vector<std::vector<int> > eucDist(h1.bins.size(), std::vector<int>(h2.bins.size()));

		const int THRESHOLD = 5;

		for (auto i = 0; i < h1.bins.size(); i++) {
			for (auto j = 0; j < h2.bins.size(); j++) {
				eucDist[i][j] = std::min(THRESHOLD, std::abs(i - j));
			}
		}

		if (totVoxels < 1) {
			return static_cast<double>(emd_hat<int>()(h1.bins, h2.bins, eucDist, THRESHOLD));
		}

		return static_cast<double>(emd_hat<int>()(h1.bins, h2.bins, eucDist, THRESHOLD)) / totVoxels;
	}

	double Histogram4d::histDiffEMDHatGdMetric(const Hist & h1, const Hist & h2, const int totVoxels)
	{
		std::vector<std::vector<int> > eucDist(h1.bins.size(), std::vector<int>(h2.bins.size()));

		const int THRESHOLD = 5;

		for (auto i = 0; i < h1.bins.size(); i++) {
			for (auto j = 0; j < h2.bins.size(); j++) {
				eucDist[i][j] = std::min(THRESHOLD, std::abs(i - j));
			}
		}

		if (totVoxels < 1) {
			return static_cast<double>(emd_hat_gd_metric<int>()(h1.bins, h2.bins, eucDist, THRESHOLD));
		}

		return static_cast<double>(emd_hat_gd_metric<int>()(h1.bins, h2.bins, eucDist, THRESHOLD)) / totVoxels;

	}

	//smaller cost implies better match
	void Histogram4d::DTWres(const std::vector<std::vector<Region> >& sequences, const std::vector< std::vector<Hist> >& histsAll,
		std::vector<std::vector<double> >& costMatrix, const int totVoxels, const HistDiffType typ) {

		size_t n = sequences.size();
		costMatrix.clear();
		costMatrix.resize(n, std::vector<double>(n, std::numeric_limits<double>::max()));

		for (size_t i = 0; i < n; i++) {

			costMatrix[i][i] = 0;

			for (size_t j = i + 1; j < n; j++) {

				size_t l1 = sequences[i].size();
				size_t l2 = sequences[j].size();

				std::vector<std::vector<double> > dp(l1 + 1, std::vector<double>(l2 + 1, 0));

				for (size_t i1 = 1; i1 <= l1; i1++) {
					dp[i1][0] = std::numeric_limits<double>::max();
				}
				for (size_t i2 = 1; i2 <= l2; i2++) {
					dp[0][i2] = std::numeric_limits<double>::max();
				}

				for (size_t i1 = 1; i1 <= l1; i1++) {

					for (size_t i2 = 1; i2 <= l2; i2++) {
						double cost = histDiff(
							histsAll[sequences[i][i1 - 1].timestep][sequences[i][i1 - 1].index],
							histsAll[sequences[j][i2 - 1].timestep][sequences[j][i2 - 1].index],
							totVoxels, typ
						);

						dp[i1][i2] = cost + std::min(
							dp[i1 - 1][i2 - 1], //match
							std::min(
								dp[i1 - 1][i2],	//insertion
								dp[i1][i2 - 1]	//deletion
							)
						);

					}
				}

				costMatrix[i][j] = costMatrix[j][i] = dp[l1][l2] / std::max(l1, l2);
			}
		}

	}

	double Histogram4d::DTW(const std::vector<Hist>& hist1, const std::vector<Hist>& hist2, const int totVoxels, const HistDiffType typ)
	{
		int l1 = static_cast<int>(hist1.size());
		int l2 = static_cast<int>(hist2.size());

		std::vector<std::vector<double> > dp(l1 + 1, std::vector<double>(l2 + 1, 0));

		for (size_t i1 = 1; i1 <= l1; i1++) {
			dp[i1][0] = std::numeric_limits<double>::max();
		}
		for (size_t i2 = 1; i2 <= l2; i2++) {
			dp[0][i2] = std::numeric_limits<double>::max();
		}

		for (int i1 = 1; i1 <= l1; i1++) {

			for (size_t i2 = 1; i2 <= l2; i2++) {

				//for (int i2 = std::max(1, i1 - W); i2 <= std::min(l2, i1 + W); i2++) {

				double cost = histDiff(hist1[i1 - 1], hist2[i2 - 1], totVoxels, typ);

				dp[i1][i2] = cost + std::min(
					dp[i1 - 1][i2 - 1], //match
					std::min(
						dp[i1 - 1][i2],	//insertion
						dp[i1][i2 - 1]	//deletion
					)
				);

			}
		}

		//return 1.0 - (dp[l1][l2] / std::max(l1, l2)); // more is better
		return dp[l1][l2]; // / std::max(l1, l2); // less is better
	}

	void Histogram4d::DTW_matrix(const std::vector<Hist>& hist1, const std::vector<Hist>& hist2, std::vector<std::vector<double>>& dp, const int totVoxels, const HistDiffType typ)
	{
		int l1 = (int)hist1.size();
		int l2 = (int)hist2.size();

		dp.clear();
		dp.resize(l1 + 1, std::vector<double>(l2 + 1, 0));

		for (size_t i1 = 1; i1 <= l1; i1++) {
			dp[i1][0] = std::numeric_limits<double>::max();
		}
		for (size_t i2 = 1; i2 <= l2; i2++) {
			dp[0][i2] = std::numeric_limits<double>::max();
		}

		for (int i1 = 1; i1 <= l1; i1++) {

			for (size_t i2 = 1; i2 <= l2; i2++) {

				//for (int i2 = std::max(1, i1 - W); i2 <= std::min(l2, i1 + W); i2++) {

				double cost = histDiff(hist1[i1 - 1], hist2[i2 - 1], totVoxels, typ);

				dp[i1][i2] = cost + std::min(
					dp[i1 - 1][i2 - 1], //match
					std::min(
						dp[i1 - 1][i2],	//insertion
						dp[i1][i2 - 1]	//deletion
					)
				);

			}
		}

		dp.erase(dp.begin());
		for (auto& arr : dp) {
			arr.erase(arr.begin());
		}

	}

	double Histogram4d::histDiff(const Hist & h1, const Hist & h2, const int totVoxels, HistDiffType typ)
	{
		switch (typ) {
		case L1_NORM: return histDiffL1(h1, h2, totVoxels);
		case L2_NORM: return histDiffL2(h1, h2, totVoxels);
		case CHI_SQUARED: return histDiffChiSquared(h1, h2, totVoxels);
		case INTERSECTION: return histDiffIntersection(h1, h2);
		case EMD_HAT: return histDiffEMDHat(h1, h2, totVoxels);
		case EMD_HAT_GD_METRIC: return histDiffEMDHatGdMetric(h1, h2, totVoxels);
		default: return 0.0;
		}
	}

	std::string Histogram4d::getStringIdentifier(HistDiffType type)
	{
		switch (type) {
		case L1_NORM: return "l1_norm";
		case L2_NORM: return "l2_norm";
		case CHI_SQUARED: return "chi_sq";
		case INTERSECTION: return "intersection";
		case EMD_HAT: return "emd_hat";
		case EMD_HAT_GD_METRIC: return "emd_hat_gd_metric";
		default: return "not_listed";
		}
	}

}





