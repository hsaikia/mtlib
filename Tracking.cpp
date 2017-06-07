/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "Tracking.h"
#include "Histogram4d.h"
#include <fstream>
#include <sstream>
#include <string>
#include <queue>

namespace mtlib {

	void get_overlap_1_with_cumul_2(
		const int& a1,
		const int& a2,
		const std::vector<std::vector<int> >& overlap,
		std::vector<std::vector<int> >& overlap_1_with_cumul_2,
		const std::map <int, std::vector<int> >& itree2
	)
	{
		if (overlap_1_with_cumul_2[a1][a2] != -1) {
			return;
		}

		int ret = overlap[a1][a2];

		if (itree2.find(a2) != itree2.end()) {
			for (auto j = 0; j < itree2.find(a2)->second.size(); j++) {
				int b2 = itree2.find(a2)->second[j];
				get_overlap_1_with_cumul_2(a1, b2, overlap, overlap_1_with_cumul_2, itree2);
				ret += overlap_1_with_cumul_2[a1][b2];
			}
		}

		overlap_1_with_cumul_2[a1][a2] = ret;

	}

	void get_overlap_cumul_1_with_2(
		const int& a1,
		const int& a2,
		const std::vector<std::vector<int> >& overlap,
		std::vector<std::vector<int> >& overlap_cumul_1_with_2,
		const std::map <int, std::vector<int> >& itree1
	)
	{
		if (overlap_cumul_1_with_2[a1][a2] != -1) {
			return;
		}

		int ret = overlap[a1][a2];

		if (itree1.find(a1) != itree1.end()) {
			for (auto i = 0; i < itree1.find(a1)->second.size(); i++) {
				int b1 = itree1.find(a1)->second[i];
				get_overlap_cumul_1_with_2(b1, a2, overlap, overlap_cumul_1_with_2, itree1);
				ret += overlap_cumul_1_with_2[b1][a2];
			}
		}

		overlap_cumul_1_with_2[a1][a2] = ret;

	}

	void getOverlap12(
		const int& a1,
		const int& a2,
		const std::vector<std::vector<int> >& overlap,
		std::vector<std::vector<int> >& overlap_1_with_cumul_2,
		std::vector<std::vector<int> >& overlap_cumul_1_with_2,
		std::vector<std::vector<int> >& cumulativeOverlap,
		const std::map <int, std::vector<int> >& itree1,
		const std::map <int, std::vector<int> >& itree2
	) {

		//std::cout << "getOverlap for " << a1 << "," << a2 << "\n";

		if (cumulativeOverlap[a1][a2] != -1) {
			return;
		}

		int ret = overlap[a1][a2];

		if (itree1.find(a1) != itree1.end()) {
			for (auto i = 0; i < itree1.find(a1)->second.size(); i++) {
				int b1 = itree1.find(a1)->second[i];
				get_overlap_cumul_1_with_2(b1, a2, overlap, overlap_cumul_1_with_2, itree1);
				ret += overlap_cumul_1_with_2[b1][a2];
			}
		}

		if (itree2.find(a2) != itree2.end()) {
			for (auto j = 0; j < itree2.find(a2)->second.size(); j++) {
				int b2 = itree2.find(a2)->second[j];
				get_overlap_1_with_cumul_2(a1, b2, overlap, overlap_1_with_cumul_2, itree2);
				ret += overlap_1_with_cumul_2[a1][b2];
			}
		}

		if (itree1.find(a1) != itree1.end() && itree2.find(a2) != itree2.end()) {
			for (auto i = 0; i < itree1.find(a1)->second.size(); i++) {
				for (auto j = 0; j < itree2.find(a2)->second.size(); j++) {

					int b1 = itree1.find(a1)->second[i];
					int b2 = itree2.find(a2)->second[j];
					getOverlap12(b1, b2, overlap, overlap_1_with_cumul_2, overlap_cumul_1_with_2, cumulativeOverlap, itree1, itree2);
					ret += cumulativeOverlap[b1][b2];

				}
			}
		}
		cumulativeOverlap[a1][a2] = ret;
		return;
	}

	void getCumulVox(
		const int& root,
		const std::vector<int>& numVoxelsForLS,
		std::vector<int>& cumulNumVoxelsForLS,
		const std::map <int, std::vector<int> >& itree
	) {

		//std::cout << "getCumulVox for " << root << "\n";

		if (cumulNumVoxelsForLS[root] != -1) {
			return;
		}

		int ret = numVoxelsForLS[root];

		if (itree.find(root) != itree.end()) {
			for (auto i = 0; i < itree.find(root)->second.size(); i++) {

				int child = itree.find(root)->second[i];

				getCumulVox(child, numVoxelsForLS, cumulNumVoxelsForLS, itree);

				ret += cumulNumVoxelsForLS[child];
			}
		}
		cumulNumVoxelsForLS[root] = ret;
		return;
	}

	void Tracking::findOverlapRecursive(
		const std::vector<FeatureIdx>& belongs1,
		const std::vector<FeatureIdx>& belongs2,
		const std::map <FeatureIdx, std::vector<FeatureIdx> >& invHierarchyTree1,
		const std::map <FeatureIdx, std::vector<FeatureIdx> >& invHierarchyTree2,
		const FeatureIdx& invTreeRoot1,
		const FeatureIdx& invTreeRoot2,
		TrackInfo& tri,
		clock_t& time_taken
	)
	{

		clock_t start = std::clock();

		std::vector<std::vector<int> > overlap(tri.numFeatures1, std::vector<int>(tri.numFeatures2, 0));
		std::vector<int> numVoxelsForLS1(tri.numFeatures1, 0); // number of voxels contained in this level set
		std::vector<int> numVoxelsForLS2(tri.numFeatures2, 0); // number of voxels contained in this level set

																   //set1 and set2 are from same domain, hence they must have same number of voxels
		for (auto i = 0; i < belongs1.size(); i++) {

			overlap[belongs1[i]][belongs2[i]]++;

			numVoxelsForLS1[belongs1[i]]++;

			numVoxelsForLS2[belongs2[i]]++;

		}

		//std::cout << "Accumulating done\n";

		/* Initially overlap only contains the overlap of two merge tree edges
		* We need to extend this overlap to overlap of subtrees
		* The overlap of any child with X is added to the overlap of all ancestors of that child with X
		* Let us assume (I) --> Intersection function
		* If M = M_1 + M_2 + ... + M_n where all M_i are mutually exclusive then
		* M (I) A = M_1 (I) A + M_2 (I) A + ... + M_n (I) A
		* i.e. the intersection of sums is the sum of the intersections.
		*/

		std::vector<std::vector<int> > overlap_1_with_cumul_2(tri.numFeatures1, std::vector<int>(tri.numFeatures2, -1));
		std::vector<std::vector<int> > overlap_cumul_1_with_2(tri.numFeatures1, std::vector<int>(tri.numFeatures2, -1));

		std::vector<std::vector<int> > cumulativeOverlap(tri.numFeatures1, std::vector<int>(tri.numFeatures2, -1));;
		std::vector<int> cumulNumVoxelsForLS1(tri.numFeatures1, -1);
		std::vector<int> cumulNumVoxelsForLS2(tri.numFeatures2, -1);

		for (auto i = 0; i < tri.numFeatures1; i++) {
			for (auto j = 0; j < tri.numFeatures2; j++) {
				getOverlap12(i, j,
					overlap, overlap_1_with_cumul_2, overlap_cumul_1_with_2,
					cumulativeOverlap, invHierarchyTree1, invHierarchyTree2);
			}
		}

		getCumulVox(invTreeRoot1, numVoxelsForLS1, cumulNumVoxelsForLS1, invHierarchyTree1);
		getCumulVox(invTreeRoot2, numVoxelsForLS2, cumulNumVoxelsForLS2, invHierarchyTree2);

		for (auto i = 0; i < tri.numFeatures1; i++) {
			for (auto j = 0; j < tri.numFeatures2; j++) {
				tri.normOverlap[i][j] = static_cast<double>(cumulativeOverlap[i][j]) / (cumulNumVoxelsForLS1[i] + cumulNumVoxelsForLS2[j] - cumulativeOverlap[i][j]);
				tri.normOverlapOld[i][j] = static_cast<double>(std::abs(cumulNumVoxelsForLS1[i] - cumulNumVoxelsForLS2[j]))
					/ std::max(cumulNumVoxelsForLS1[i], cumulNumVoxelsForLS2[j]);
			}
		}

		//std::cout << "Overlap Computed.\n";

		time_taken = std::clock() - start;

	}

	void Tracking::findOverlapFlat(
		const std::vector<FeatureIdx>& belongs1,
		const std::vector<FeatureIdx>& belongs2,
		TrackInfo& tri,
		clock_t& time_taken
	)
	{
		clock_t start = std::clock();

		std::vector<std::vector<int> > overlap(tri.numFeatures1, std::vector<int>(tri.numFeatures2, 0));
		std::vector<int> numVoxelsForLS1(tri.numFeatures1, 0); // number of voxels contained in this level set
		std::vector<int> numVoxelsForLS2(tri.numFeatures2, 0); // number of voxels contained in this level set

																 //set1 and set2 are from same domain, hence they must have same number of voxels
		for (auto i = 0; i < belongs1.size(); i++) {

			overlap[belongs1[i]][belongs2[i]]++;

			numVoxelsForLS1[belongs1[i]]++;

			numVoxelsForLS2[belongs2[i]]++;

		}

		for (auto i = 0; i < tri.numFeatures1; i++) {
			for (auto j = 0; j < tri.numFeatures2; j++) {
				tri.normOverlap[i][j] = static_cast<double>(overlap[i][j]) / (numVoxelsForLS1[i] + numVoxelsForLS2[j] - overlap[i][j]);
				tri.normOverlapOld[i][j] = static_cast<double>(std::max(numVoxelsForLS1[i] - overlap[i][j], numVoxelsForLS2[j] - overlap[i][j])) 
					/ std::max(numVoxelsForLS1[i], numVoxelsForLS2[j]);
			}
		}

		time_taken = std::clock() - start;

	}

	void Tracking::findHistDiff(
		const std::vector < Hist >& hists1,
		const std::vector < Hist >& hists2,
		const HistDiffType& type,
		const int& normFactor,
		TrackInfo& tri,
		clock_t& time_taken
	)
	{
		// computing histogram diffs

		assert(tri.numFeatures1 == hists1.size() && tri.numFeatures2 == hists2.size());

		clock_t start = std::clock();

		for (auto i = 0; i < hists1.size(); i++) {
			for (auto j = 0; j < hists2.size(); j++) {
				tri.histogramDist[i][j] = Histogram4d::histDiff(hists1[i], hists2[j], normFactor, type);
				tri.volDist[i][j] = ( static_cast<double>(hists2[j].vol - hists1[i].vol) / normFactor);
				tri.cmDist[i][j] = Vertex3Df::EucDist(hists2[j].cm, hists1[i].cm);
			}
		}

		time_taken = std::clock() - start;
	}

	void Tracking::createDAG(const std::vector<std::vector<Hist>>& histsAll, const std::vector<TrackInfo>& tracks, Dag & dag, bool ignoreLastFeature)
	{
		dag.clear();

		int tot_poss_edge = 0;

		int withLast = ignoreLastFeature ? -1 : 0;

		for (auto t = 0; t < tracks.size(); t++) {
			for (auto s1 = 0; s1 < tracks[t].numFeatures1 + withLast; s1++) {
				for (auto s2 = 0; s2 < tracks[t].numFeatures2 + withLast; s2++) {

					tot_poss_edge++;

					auto edgeKey = DagEdge::makeKey(t, s1, s2);
					auto nodeKey1 = DagNode::makeKey(t, s1);
					auto nodeKey2 = DagNode::makeKey(t + 1, s2);

					dag.addEdge(edgeKey, nodeKey1, nodeKey2, 1);

				}
			}
		}

		std::vector<double> overlap_dists(dag.numEdges());
		std::vector<double> overlap_dists_old(dag.numEdges());
		std::vector<double> hist_dists(dag.numEdges());
		std::vector<double> vol_dists(dag.numEdges());
		std::vector<double> cm_dists(dag.numEdges());
		std::vector<double> node_vols(dag.numNodes());
		std::vector<double> node_cmX(dag.numNodes());
		std::vector<double> node_cmY(dag.numNodes());
		std::vector<double> node_cmZ(dag.numNodes());

		for (auto t = 0; t < tracks.size(); t++) {
			for (auto s1 = 0; s1 < tracks[t].numFeatures1 + withLast; s1++) {
				for (auto s2 = 0; s2 < tracks[t].numFeatures2 + withLast; s2++) {

					auto edgeKey = DagEdge::makeKey(t, s1, s2);
					overlap_dists_old[dag.getEIdx(edgeKey)] = tracks[t].normOverlapOld[s1][s2];
					overlap_dists[dag.getEIdx(edgeKey)] = 1.0 - tracks[t].normOverlap[s1][s2];
					hist_dists[dag.getEIdx(edgeKey)] = tracks[t].histogramDist[s1][s2];
					vol_dists[dag.getEIdx(edgeKey)] = tracks[t].volDist[s1][s2];
					cm_dists[dag.getEIdx(edgeKey)] = tracks[t].cmDist[s1][s2];
				}
			}
		}

		for (auto t = 0; t < histsAll.size(); t++) {
			for (auto i = 0; i < histsAll[t].size() + withLast; i++) {
				auto nodeKey = DagNode::makeKey(t, i);
				node_vols[dag.getNIdx(nodeKey)] = histsAll[t][i].vol;
				node_cmX[dag.getNIdx(nodeKey)] = histsAll[t][i].cm.x;
				node_cmY[dag.getNIdx(nodeKey)] = histsAll[t][i].cm.y;
				node_cmZ[dag.getNIdx(nodeKey)] = histsAll[t][i].cm.z;
			}
		}

		dag.addEdgesData(OVRLP_DIST, overlap_dists);
		dag.addEdgesData(OVRLP_DIST_OLD, overlap_dists_old);
		dag.addEdgesData(HIST_DIST, hist_dists);
		dag.addEdgesData(VOL_DIST, vol_dists);
		dag.addEdgesData(CM_DIST, cm_dists);
		dag.addNodesData(VOLUME, node_vols);
		dag.addNodesData(POS_CM_X, node_cmX);
		dag.addNodesData(POS_CM_Y, node_cmY);
		dag.addNodesData(POS_CM_Z, node_cmZ);

		std::cout << "Dag Created. Nodes " << dag.numNodes() << ".Edges " << dag.numEdges()
			<< " (" << std::setprecision(5) << 100.0 * dag.numEdges() / tot_poss_edge << "% out of " << tot_poss_edge << ")\n";


	}

	double stdDev(const std::vector<double>& nums) {
		double mean = 0;
		for (auto i = 0; i < nums.size(); i++) {
			mean = (mean * i + nums[i]) / (i + 1);
		}
		double stdD2 = 0;
		for (auto num : nums) {
			stdD2 += pow((num - mean), 2);
		}
		return sqrt(stdD2 / nums.size());
	}

}
