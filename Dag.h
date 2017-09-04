/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#pragma once
#include "globals.h"
#include "Util.h"
#include <vector>
#include <map>
#include <memory>
#include <string>

namespace mtlib {

	typedef long long dagkey;

	const std::string OVRLP_DIST = "OVERLAP_DISTANCE";
	const std::string OVRLP_DIST_OLD = "OVERLAP_DISTANCE_OLD";
	const std::string HIST_DIST  = "HISTOGRAM_DISTANCE";
	const std::string VOL_DIST   = "VOLUME_DISTANCE";
	const std::string CM_DIST    = "CENTER_OF_MASS_DISTANCE";
	const std::string VOLUME     = "VOLUME";
	const std::string POS_CM_X   = "POSITION_CENTER_OF_MASS_X";
	const std::string POS_CM_Y   = "POSITION_CENTER_OF_MASS_Y";
	const std::string POS_CM_Z   = "POSITION_CENTER_OF_MASS_Z";
	const std::string PHEROMONES = "PHEROMONES";
	const std::string PHEROMONES_INV = "PHEROMONES_INVERTED";

	struct Sequence {
		std::vector<dagkey> seq;
		double score;
	};

	struct DagEdgeLabeledData {
		std::string label;
		std::vector<double> w;

	};

	struct DagNodeLabeledData {
		std::string label;
		std::vector<double> w;
	};

	struct DagNode {

		dagkey key_;
		std::vector<dagkey> E[2]; // E[0] -> incoming edge keys, E[1] -> outgoing

		DagNode() {

		}

		long long getMemory() const {
			return (1 + E[0].size() + E[1].size()) * sizeof(long long);
		}

		DagNode(dagkey key) {
			key_ = key;
		}

		void addIncEdge(dagkey key) {
			E[0].push_back(key);
		}

		void addOutEdge(dagkey key) {
			E[1].push_back(key);
		}

		static bool compareByTime(dagkey n1, dagkey n2) {
			return DagNode::getNodeTime(n1) < DagNode::getNodeTime(n2);
		}

		long long getNodeTime() const;

		long long getNodeIndex() const;

		static long long getNodeTime(dagkey key);

		static long long getNodeIndex(dagkey key);

		static dagkey makeKey(long long t, long long i);

	};

	struct DagEdge {

		dagkey key_;
		dagkey idx1_;
		dagkey idx2_;

		DagEdge() {

		}

		DagEdge(dagkey key, dagkey curr, dagkey other, int incOrOut) {

			key_ = key;

			if (incOrOut == 0) { // incoming
				idx1_ = other;
				idx2_ = curr;
			}
			else { // outgoing
				idx1_ = curr;
				idx2_ = other;
			}

		}

		dagkey getOtherNode(dagkey node_) const {
			return (node_ == idx1_) ? idx2_ : idx1_;
		}

		static bool compareByTime(DagEdge e1, DagEdge e2) {
			return e1.getEdgeTime() < e2.getEdgeTime();
		}

		long long getEdgeTime() const;

		static long long getEdgeTime(dagkey key);

		static dagkey makeKey(long long t1, long long i1, long long i2);

		static dagkey makeKey(const dagkey nkey1, const dagkey nkey2);

	};

	class Dag {
	public:

		Dag();

		//functions

		// modifiers

		void clear();

		void addNode(dagkey key);

		void addEdge(dagkey edgeKey, dagkey curr, dagkey other, int incOrOut);

		void addNodesData(const std::string& label, const std::vector<double>& arr);

		void addEdgesData(const std::string& label, const std::vector<double>& arr);

		void readDAGFromFile(const std::string& filename);

		// const functions

		void analyzeLambda() const;

		bool isEmpty() const;

		void getDAGViaNode(dagkey nodeKey, Dag& dag) const;

		void writeDAGToFile(const std::string& filename) const;

		void sortByNodeTime(std::vector<size_t>& ordering) const;

		void filterEdges(Dag& filteredDag, const int edgeDataIdx, const bool greaterBetter, const bool manual, const double thres) const;

		void filterNodes(Dag& filteredDag, const int nodeDataIdx, const double atleast, const double atmost) const;

		// subdag is created from 'this' DAG by adding node keys and edge keys
		// This function adds as many node data/edge data arrays as in 'this' dag 
		// and transfers the node data and edge data values to their correct positions.
		void transferNodeEdgeDataToSubDag(Dag& subdag) const;

		void findNodeMaxMin(const int nodeDataIdx, double& minVal, double& maxVal) const;

		double scorePath(const size_t& edgeDataIdx, const ObjectiveFunction& f) const;

		/* Finds best paths through all nodes using Djikstra - TopoInVis 2017 */
		void getAllPaths(Dag& trackingGraph, const size_t edgeDataIdx, const ObjectiveFunction& fun, const std::string& statsFile, double& totMem) const;

		/* Naive version of getAllPaths - calls findShortestPathViaNode from all nodes */
		void getAllPathsNaive(std::vector<Sequence>& seqAndScores, const size_t& edgeDataIdx, const ObjectiveFunction& fun, double& totMem) const;

		/* Finds the best paths through all nodes in a given timestep - calls findShortestPathViaNode from all nodes in that timestep */
		void getAllPathsAtTime(const int time, std::vector<std::vector<dagkey> >& sequences, const size_t edgeDataIdx, const ObjectiveFunction& fun, double& totMem) const;

		/* Finds the best path through a node using Djikstra - EuroVis 2017 */
		void findShortestPathViaNode(const dagkey nodeKey, Dag& path, const size_t edgeDataIdx, const ObjectiveFunction& fun, double& totMem) const;

		/* Finds all best matching paths to a path through a node using Djikstra - EuroVis 2017 */
		void findAllMatchesToPathViaNode(const dagkey nodeKey, const int numCandidates, const int numMatches, 
			Dag& matchedPaths, const size_t edgeDataIdx, const ObjectiveFunction& fun, const std::map<dagkey, Region>& dagNodes) const;

		/* Greedy method */
		void getLocallyBestPathsBothWays(Dag& greedy, const size_t edgeDataIdx, const bool greaterBetter) const;

		/* Silver et al. */
		void getLocallyBestPathsOneWay(Dag& greedy, const size_t edgeDataIdx, const bool greaterBetter, bool pastToFuture) const;

		std::vector<Hist> getHistsForPath(const std::map<dagkey, Region>& dagNodes) const;
		
		void getDPTableForMatch(const dagkey nodeKey1, const dagkey nodeKey2, const size_t edgeDataIdx, const ObjectiveFunction& fun,
			const std::map<dagkey, Region>& dagNodes, std::vector< std::vector<double> >& dp) const;

		void filterSequences(std::vector<Sequence>& seqAndScores, std::vector<std::vector<dagkey> >& sequences, const double filterRate) const;

		long long getMemory() const;

		//double findSD(const size_t edgeDataIdx) const;

		bool keyExists(const dagkey key) const;

		//getters

		const Vertex3Df getNodeCM(const dagkey key) const;

		const DagNode& getN(const dagkey id) const;
		const DagEdge& getE(const dagkey id) const;

		const size_t getNIdx(const dagkey id) const;
		const size_t getEIdx(const dagkey id) const;

		const size_t numNodes() const;
		const size_t numEdges() const;

		//const double getEdgeData(const dagkey eKey) const;
		const double getEdgeData(const dagkey eKey, const size_t idx) const;

		//static functions

		static void combineDags(Dag& dag1, const Dag& dag2);

		static double incMean(double oldMean, long oldLen, double newVal);
		static double incRootSquaredMean(double oldMean, long oldLen, double newVal);
		static double combinedMean(double m1, double m2, long l1, long l2);
		static double combinedRootSquaredMean(double m1, double m2, long l1, long l2);

		std::vector<DagNode> nodes_;
		std::vector<DagEdge> edges_;

		std::vector<DagNodeLabeledData> nodeData_;
		std::vector<DagEdgeLabeledData> edgeData_;

		std::map<std::string, size_t> nodeDataMap_;
		std::map<std::string, size_t> edgeDataMap_;

		std::map<dagkey, size_t> nodeMap_;
		std::map<dagkey, size_t> edgeMap_;

		long long numTimeSteps_;

	};

}