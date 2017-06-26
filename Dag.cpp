/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "Dag.h"
#include <iostream>
#include <queue>
#include <fstream>
#include <iomanip>
#include "Histogram4d.h"

namespace mtlib {

	void Dag::clear()
	{
		nodes_.clear();
		edges_.clear();
		nodeMap_.clear();
		edgeMap_.clear();
		nodeData_.clear();
		edgeData_.clear();
		nodeDataMap_.clear();
		edgeDataMap_.clear();
		numTimeSteps_ = 1;
	}

	Dag::Dag()
	{
		numTimeSteps_ = 1;
	}

	void Dag::addNode(dagkey key)
	{
		if (nodeMap_.find(key) == nodeMap_.end()) {
			nodeMap_[key] = nodes_.size();
			nodes_.push_back(DagNode(key));

			numTimeSteps_ = std::max(numTimeSteps_, 1 + DagNode::getNodeTime(key));

		}
	}

	void Dag::addEdge(dagkey edgeKey, dagkey curr, dagkey other, int incOrOut)
	{
		if (edgeMap_.find(edgeKey) != edgeMap_.end()) {
			//std::cout << "Edge Key already in DAG. Exiting!\n";
			return;
		}

		edgeMap_[edgeKey] = edges_.size();
		edges_.push_back(DagEdge(edgeKey, curr, other, incOrOut));

		addNode(curr);
		addNode(other);

		nodes_[nodeMap_[curr]].addOutEdge(edgeKey);
		nodes_[nodeMap_[other]].addIncEdge(edgeKey);

	}

	void Dag::addNodesData(const std::string & label, const std::vector<double>& arr)
	{
		if (arr.size() != nodes_.size()) {
			std::cout << "Wrong number of array elements passed for DAG nodes!\n";
			std::cout << "Nodes " << nodes_.size() << " Array size " << arr.size() << "\n";
			return;
		}

		if (nodeDataMap_.find(label) == nodeDataMap_.end()) {

			nodeDataMap_[label] = nodeData_.size();

			DagNodeLabeledData newSet;
			newSet.label = label;
			newSet.w = arr;
			nodeData_.push_back(newSet);
		}
		else {
			nodeData_[nodeDataMap_[label]].w = arr;
		}
	}

	void Dag::addEdgesData(const std::string & label, const std::vector<double>& arr)
	{
		if (arr.size() != edges_.size()) {
			std::cout << "Wrong number of array elements passed for DAG edges!\n";
			std::cout << "Edges " << edges_.size() << " Array size " << arr.size() << "\n";
			return;
		}

		//std::cout << "Adding new edge data array : " << label << "\n";

		if (edgeDataMap_.find(label) == edgeDataMap_.end()) {

			edgeDataMap_[label] = edgeData_.size();

			DagEdgeLabeledData newSet;
			newSet.label = label;
			newSet.w = arr;
			edgeData_.push_back(newSet);
		}
		else {
			edgeData_[edgeDataMap_[label]].w = arr;
		}
	}

	bool Dag::isEmpty() const
	{
		return nodes_.size() == 0;
	}

	void Dag::getDAGViaNode(dagkey nodeKey, Dag & dag) const
	{
		dag.clear();

		if (nodeMap_.find(nodeKey) == nodeMap_.end()) {
			return;
		}

		std::queue<dagkey> paths[2];

		for (auto g = 0; g < 2; g++) {

			paths[g].push(nodeKey);

			std::map<const dagkey, bool> seen;

			for (const auto& key : nodeMap_) {
				seen[key.first] = false;
			}

			while (!paths[g].empty()) {
				auto key = paths[g].front();
				paths[g].pop();

				if (seen[key]) {
					continue;
				}

				seen[key] = true;

				auto i = nodeMap_.find(key)->second;

				for (auto& edgeKey : nodes_[i].E[g]) {

					auto e = edgeMap_.find(edgeKey)->second;

					auto otherKey = edges_[e].getOtherNode(key);

					dag.addEdge(edgeKey, key, otherKey, g);

					paths[g].push(edges_[e].getOtherNode(key));

				}
			}

		}
	}

	void Dag::getAllPathsNaive(std::vector<Sequence>& seqAndScores, const size_t& edgeDataIdx, const ObjectiveFunction& fun) const
	{
		seqAndScores.clear();

		for (const auto& node : nodes_) {

			Dag path;
			findShortestPathViaNode(node.key_, path, edgeDataIdx, fun);

			Sequence s;

			s.score = path.scorePath(edgeDataIdx, fun);

			for (const auto& bnode : path.nodes_) {
				s.seq.push_back(bnode.key_);
			}

			seqAndScores.push_back(s);
		}
	}

	void Dag::findAllMatchesToPathViaNode(const dagkey nodeKey, const int numCandidates, const int numMatches, 
		Dag & matchedPaths, const size_t edgeDataIdx, const ObjectiveFunction& fun, const std::map<dagkey, Region>& dagNodes) const
	{

		if (dagNodes.find(nodeKey) == dagNodes.end()) {
			std::cout << "[findAllMatchesToPathViaNode] Could not find nodeKey in dagNodes!\n";
			return;
		}

		if (numCandidates > dagNodes.size()) {
			std::cout << "Too large number of Candidates\n";
			return;
		}

		matchedPaths.clear();

		struct CandidateMatch {
			long long key;
			double histSim;
			double pathDtwCostWithSel;
			size_t dtwRank;

			static bool compareCM(const CandidateMatch& cm1, const CandidateMatch& cm2) {
				return cm1.histSim < cm2.histSim;
			}

		};

		std::vector<CandidateMatch> CMs;

		auto selHist = dagNodes.find(nodeKey)->second.hist;

		for (const auto& node : dagNodes) {

			auto sim = Histogram4d::histDiff(node.second.hist, selHist, 0 , CHI_SQUARED);
			
			CandidateMatch cm;
			cm.histSim = sim;
			cm.key = node.first;
			CMs.push_back(cm);
			
		}

		sort(CMs.begin(), CMs.end(), CandidateMatch::compareCM);

		struct DTWMatches {
			size_t idx;
			double dtwScore;

			static bool compareDM(const DTWMatches& dm1, const DTWMatches& dm2) {
				return dm1.dtwScore < dm2.dtwScore;
			}

		};

		std::vector<DTWMatches> DTWMs;

		Dag path;

		findShortestPathViaNode(nodeKey, path, edgeDataIdx, fun);

		auto selHists = path.getHistsForPath(dagNodes);

		for (auto i = 0; i < numCandidates; i++) {

			Dag candPath;

			findShortestPathViaNode(CMs[i].key, candPath, edgeDataIdx, fun);

			auto candHists = candPath.getHistsForPath(dagNodes);

			auto score = Histogram4d::DTW(selHists, candHists, 0, CHI_SQUARED);

			DTWMatches dtwm;
			dtwm.idx = i;
			dtwm.dtwScore = score;
			DTWMs.push_back(dtwm);

		}

		sort(DTWMs.begin(), DTWMs.end(), DTWMatches::compareDM);

		for (auto i = 0; i < DTWMs.size(); i++) {

			CMs[DTWMs[i].idx].dtwRank = i;
			CMs[DTWMs[i].idx].pathDtwCostWithSel = DTWMs[i].dtwScore;

		}

		//print stats
		std::ofstream f, fcsv;
		f.open("dtw_stats.txt");
		fcsv.open("dtw_stats.csv");

		f << "Node Time  " << dagNodes.find(nodeKey)->second.timestep << "\n";
		f << "Node Index " << dagNodes.find(nodeKey)->second.index << "\n\n";

		for (auto i = 0; i < numCandidates; i++) {
			
			fcsv << i << "," << CMs[i].histSim << "," << CMs[i].dtwRank << "," << CMs[i].pathDtwCostWithSel << ","
				<< dagNodes.find(CMs[i].key)->second.timestep << "," << dagNodes.find(CMs[i].key)->second.index << "\n";

			f << "Candidate Rank " << i << "\n";
			f << "DTW Rank       " << CMs[i].dtwRank << "\n";
			f << "DTW Score      " << CMs[i].pathDtwCostWithSel << "\n";
			f << "Time           " << dagNodes.find(CMs[i].key)->second.timestep << "\n";
			f << "Subtree       #" << dagNodes.find(CMs[i].key)->second.index << "\n\n";

		}
		f.close();
		fcsv.close();
	}

	void Dag::getLocallyBestPathsBothWays(Dag & greedy, const size_t edgeDataIdx, const bool greaterBetter) const
	{
		greedy.clear();

		for (const auto& n : nodes_) {
			for (auto g = 0; g < 2; g++) {

				double bestEdgeVal;
				dagkey bestEdgeKey = -1;

				if (greaterBetter) {
					bestEdgeVal = std::numeric_limits<double>::lowest();
				}
				else {
					bestEdgeVal = std::numeric_limits<double>::max();
				}

				for (const auto& e : n.E[g]) {
					auto val = std::fabs(getEdgeData(e, edgeDataIdx));
					if ((greaterBetter && (val > bestEdgeVal)) || (!greaterBetter && (val < bestEdgeVal))) {
						
						bestEdgeVal = val;
						bestEdgeKey = e;
						
					}
				}

				if (bestEdgeKey != -1) {
					greedy.addEdge(bestEdgeKey, n.key_, edges_[edgeMap_.find(bestEdgeKey)->second].getOtherNode(n.key_), g);
				}

			}
		}

		this->transferNodeEdgeDataToSubDag(greedy);

	}

	void Dag::getLocallyBestPathsOneWay(Dag & greedy, const size_t edgeDataIdx, const bool greaterBetter, bool pastToFuture) const
	{
		greedy.clear();

		for (auto p = 0; p < 2; p++) {

			int g = 1 - p;

			if (!pastToFuture) {
				g = p;
			}

			for (const auto& node : nodes_) {

				//only deal with leftover nodes i.e., nodes that did not match up before
				if (p == 1) {

					if (greedy.keyExists(node.key_) && greedy.getN(node.key_).E[1 - g].size() > 0) {
						continue;
					}
				}

				double bestEdgeVal;
				dagkey bestEdgeKey = -1;

				if (greaterBetter) {
					bestEdgeVal = std::numeric_limits<double>::lowest();
				}
				else {
					bestEdgeVal = std::numeric_limits<double>::max();
				}

				for (const auto& e : node.E[g]) {
					auto val = std::fabs(getEdgeData(e, edgeDataIdx));
					if ((greaterBetter && (val > bestEdgeVal)) || (!greaterBetter && (val < bestEdgeVal))) {

						bestEdgeVal = val;
						bestEdgeKey = e;

					}
				}

				if (bestEdgeKey != -1) {
					greedy.addEdge(bestEdgeKey, node.key_, edges_[edgeMap_.find(bestEdgeKey)->second].getOtherNode(node.key_), g);
				}

			}

		}

		this->transferNodeEdgeDataToSubDag(greedy);

	}

	std::vector<Hist> Dag::getHistsForPath(const std::map<dagkey, Region>& dagNodes) const
	{
		std::vector<Hist> ret;

		std::vector<size_t> ord;
		sortByNodeTime(ord);

		for (auto i = 0; i < ord.size(); i++) {
			ret.push_back(dagNodes.find( nodes_[ord[i]].key_)->second.hist);
		}

		return ret;
	}

	void Dag::getDPTableForMatch(const dagkey nodeKey1, const dagkey nodeKey2, const size_t edgeDataIdx, const ObjectiveFunction& fun,
		const std::map<dagkey, Region>& dagNodes, std::vector<std::vector<double>>& dp) const
	{
		if (dagNodes.find(nodeKey1) == dagNodes.end() || dagNodes.find(nodeKey2) == dagNodes.end()) {
			std::cout << "[findAllMatchesToPathViaNode] Could not find nodeKeys in dagNodes!\n";
			return;
		}

		Dag path1, path2;

		findShortestPathViaNode(nodeKey1, path1, edgeDataIdx, fun);
		findShortestPathViaNode(nodeKey2, path2, edgeDataIdx, fun);

		auto hists1 = path1.getHistsForPath(dagNodes);
		auto hists2 = path2.getHistsForPath(dagNodes);

		Histogram4d::DTW_matrix(hists1, hists2, dp, 0, CHI_SQUARED);

	}

	void Dag::filterSequences(std::vector<Sequence>& seqAndScores, std::vector<std::vector<dagkey>>& sequences, const double filterRate) const
	{

		std::sort(seqAndScores.begin(), seqAndScores.end(), [](Sequence s1, Sequence s2) {return s1.score < s2.score; });

		// filtering based on % match
		// filterRate works like this
		// 0 -> everything with 0% match or above is rejected -> outputs only 1 sequence
		// 1.00 -> everything with 100% match or above is rejected -> outputs only unique sequences
		// 0.90 -> everything with 90% match or above is rejected -> outputs cluster centers with < 90% similarity

		sequences.clear();

		for (auto& sQ : seqAndScores) {

			auto s1 = DagNode::getNodeTime(sQ.seq[0]);
			auto l1 = (long long)sQ.seq.size();

			bool canBeAdded = true;

			for (auto seqInList : sequences) {

				int numMatches = 0;

				// first we need to align the sequences by time

				auto s2 = DagNode::getNodeTime(seqInList[0]);
				auto l2 = (long long)seqInList.size();

				if (s1 + l1 <= s2 || s2 + l2 <= s1) {
					continue;
				}

				for (auto k = 0; k < std::min(l1 - std::max(s2 - s1, 0LL), l2 - std::max(s1 - s2, 0LL)); k++) {
					if (sQ.seq[k + std::max(s2 - s1, 0LL)] == seqInList[k + std::max(s1 - s2, 0LL)]) {
						numMatches++;
					}
				}

				if ((static_cast<double>(numMatches) / std::max(l1, l2)) > filterRate) {
					canBeAdded = false;
					break;
				}
			}

			if (canBeAdded) {
				sequences.push_back(sQ.seq);
			}
		}
	}

	double diff(const std::vector<double>& h1, const std::vector<double>& h2) {

		double ret = 0.0;

		for (auto i = 0; i < h1.size(); i++) {
			if (h1[i] + h2[i] == 0) {
				continue;
			}
			double dist = std::fabs(h1[i] - h2[i]);
			ret += ((dist * dist) / (h1[i] + h2[i]));
		}

		return ret;

	}

	void Dag::analyzeLambda() const
	{
		// analyze the weights and suggest a good value for lambda
		std::vector<double> e_do = edgeData_[edgeDataMap_.find(OVRLP_DIST)->second].w;
		std::vector<double> e_ds = edgeData_[edgeDataMap_.find(HIST_DIST)->second].w;
		std::vector<double> e_de;

		double l = 0;
		double r = 1;

		e_de.resize(e_do.size());

		double lambda, diff_o, diff_s;

		while (true) {

			lambda = (l + r) / 2;



			for (auto i = 0; i < e_do.size(); i++) {
				e_de[i] = lambda * e_do[i] + (1 - lambda) * e_ds[i];
			}

			diff_o = diff(e_do, e_de);
			diff_s = diff(e_ds, e_de);

			if (diff_o + 1e-6 < diff_s) {
				r = lambda;
			}
			else if (diff_o > diff_s + 1e-6) {
				l = lambda;
			}
			else {
				break;
			}
		}

		std::cout << "\nAnalyzed Lambda = " << lambda << "\n";
		std::cout << "Divergence from Overlap " << diff_o << "\n";
		std::cout << "Divergence from Signature " << diff_s << "\n";

	}

	void Dag::findShortestPathViaNode(const dagkey nodeKey, Dag& path, const size_t edgeDataIdx, const ObjectiveFunction& fun) const
	{
		// Djikstra Shortest Path - forward and back
		// we minimize sqrt(sum of squared distances of edge weights / number of edges) 

		std::cout << "Shortest Path through " << nodeKey << " called\n";

		path.clear();

		if (nodeMap_.find(nodeKey) == nodeMap_.end()) {
			std::cout << "nodeKey not found in DAG. Returning.\n";
			return;
		}

		struct TraversalInfo {
			dagkey prevNode;
			dagkey prevEdge;
			double best_cost_to_root;
			int length_to_root;
		};

		for (auto g = 0; g < 2; g++) {
			std::map<dagkey, TraversalInfo> info;

			info[nodeKey].prevNode = -1;
			info[nodeKey].prevEdge = -1;
			info[nodeKey].best_cost_to_root = 0;
			info[nodeKey].length_to_root = 0;

			std::queue<dagkey> Q;

			Q.push(nodeKey);

			dagkey bestKey = -1;
			double bestCost = std::numeric_limits<double>::max();

			while (!Q.empty()) {

				auto Key = Q.front();

				//std::cout << "Key " << Key << "\n";

				Q.pop();

				auto idx = getNIdx(Key);

				if (nodes_[idx].E[g].size() == 0) { // end reached, record idx
					if (info[Key].best_cost_to_root < bestCost) {
						bestCost = info[Key].best_cost_to_root;
						bestKey = Key;
					}
				}

				for (auto edgeKey : nodes_[idx].E[g]) {

					auto eIdx = getEIdx(edgeKey);

					auto w = getEdgeData(edgeKey, edgeDataIdx);

					auto nextNodeKey = edges_[eIdx].getOtherNode(Key);

					if (info.find(nextNodeKey) == info.end()) { // first time we see this node
						info[nextNodeKey].best_cost_to_root = w;
						info[nextNodeKey].prevNode = Key;
						info[nextNodeKey].prevEdge = edgeKey;
						info[nextNodeKey].length_to_root = 1;

						Q.push(nextNodeKey); // push into queue
					}
					else {
						//calculate the cost from ... -> Key -> nextNodeKey

						double newCost = 0;

						if (fun == SUM) {
							newCost = info[Key].best_cost_to_root + w;
						}
						else if (fun == SQUARED_AVERAGE) {
							newCost = incRootSquaredMean(info[Key].best_cost_to_root, info[Key].length_to_root, w);
						}

						//if this cost is lower, update the info for nextNodeKey
						if (newCost < info[nextNodeKey].best_cost_to_root) {
							info[nextNodeKey].best_cost_to_root = newCost;
							info[nextNodeKey].prevNode = Key;
							info[nextNodeKey].prevEdge = edgeKey;
							info[nextNodeKey].length_to_root = info[Key].length_to_root + 1;
						}
					}

				}
			}

			//std::cout << "Computed till end. Now recording best path.\n";

			// push this part of the path

			auto currNode = bestKey;

			while (true) {
				if (currNode == nodeKey) {
					break;
				}

				path.addEdge(info[currNode].prevEdge, info[currNode].prevNode, currNode, g);
				currNode = info[currNode].prevNode;

			}

		}

		//debug - check if path is consistent

		std::vector<size_t> ord;
		path.sortByNodeTime(ord);

		for (auto i = 0; i < path.nodes_.size() - 1; i++) {
			auto t1 = path.nodes_[ord[i]].getNodeTime();
			auto t2 = path.nodes_[ord[i + 1]].getNodeTime();

			if (t1 + 1 == t2) {
				//good!
			}
			else {
				std::cout << "Wrong path ordering! t1 = " << t1 << " t2 = " << t2 << "\n";
				auto l = path.nodes_.size();
				std::cout << "Path length " << l << "\n";
				std::cout << "First Node time " << path.nodes_[ord[0]].getNodeTime() << "\n";
				std::cout << "Last Node time " << path.nodes_[ord[l - 1]].getNodeTime() << "\n";
			}

		}

		transferNodeEdgeDataToSubDag(path);

		std::cout << "Path of length " << path.nodes_.size() << "\n\n";

	}

	// works only for arguments < 1,048,576
	dagkey DagNode::makeKey(long long t, long long i)
	{
		return (t << 20) | i;
	}

	long long DagNode::getNodeTime() const
	{
		return key_ >> 20;
	}

	long long DagNode::getNodeIndex() const
	{
		long long mask = 0;
		for (auto i = 0; i < 20; i++) {
			mask |= 1LL << i;
		}
		return key_ & mask;
	}

	long long DagNode::getNodeTime(dagkey key)
	{
		return key >> 20;
	}

	long long DagNode::getNodeIndex(dagkey key_)
	{
		long long mask = 0;
		for (auto i = 0; i < 20; i++) {
			mask |= 1LL << i;
		}
		return key_ & mask;
	}

	long long DagEdge::getEdgeTime() const
	{
		return key_ >> 40;
	}

	long long DagEdge::getEdgeTime(dagkey key)
	{
		return key >> 40;
	}

	// works only for arguments < 1,048,576
	dagkey DagEdge::makeKey(long long t1, long long i1, long long i2)
	{
		return (t1 << 40) | (i1 << 20) | i2;
	}

	dagkey DagEdge::makeKey(const dagkey nkey1, const dagkey nkey2)
	{
		auto t1 = DagNode::getNodeTime(nkey1);
		auto t2 = DagNode::getNodeTime(nkey2);
		
		if (t1 + 1 == t2) {
			return DagEdge::makeKey(t1, DagNode::getNodeIndex(nkey1), DagNode::getNodeIndex(nkey2));
		}
		else if (t2 + 1 == t1) {
			return DagEdge::makeKey(t2, DagNode::getNodeIndex(nkey2), DagNode::getNodeIndex(nkey1));
		}
		else {
			std::cout << "Incorrect Dag Node keys passed to create Edge!\n";
			return 0;
		}
	}

	void Dag::readDAGFromFile(const std::string & filename)
	{

		clear();

		std::ifstream f;
		f.open(filename);

		size_t numN;

		f >> numN;

		for (auto i = 0; i < numN; i++) {
			dagkey key;
			f >> key;
			addNode(key);
		}

		size_t numNData;

		f >> numNData;

		for (auto i = 0; i < numNData; i++) {
			std::string label;
			f >> label;
			std::vector<double> data(numN);
			for (auto j = 0; j < numN; j++) {
				f >> data[j];
			}
			addNodesData(label, data);
		}

		size_t numE;

		f >> numE;

		for (auto i = 0; i < numE; i++) {
			dagkey key, n1, n2;
			f >> key >> n1 >> n2;
			this->addEdge(key, n1, n2, 1);
		}

		size_t numEData;

		f >> numEData;

		for (auto i = 0; i < numEData; i++) {
			std::string label;
			f >> label;
			std::vector<double> data(numE);
			for (auto j = 0; j < numE; j++) {
				f >> data[j];
			}
			this->addEdgesData(label, data);
		}

		f.close();

	}

	void Dag::writeDAGToFile(const std::string & filename) const
	{
		std::ofstream f;
		f.open(filename);

		f << nodes_.size() << "\n\n";

		for (auto i = 0; i < nodes_.size(); i++) {
			f << nodes_[i].key_ << "\n";
		}

		f << "\n";

		f << nodeData_.size() << "\n\n";

		//std::cout << "Node data arrays " << nodeData_.size() << "\n";

		for (auto i = 0; i < nodeData_.size(); i++) {

			f << Util::space2underscore(nodeData_[i].label) << "\n";

			for (auto j = 0; j < nodes_.size(); j++) {
				f << nodeData_[i].w[j] << "\n";
			}
		}

		f << "\n";

		f << edges_.size() << "\n\n";

		for (auto i = 0; i < edges_.size(); i++) {
			f << edges_[i].key_ << " " << edges_[i].idx1_ << " " << edges_[i].idx2_ << "\n";
		}

		f << "\n";

		f << edgeData_.size() << "\n\n";

		//std::cout << "Edge data arrays " << edgeData_.size() << "\n";

		for (auto i = 0; i < edgeData_.size(); i++) {

			f << Util::space2underscore(edgeData_[i].label) << "\n";

			for (auto j = 0; j < edges_.size(); j++) {
				f << edgeData_[i].w[j] << "\n";
			}
		}

		f.close();
	}

	void Dag::filterEdges(Dag& ret, const int edgeDataIdx, const bool greaterBetter, const bool manual, double thres) const
	{
		ret.clear();

		if (!manual) {
			double maxW = std::numeric_limits<double>::lowest();
			double minW = std::numeric_limits<double>::max();

			for (const auto& e : edgeData_[edgeDataIdx].w) {
				maxW = std::max(maxW, e);
				minW = std::min(minW, e);
			}

			thres = pow(10, (log10(maxW) + log10(minW)) / 2.0);
		}

		for (auto i = 0; i < edges_.size(); i++) {
			if ((greaterBetter && edgeData_[edgeDataIdx].w[i] < thres) || (!greaterBetter && edgeData_[edgeDataIdx].w[i] > thres)) {
				continue;
			}

			ret.addEdge(edges_[i].key_, edges_[i].idx1_, edges_[i].idx2_, 1);
		}

		transferNodeEdgeDataToSubDag(ret);
	}

	void Dag::filterNodes(Dag & ret, const int nodeDataIdx, const double atleast, const double atmost) const
	{
		ret.clear();

		for (auto i = 0; i < nodes_.size(); i++) {

			if (nodeData_[nodeDataIdx].w[i] < atleast || nodeData_[nodeDataIdx].w[i] > atmost) {
				continue;
			}

			ret.addNode(nodes_[i].key_);

		}

		for (auto i = 0; i < edges_.size(); i++) {

			if (ret.nodeMap_.find(edges_[i].idx1_) != ret.nodeMap_.end() && ret.nodeMap_.find(edges_[i].idx2_) != ret.nodeMap_.end()) {

				ret.addEdge(edges_[i].key_, edges_[i].idx1_, edges_[i].idx2_, 1);

			}
		}

		transferNodeEdgeDataToSubDag(ret);
	}

	void Dag::transferNodeEdgeDataToSubDag(Dag & subdag) const
	{
		subdag.nodeData_.clear();
		subdag.nodeDataMap_.clear();
		subdag.edgeData_.clear();
		subdag.edgeDataMap_.clear();

		for (auto i = 0; i < nodeData_.size(); i++) {
			std::vector<double> data(subdag.numNodes());
			for (auto j = 0; j < subdag.numNodes(); j++) {
				data[j] = nodeData_[i].w[nodeMap_.find(subdag.nodes_[j].key_)->second];
			}
			subdag.addNodesData(nodeData_[i].label, data);
		}

		for (auto i = 0; i < edgeData_.size(); i++) {
			std::vector<double> data(subdag.numEdges());
			for (auto j = 0; j < subdag.numEdges(); j++) {
				data[j] = edgeData_[i].w[edgeMap_.find(subdag.edges_[j].key_)->second];
			}
			subdag.addEdgesData(edgeData_[i].label, data);
		}
	}

	void Dag::findNodeMaxMin(const int nodeDataIdx, double & minVal, double & maxVal) const
	{
		minVal = std::numeric_limits<double>::max();
		maxVal = std::numeric_limits<double>::lowest();

		for (const auto& w : nodeData_[nodeDataIdx].w) {
			minVal = std::min(minVal, w);
			maxVal = std::max(maxVal, w);
		}
	}

	void Dag::sortByNodeTime(std::vector<size_t>& ordering) const
	{
		ordering.clear();
		ordering.resize(nodes_.size());

		std::vector<dagkey> nodeKeys(nodes_.size());

		for (auto i = 0; i < nodes_.size(); i++) {
			nodeKeys[i] = nodes_[i].key_;
		}

		sort(nodeKeys.begin(), nodeKeys.end(), DagNode::compareByTime);

		for (auto i = 0; i < nodes_.size(); i++) {
			ordering[i] = nodeMap_.find(nodeKeys[i])->second;
		}
	}

	void Dag::getAllPathsAtTime(const int time, std::vector<std::vector<dagkey>>& sequences, const size_t edgeDataIdx, const ObjectiveFunction& fun) const
	{
		sequences.clear();

		//std::ofstream f;

		//f.open("edge_weights_of_all_sequences.txt");


		for (const auto& node : nodes_) {
			if (node.getNodeTime() != time) {
				continue;
			}

			Dag path;
			findShortestPathViaNode(node.key_, path, edgeDataIdx, fun);

			std::vector<dagkey> sequence;

			for (const auto& bnode : path.nodes_) {
				sequence.push_back(bnode.key_);
			}

			sequences.push_back(sequence);

		}

		//f.close();
		//std::cout << "Added " << sequences.size() << " Sequences.\n";
	}

	long long Dag::getMemory() const
	{
		//TODO : Need to be rewritten
		auto ret = edges_.size() * sizeof(DagEdge) + nodeMap_.size() * sizeof(nodeMap_) + edgeMap_.size() * sizeof(edgeMap_);

		for (const auto& n : nodes_) {
			ret += n.getMemory();
		}

		return ret;

	}

	//double Dag::findSD(const size_t edgeDataIdx) const
	//{
	//	double mu = 0;
	//	double M2 = 0;
	//	size_t n = 0;

	//	for (const auto& e : edges_) {
	//		n = n + 1;
	//		double delta = getEdgeData(e.key_, edgeDataIdx) - mu;
	//		mu += delta / n;
	//		M2 += delta * (getEdgeData(e.key_, edgeDataIdx) - mu);
	//	}

	//	return std::sqrt(M2 / edges_.size());
	//}

	bool Dag::keyExists(const dagkey key) const
	{
		return nodeMap_.find(key) != nodeMap_.end();
	}

	const Vertex3Df Dag::getNodeCM(const dagkey key) const
	{
		Vertex3Df ret(0, 0, 0);
		if ((nodeDataMap_.find(POS_CM_X) == nodeDataMap_.end()) ||
			(nodeDataMap_.find(POS_CM_Y) == nodeDataMap_.end()) ||
			(nodeDataMap_.find(POS_CM_Z) == nodeDataMap_.end())
			) {
			return ret;
		}

		ret.x = nodeData_[nodeDataMap_.find(POS_CM_X)->second].w[nodeMap_.find(key)->second];
		ret.y = nodeData_[nodeDataMap_.find(POS_CM_Y)->second].w[nodeMap_.find(key)->second];
		ret.z = nodeData_[nodeDataMap_.find(POS_CM_Z)->second].w[nodeMap_.find(key)->second];
		return ret;
	}

	const DagNode & Dag::getN(const dagkey id) const
	{
		if (nodeMap_.find(id) != nodeMap_.end()) {
			return nodes_[nodeMap_.find(id)->second];
		}
		else {
			std::cout << "DagNode not found!\n";
			return DagNode();
		}
	}

	const DagEdge & Dag::getE(const dagkey id) const
	{
		if (edgeMap_.find(id) != edgeMap_.end()) {
			return edges_[edgeMap_.find(id)->second];
		}
		else {
			std::cout << "DagEdge not found!\n";
			return DagEdge();
		}
	}

	const size_t Dag::getNIdx(const dagkey id) const
	{
		return nodeMap_.find(id)->second;
	}

	const size_t Dag::getEIdx(const dagkey id) const
	{
		return edgeMap_.find(id)->second;
	}

	const size_t Dag::numNodes() const
	{
		return nodes_.size();
	}

	const size_t Dag::numEdges() const
	{
		return edges_.size();
	}

	const double Dag::getEdgeData(const dagkey eKey, const size_t idx) const
	{
		return edgeData_[idx].w[edgeMap_.find(eKey)->second];
	}

	double Dag::incMean(double oldMean, long oldLen, double newVal)
	{
		return (oldMean * oldLen + newVal) / (oldLen + 1);
	}

	double Dag::combinedMean(double m1, double m2, long l1, long l2)
	{
		return (m1 * l1 + m2 * l2) / (l1 + l2);
	}

	double Dag::combinedRootSquaredMean(double m1, double m2, long l1, long l2)
	{
		return std::sqrt((l1 * m1 * m1 + l2 * m2 * m2) / (l1 + l2));
	}

	double Dag::incRootSquaredMean(double oldMean, long oldLen, double newVal)
	{
		return std::sqrt((oldLen * oldMean * oldMean + newVal * newVal) / (oldLen + 1));
	}


	//assuming the object is a straight path
	double Dag::scorePath(const size_t& edgeDataIdx, const ObjectiveFunction& f) const
	{
		std::vector<double> values;

		for (const auto& e : edges_) {
			values.push_back(getEdgeData(e.key_, edgeDataIdx));
		}

		return Util::getCost(values, f);
	}
}


