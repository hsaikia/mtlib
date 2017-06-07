/*********************************************************************
*  Author  : Himangshu Saikia
*  Email   : saikia@kth.se
*  Project : Merge Tree Library
*
*********************************************************************
*/

#include "MergeTree.h"
#include <algorithm>
#include <set>
#include <queue>

namespace mtlib {

	MergeTree::MergeTree(ScalarField& sf, std::vector<VertexIdx>& vertexOrder, const HierarchyType& htyp)
		:sf_(sf),
		vertexOrder_(vertexOrder),
		hierarchyType_(htyp)
	{
		sortedOrder_.clear();

		numVoxels_ = vertexOrder.size();

		sortedOrder_.resize(vertexOrder.size());
		for (auto i = 0; i < vertexOrder.size(); i++) {
			sortedOrder_[vertexOrder[i]] = i;
		}
	}

	void MergeTree::createTree() {

		levelSetReprMap_.clear();
		levelSetReprMap_.resize(numVoxels_);
		activeLevelSetReprMap_.clear();
		activeLevelSetReprMap_.resize(numVoxels_, std::numeric_limits<size_t>::max());

		brIdxMap_.clear();
		brs_.clear();
		BDTree_.clear();
		simpMap_.clear();

		std::vector<bool> seen;
		seen.resize(numVoxels_, false);

		//set global Extremum
		rootCpIdx_ = vertexOrder_[numVoxels_ - 1];


		for (auto i = 0; i < numVoxels_; i++) {

			auto idx = vertexOrder_[i];

			seen[idx] = true;

			std::vector<size_t> Ns;
			ScalarField::getNeighbors(idx, sf_.getDims(), Ns);

			std::set<size_t> components;

			for (size_t ns = 0; ns < Ns.size(); ns++) {
				if (seen[Ns[ns]]) {
					size_t component = findComponent(Ns[ns]);
					components.insert(component);
				}
			}

			std::vector<size_t> componentsV(components.begin(), components.end());

			if (componentsV.size() == 0) {
				//create New
				createNewComponent(idx, EXTREMUM);
			}
			else if (componentsV.size() == 1) {
				// merge one
				mergeWithOneComponent(idx, componentsV[0]);

				//add a last critical point - the global min(max)
				if (i == numVoxels_ - 1) {
					cpIdxMap_[idx] = (int)cps_.size();

					CriticalPoint cp;
					cp.idx = idx;
					cp.typ = GBL_EXT;
					cps_.push_back(cp);

					unAugTree_[componentsV[0]] = idx;

				}
			}
			else {
				//merge many - sort the components according to their vertex order
				std::sort(componentsV.begin(), componentsV.end(), doCompare(*this));
				mergeWithManyComponents(idx, componentsV);
			}

		}

		// check for consistency of simpMap

		int numGlobalMax = 0;

		for (const auto& sm : simpMap_) {
			if (sm.second == -1) {
				numGlobalMax++;
				simpMap_[sm.first] = sm.first;
			}
		}

		if (numGlobalMax == 1) {
			//std::cout << "Simp Map Consistent.\n";
		}
		else {
			std::cout << "ERROR! Simp Map Inconsistent!! " << numGlobalMax << " Global Maximimas found!\n";
		}

		//populating height differences

		if (hierarchyType_ == PERSISTENCE) {
			for (Branch& br : brs_) {
				br.weight = sf_.getNormalizedHeightDifference(br.brReprExtIdx, br.saddleIdx);
			}
		}
		else if (hierarchyType_ == INTEGRATED_INTENSITY_DIFFERENCE) {

			double maxWt = std::numeric_limits<double>::lowest();

			for (Branch& br : brs_) {
				br.weight -= (sf_.F(br.saddleIdx) * br.vol);

				//taking the cube root of the weight as the variance between high and low weights is proportional to the volume (N^3)

				br.weight = pow(br.weight, 1.0 / 3);

				maxWt = std::max(maxWt, br.weight);
			}

			for (Branch& br : brs_) {
				br.weight /= maxWt;
			}

		}

		//fix for a completely flat scalar field
		if (brs_.size() == 1) {
			brs_[0].weight = 1.0;
		}

		doConsistencyCheck();

		//createSimplificationOrder();

	}

	void MergeTree::createNewComponent(const size_t& idx, const CPType& typ) {

		levelSetReprMap_[idx] = idx;

		cpIdxMap_[idx] = (FeatureIdx)cps_.size();

		CriticalPoint cp;
		cp.idx = idx;
		cp.typ = typ;
		cps_.push_back(cp);

		simpMap_[(int)idx] = -1;

		if (typ == EXTREMUM) {
			Branch b;
			b.brReprExtIdx = idx;
			b.level = 0;
			b.weight = sf_.F(idx);
			b.vol = 1;
			brIdxMap_[idx] = (FeatureIdx)brs_.size();
			brs_.push_back(b);
		}

	}

	void MergeTree::mergeWithOneComponent(const size_t& idx, const size_t& belongsTo) {

		levelSetReprMap_[idx] = belongsTo;
		activeLevelSetReprMap_[idx] = belongsTo;

		if (hierarchyType_ == INTEGRATED_INTENSITY_DIFFERENCE) {
			brs_[brIdxMap_[belongsTo]].weight += sf_.F(idx);
		}

		brs_[brIdxMap_[belongsTo]].saddleIdx = idx;
		brs_[brIdxMap_[belongsTo]].vol++;

	}

	void MergeTree::setSimpMap(int idx, int master)
	{
		if (simpMap_[idx] != -1) {
			setSimpMap(simpMap_[idx], master);
		}
		else {
			simpMap_[idx] = master;
		}
	}

	void MergeTree::mergeWithManyComponents(const size_t& idx, const std::vector<size_t>& saddleOf) {

		createNewComponent(idx, SADDLE);

		size_t dominantBranchId = 0;

		double maxWeight = std::numeric_limits<double>::lowest(); // VERY IMPORTANT, using DBL_MIN resulted in a massive bug!

		for (auto& compIdx : saddleOf) {

			unAugTree_[compIdx] = idx;

			// All neighbors belonging to these components will now refer to the new merged component
			activeLevelSetReprMap_[compIdx] = idx;

			brs_[brIdxMap_[compIdx]].saddleIdx = idx;
			brs_[brIdxMap_[compIdx]].level++;


			if (maxWeight < brs_[brIdxMap_[compIdx]].weight) {
				maxWeight = brs_[brIdxMap_[compIdx]].weight;
				dominantBranchId = compIdx;
			}

		}

		brIdxMap_[idx] = brIdxMap_[dominantBranchId]; // assigning new saddle as a member of the dominant branch

		setSimpMap((int)idx, (int)dominantBranchId);

		if (hierarchyType_ == INTEGRATED_INTENSITY_DIFFERENCE) {
			brs_[brIdxMap_[idx]].weight += sf_.F(idx);
		}

		brs_[brIdxMap_[idx]].vol++;

		for (auto& compIdx : saddleOf) {
			if (compIdx == dominantBranchId)
				continue;

			setSimpMap((int)compIdx, (int)dominantBranchId);

			if (hierarchyType_ == INTEGRATED_INTENSITY_DIFFERENCE) {
				brs_[brIdxMap_[dominantBranchId]].weight += brs_[brIdxMap_[compIdx]].weight;
			}
			brs_[brIdxMap_[dominantBranchId]].vol += brs_[brIdxMap_[compIdx]].vol;
			BDTree_[brIdxMap_[compIdx]] = brIdxMap_[dominantBranchId];
		}

	}

	size_t MergeTree::findComponent(const size_t& idx) {

		if (activeLevelSetReprMap_[idx] == std::numeric_limits<size_t>::max()) {
			return idx;
		}
		else {
			activeLevelSetReprMap_[idx] = findComponent(activeLevelSetReprMap_[idx]);
			return activeLevelSetReprMap_[idx];
		}

	}

	size_t MergeTree::getTreeSize() {
		return cps_.size();
	}

	size_t MergeTree::getBDTSize() {
		return brs_.size();
	}

	void MergeTree::doConsistencyCheck()
	{
		//check whether the BDT is consistent
		bool consistent = true;
		for (auto& br : BDTree_) {
			if (brs_[br.first].weight > brs_[br.second].weight) {
				consistent = false;
				std::cout << "Inconsistency found in BDT weights\n";
				break;
			}
			if (brs_[br.first].vol > brs_[br.second].vol) {
				consistent = false;
				std::cout << "Inconsistency found in BDT volumes\n";
				break;
			}
		}

		if (consistent) {
			//std::cout << "All checks consistent. BDT is valid.\n";
		}

	}

	int MergeTree::createSimplificationOrderForPath(const int& startCpIdx, const std::map<size_t, bool>& survivedCPs, std::map<int, int>& simpMap) const {
		if (!survivedCPs.find(startCpIdx)->second) {

			// CP didn't survive, let's look into the new simpMap if it's already seen before
			// if yes, simply return this value
			if (simpMap.find(startCpIdx) != simpMap.end()) {
				return simpMap[startCpIdx];
			}

			// else let's see who the parent cp is in simpMap_ for the original tree
			int nextCpIdx = simpMap_.find(startCpIdx)->second;
			simpMap[startCpIdx] = createSimplificationOrderForPath(nextCpIdx, survivedCPs, simpMap);
		}
		else {
			simpMap[startCpIdx] = startCpIdx;
		}
		return simpMap[startCpIdx];
	}

	/* This does not prune the original tree, but creates all necessary DSs of the simplified one and populates the 'SimplifiedMergeTree' object */
	/* This method follows pruning leaves and vertex reduction */
	/* It removes all branches which have a height difference lower than the threshold and recreates the tree from the remaining branches */

	void MergeTree::createSimplifiedMergeTree(double threshold, size_t numToRemain, SimplifiedMergeTree& smt, SimplificationMethod method) const {

		smt.dims = sf_.getDims();

		std::map<size_t, bool> survivedCPs;
		std::map<size_t, size_t> unAugTree;
		std::map<size_t, size_t> cpIdxMap;

		smt.unAugIdxTree.clear();
		smt.iUnAugIdxTree.clear();
		smt.cps.clear();
		smt.levelSetReprIdxMap.clear();
		smt.levelSetReprIdxMap.resize(numVoxels_);
		smt.simplifiedScalarField.clear();
		smt.simplifiedScalarField.resize(numVoxels_);

		// Set all cps to NOT survived

		for (const CriticalPoint& cp : cps_) {
			survivedCPs[cp.idx] = false;
		}

		// Now add those critical points which belong to the branches which survived

		if (method == SIM_THRESHOLD) {

			for (const Branch& br : brs_) {
				if (br.weight > threshold) {
					survivedCPs[br.brReprExtIdx] = true;
					survivedCPs[br.saddleIdx] = true;
				}
			}
		}
		else if (method == SIM_NUMBER) {

			if (numToRemain >= brs_.size()) {
				numToRemain = brs_.size();
			}

			struct br {
				double htdiff;
				size_t level;
				size_t idx;

				static bool sortBr(br a, br b) {
					if (std::fabs(a.htdiff - b.htdiff) < 0.0000001) {
						return a.level < b.level;
					}
					return a.htdiff < b.htdiff;
				}

			};

			std::vector<br> brOrder;
			for (size_t i = 0; i < brs_.size(); i++) {
				br brnch;
				brnch.htdiff = brs_[i].weight;
				brnch.level = brs_[i].level;
				brnch.idx = i;
				brOrder.push_back(brnch);
			}
			sort(brOrder.begin(), brOrder.end(), br::sortBr);

			for (size_t i = brOrder.size() - numToRemain; i < brOrder.size(); i++) {
				survivedCPs[brs_[brOrder[i].idx].brReprExtIdx] = true;
				survivedCPs[brs_[brOrder[i].idx].saddleIdx] = true;
			}
		}

		int countCps = 0;

		for (const auto& Map : survivedCPs) {
			if (Map.second) {
				countCps++;
			}
		}

		//std::cout << countCps << " CPs survived. Changing simpMap.\n";

		// Now change the simpMap for all those cps which did not survive

		std::map<int, int> simpMap;

		for (const CriticalPoint& cp : cps_) {
			createSimplificationOrderForPath((int)cp.idx, survivedCPs, simpMap);
		}

		//std::cout << "Creating new unAugTree..\n";

		// Finally,
		// 1. populate the Simplified Un-Augmented Tree
		for (const CriticalPoint& cp : cps_) {
			if (survivedCPs[cp.idx] && cp.typ != GBL_EXT) {
				//start an arc
				size_t sadd = unAugTree_.find(cp.idx)->second;

				while (!survivedCPs[sadd]) {
					sadd = unAugTree_.find(sadd)->second;
				}

				unAugTree[cp.idx] = sadd;
			}
		}

		//std::cout << "Creating new CP list\n";

		// 2. Populate the survived CP to idx map
		for (const auto& cp : cps_) {
			if (survivedCPs[cp.idx]) {
				cpIdxMap[cp.idx] = smt.cps.size();
				smt.cps.push_back(cp);
			}
		}

		//std::cout << "Creating simpUnAugTreeIdx\n";

		// 3. populate the un-augmented tree with cp indices
		for (const auto& vidx : unAugTree) {

			// edge a --> b
			size_t a = vidx.first;
			size_t b = vidx.second;

			int idx_a = (int)cpIdxMap[a];
			int idx_b = (int)cpIdxMap[b];

			smt.unAugIdxTree[idx_a] = idx_b;  // unAugIdxTree
		}

		//std::cout << "Creating new levelSetReprMap\n";

		// 4. re-create the levelSetReprIdxMap for the simplified tree
		for (size_t i = 0; i < numVoxels_; i++) {
			// i											is the voxel index
			// levelSetReprMap_[i]							is the subtree representative for the unsimplified case
			// simpMap[levelSetReprMap_[i]]					is the subtree representative for the simplified case
			// smt.cpIdxMap[simpMap[levelSetReprMap_[i]]]	is the index of the subtree representative for the simplified case
			smt.levelSetReprIdxMap[i] = (FeatureIdx)cpIdxMap[simpMap[(FeatureIdx)levelSetReprMap_[i]]];

			//if CP was not simplified
			if (simpMap[(FeatureIdx)levelSetReprMap_[i]] == (FeatureIdx)levelSetReprMap_[i]) {
				smt.simplifiedScalarField[i] = sf_.F(i);
			}
			else {
				// assign the value of the saddle of the edge to which this edge eventually merges
				// let's assume this edge is a -> b, and voxel indices are va -> vb
				// so voxel value will be F(vb)

				auto a = smt.levelSetReprIdxMap[i];
				auto b = smt.unAugIdxTree[a];
				auto vb = smt.cps[b].idx;

				smt.simplifiedScalarField[i] = sf_.F(vb);
			}

		}

		//5. create Inverted Tree
		for (const auto& edge : smt.unAugIdxTree) {
			int e1 = edge.first;
			int e2 = edge.second;

			if (e2 == smt.unAugIdxTree.size()) {
				smt.rootIdxTree = e1;
				continue;
			}

			smt.iUnAugIdxTree[e2].push_back(e1);
		}
	}

	/* Returns true even if nodeIdx == ancestorIdx*/
	bool MergeTree::isAncestor(const int& nodeIdx, const int& ancestorIdx, const std::map<int, int>& unAugTree) {
		int currNode = nodeIdx;
		while (currNode >= 0 && currNode < unAugTree.size()) {
			if (currNode == ancestorIdx) {
				return true;
			}
			if (unAugTree.find(currNode) == unAugTree.end()) {
				return false;
			}
			currNode = unAugTree.find(currNode)->second;
		}
		return false;
	}

	//needs to pass the unAugTree as it can also be a simplified one, sorted Map is obtained from create Tree

	void MergeTree::createAugTree() {

		augTree_.clear();
		augTree_.resize(levelSetReprMap_.size());

		for (size_t i = 0; i < augTree_.size(); i++) {
			augTree_[i] = i;
		}

		//given that Map is already sorted

		for (auto it = unAugTree_.begin(); it != unAugTree_.end(); it++) {
			size_t ext = it->first;
			size_t sad = it->second;

			size_t prev = ext;

			for (size_t i = 0; i < vertexOrder_.size(); i++) {
				if (levelSetReprMap_[vertexOrder_[i]] == ext) {
					augTree_[prev] = vertexOrder_[i];
					prev = vertexOrder_[i];
				}
			}

			augTree_[prev] = sad;

		}

	}

	void MergeTree::drawUnAugTree(std::ofstream& file, const std::map<int, int>& unAugTree, const std::string& prefix) {


		file << "subgraph cluster_" << prefix << "{" << std::endl;
		file << "label=\"Tree " << prefix << " \";" << std::endl;

		std::map<int, int> connections;

		for (const auto& arc : unAugTree) {
			file << prefix << arc.first << "->" << prefix << arc.second << "[dir=none,style=\"setlinewidth(4)\"];\n";
			connections[arc.second]++;
		}

		for (const auto& arc : unAugTree) {
			if (connections.find(arc.first) == connections.end()) {
				file << prefix << arc.first << "[shape=circle,label=\"" << arc.first << "\",color=red,style=\"setlinewidth(4)\"];\n";
			}
			else {
				if (connections[arc.second] > 1) {
					file << prefix << arc.first << "[shape=circle,label=\"" << arc.first << "\",color=\"#D9D919\",style=\"setlinewidth(4)\"];\n";
				}
				else {
					file << prefix << arc.first << "[shape=circle,label=\"" << arc.first << "\",color=\"#D9D919\",style=\"setlinewidth(4)\"];\n";
					file << prefix << arc.second << "[shape=circle,label=\"" << arc.second << "\",color=blue,style=\"setlinewidth(4)\"];\n";
				}
			}
		}

		file << "}\n\n";
	}

	void MergeTree::drawUnAugTree(std::ofstream & file, const std::string & prefix)
	{
		file << "subgraph cluster_" << prefix << "{" << std::endl;
		file << "label=\"Tree " << prefix << " \";" << std::endl;

		std::map<VertexIdx, int> connections;

		for (const auto& arc : unAugTree_) {
			file << prefix << cpIdxMap_[arc.first] << "->" << prefix << cpIdxMap_[arc.second] << "[dir=none,style=\"setlinewidth(4)\"];\n";
			connections[arc.second]++;
		}

		for (const auto& arc : unAugTree_) {
			if (connections.find(arc.first) == connections.end()) {
				file << prefix << cpIdxMap_[arc.first] << "[shape=circle,label=\"" << cpIdxMap_[arc.first] << "\",color=red,style=\"setlinewidth(4)\"];\n";
			}
			else {
				file << prefix << cpIdxMap_[arc.first] << "[shape=circle,label=\"" << cpIdxMap_[arc.first] << "\",color=\"#D9D919\",style=\"setlinewidth(4)\"];\n";
			}

			if (connections[arc.second] == 1) {
				file << prefix << cpIdxMap_[arc.second] << "[shape=circle,label=\"" << cpIdxMap_[arc.second] << "\",color=blue,style=\"setlinewidth(4)\"];\n";
			}

		}

		file << "}\n\n";
	}

	void MergeTree::writeUnAugTree(const std::string & filename)
	{
		std::ofstream f;
		f.open(filename);

		for (const auto& e : unAugTree_) {
			f << e.first << "\t ->\t" << e.second << "\n";
		}

		f.close();
	}

	void MergeTree::writeSimpMap(const std::string & filename)
	{
		std::ofstream f;
		f.open(filename);

		for (const auto& s : simpMap_) {
			f << s.first << "\tsimplflies to \t" << s.second << "\n";
		}

		f.close();
	}

	void MergeTree::drawBDTree(std::ofstream & file, const std::string & prefix)
	{
		file << "subgraph cluster_" << prefix << "{" << std::endl;
		file << "label=\"Tree " << prefix << " \";" << std::endl;

		for (const auto& bdtedge : BDTree_) {
			file << prefix << bdtedge.first << "->" << prefix << bdtedge.second << "[dir=none,style=\"setlinewidth(4)\"];\n";
		}

		for (int i = 0; i < brs_.size(); i++) {
			file << prefix << i << "[shape=circle,label=\"Ext = " << cpIdxMap_[brs_[i].brReprExtIdx]
				<< "|Sad = " << cpIdxMap_[brs_[i].saddleIdx] << "|W = " << brs_[i].weight << "\",color=black,style=\"setlinewidth(4)\"];\n";
		}

		file << "}\n\n";
	}

	void MergeTree::writeStatistics(std::string filename)
	{
		std::ofstream f;
		f.open(filename);

		f << "Branch_Id, Branch_EXTREMUM, Branch_SADDLE, Branch_WEIGHT, Branch_VOLUME, Branch_LEVEL, Branch_Parent_in_BDT\n";

		int id = 0;

		for (const auto& br : brs_) {
			f << id << "," << br.brReprExtIdx << "," << br.saddleIdx << "," << br.weight << "," << br.vol << "," << br.level << ",";

			if (BDTree_.find(id) != BDTree_.end()) {
				f << BDTree_[id];
			}
			else {
				f << "-1";
			}

			f << "\n";
			id++;
		}

		f.close();
	}

	//void MergeTree::writeTree(const std::string& filename, const std::map<int, int>& unAugTree){
	//	std::ofstream file;
	//	file.open(filename);
	//
	//	file << unAugTree.size() << "\n";
	//
	//	for (const auto& edge : unAugTree){
	//		file << edge.first << "\t" << edge.second << "\n";
	//	}
	//
	//	file.close();
	//
	//}

	void SimplifiedMergeTree::getIsoLevelSets(std::vector<FeatureIdx>& belongsTo, FeatureIdx& maxFeatureIdx, const double isoValue)
	{
		std::vector<FeatureIdx> EdgeIdx;

		auto eidx = 0;
		for (const auto& edge : unAugIdxTree) {

			if (simplifiedScalarField[cps[edge.first].idx] > isoValue && simplifiedScalarField[cps[edge.second].idx] < isoValue) {
				EdgeIdx.push_back(eidx);
			}
			eidx++;
		}

		maxFeatureIdx = static_cast<int>(EdgeIdx.size());

		std::map<FeatureIdx, FeatureIdx> levelSetToCompMap;

		std::set<FeatureIdx> allChildEdgeIdx;

		auto compId = 0;

		for (const auto& e : EdgeIdx) {

			std::queue<FeatureIdx> Q;

			Q.push(e);
			
			while (!Q.empty()) {

				auto eIdx = Q.front();

				allChildEdgeIdx.insert(eIdx);
				levelSetToCompMap[eIdx] = compId;

				Q.pop();

				for (const auto& ch : iUnAugIdxTree[eIdx]) {
					Q.push(ch);
				}
			}

			compId++;

		}

		belongsTo.clear();
		belongsTo.resize(levelSetReprIdxMap.size());

		for (auto i = 0; i < simplifiedScalarField.size(); i++) {
			if (simplifiedScalarField[i] > isoValue) {
				belongsTo[i] = levelSetToCompMap[levelSetReprIdxMap[i]];
			}
			else {
				belongsTo[i] = maxFeatureIdx;
			}
		}
	}

	//void MergeTree::readTree(const std::string& filename, std::map<int, int>& unAugTree){
	//	std::ifstream file;
	//	file.open(filename);
	//
	//	unAugTree.clear();
	//
	//	int n;
	//
	//	file >> n;
	//
	//	while (n--){
	//		size_t a, b;
	//		file >> a >> b;
	//		unAugTree[a] = b;
	//	}
	//
	//	file.close();
	//}
}